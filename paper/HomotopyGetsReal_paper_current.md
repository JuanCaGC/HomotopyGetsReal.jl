# HomotopyGetsReal.jl: A Julia Package for the Numerical Decomposition of Real Algebraic Curves and Surfaces via Homotopy Continuation

**Author:** Juan Camilo Gonzalez
**Address:** Courant Institute of Mathematical Sciences, New York University, New York, NY, USA
**Email:** juanca2001gon@gmail.com *(TODO: confirm — consider an NYU-affiliated address before submission)*

**2020 MSC:** 14Q65, 14Q30, 68-04
**Keywords:** numerical algebraic geometry, homotopy continuation, real algebraic curves and surfaces, cell decomposition, isosingular deflation, Julia software

---

## Abstract

HomotopyGetsReal.jl is a Julia library for the numerical decomposition of real algebraic curves and surfaces into cell complexes via homotopy continuation, built natively on HomotopyContinuation.jl. It reimplements Bertini_real's critical-point-based decomposition pipeline and extends it in two respects: an optional generic-projection mode resolving configurations where fixed-coordinate slicing fails outright, such as a torus sliced along its hole axis, and a diagnostic isosingular deflation subsystem validated exactly against published examples including the Whitney umbrella. The architecture centers on type-stable parametric data structures and a single configuration object holding every numerical tolerance used across the pipeline. We validate the implementation on a sphere, an asymmetric ellipsoid, the Taubin heart, and a genus-one torus, reporting pointwise residual statistics and the results of a 537-assertion test suite, and state explicitly that deflation remains diagnostic only: it does not yet decompose a singular curve or alter mesh geometry, a real scope gap relative to Bertini_real.

*(149 words, self-contained, no references — verified.)*

---

## 1. Introduction

Determining and visualizing the real points of a complex algebraic curve or surface is a recurring computational task in real algebraic geometry and in the engineering fields — kinematics and mechanism design in particular — that motivate much of the field's applied work. A defining polynomial system is easy to write down and easy to solve numerically over ℂ, but the real points of its solution set carry no distinguished algebraic structure of their own: they must be extracted from the complex solution set after the fact, and their topology (how many connected components there are, how branches meet at singular points, how a surface's sheets connect across critical slices) is precisely what symbolic elimination methods struggle to compute efficiently and what purely numerical point-sampling methods cannot certify. The classical exact approach, cylindrical algebraic decomposition [Arnon 1984], resolves this in principle but at doubly-exponential worst-case cost in the number of variables, which limits it in practice to systems of low complexity.

Numerical algebraic geometry offers an alternative that trades symbolic completeness for numerical scalability: homotopy continuation [Sommese–Wampler 2005; Bates et al. 2013] is used not to certify a decomposition symbolically but to track the real points of the variety directly, slicing along a projection at its critical values and interpolating between them by path tracking. This is the approach taken by Bertini_real [Brake et al. 2014, 2017], which implements the curve algorithm of Lu, Bates, Sommese, and Wampler [2007] and the surface algorithm of Besana, Di Rocco, Hauenstein, Sommese, and Wampler [2013] on top of the Bertini numerical algebraic geometry engine, and which has been used to decompose the solution curves and workspaces of planar linkage mechanisms among other examples.

Bertini_real's dependency footprint, however, is substantial: it is an MPI-parallel program requiring a library-compiled build of Bertini, Python (sympy) for symbolic computation, and the Boost, GMP, and MPFR libraries (§6). This complexity is not incidental — it reflects real engineering requirements, in particular Python/sympy's role in isosingular deflation — but it puts a working installation out of easy reach for exploratory use, and it predates the emergence of HomotopyContinuation.jl [Breiding–Timme 2018], a homotopy continuation engine written natively in Julia that has since become a standard tool in the numerical algebraic geometry community. There is, to our knowledge, no implementation of the Bertini_real-style decomposition pipeline built directly on this engine, based on a search of the Julia General registry and of the list of related packages maintained alongside HomotopyContinuation.jl itself.

HomotopyGetsReal.jl fills this gap, and does more than reimplement Bertini_real's published algorithms as originally described. It is a native Julia library, built directly on HomotopyContinuation.jl and GLMakie, that reimplements the core of the Bertini_real pipeline — critical point computation, singularity classification, midslice-based curve tracking, and patch-based surface sweeping — as a single-process, dependency-light package installable through Julia's package manager, and extends it in two concrete respects. First, decomposition is no longer restricted to a fixed ambient coordinate: a generic (randomized) projection mode resolves configurations where axis-aligned slicing fails outright, such as a torus sliced along its hole axis, where the critical locus is genuinely one-dimensional rather than a finite set of points (§4, §6). Second, the library implements a diagnostic isosingular deflation subsystem — corank estimation, one-round deflation, and a stabilization test that verifies isosingular local dimension — validated exactly against Hauenstein and Wampler's own published examples, including the Whitney umbrella (§3.5). This is, to our knowledge, the most rigorously literature-validated piece of new technical content in this rebuild, and it brings HomotopyGetsReal.jl's theoretical treatment of singular points into closer alignment with Bertini_real's, though not yet its full practical consequence: the diagnosis is not currently used to decompose a singular curve as a first-class geometric object or to improve mesh geometry near one, a real limitation we return to in §6. Beyond these two extensions, the contribution is an architecture built around type-stable parametric data structures and a single centralized configuration object holding every numerical tolerance used by the pipeline (§3), which makes the precision–performance trade-offs of the implementation auditable from one location rather than scattered across the codebase — real engineering quality, but secondary to the two extensions above rather than the paper's headline claim.

This is obtained at a real cost in scope, which we state plainly rather than leave implicit. HomotopyGetsReal.jl does not implement numerical irreducible decomposition or witness-set generation (§2): it takes a single defining system as input and assumes it cuts out an equidimensional real curve or surface, rather than handling reducible or mixed-dimensional solution sets as Bertini_real does. And while its isosingular deflation subsystem correctly diagnoses isosingular local dimension, that diagnosis is not yet used to decompose a singular curve as a first-class geometric object or to alter mesh geometry near a singular point, the way Bertini_real's deflation-based decomposition does (§6). These are genuine scope limitations of the current implementation, not stylistic differences, and they bound the class of problems for which HomotopyGetsReal.jl is currently a substitute for Bertini_real rather than a complement to it.

The remainder of the paper is organized as follows. §2 reviews the homotopy continuation and critical-point theory the pipeline relies on. §3 describes the library's architecture: its type system, its centralized configuration, each of the six pipeline stages from critical-point solving through mesh welding and visualization, and its diagnostic isosingular deflation subsystem (§3.5). §4 demonstrates the library on the unit sphere, the Taubin heart — a surface with a topologically nontrivial singular configuration — and a genus-one torus that requires generic-projection slicing. §5 reports pointwise residual statistics across these validation surfaces and the results of the library's test suite. §6 compares the design of HomotopyGetsReal.jl to that of Bertini_real in detail, including a coordinate-verified case study on the astroid curve. §7 describes the library's license, repository, and installation.

---

## 2. Mathematical background

### 2.1 Homotopy continuation

Numerical algebraic geometry solves a target polynomial system F : ℂⁿ → ℂᵐ not by symbolic elimination but by numerical path tracking: a *start system* G with known solutions is connected to F by a homotopy H(x,t) = (1−t)γG(x) + tF(x), t ∈ [0,1], for a generic complex constant γ (the gamma trick), and each solution of G is tracked as t moves from 0 to 1 by a predictor-corrector scheme. Under genericity of γ, Bertini's theorem guarantees that the tracked paths converge to the isolated solutions of F without path crossings or divergence to infinity along the way, for t < 1. This framework, together with its extension to positive-dimensional solution sets via witness sets and numerical irreducible decomposition, is developed systematically in Sommese–Wampler [2005]; Bates et al. [2013] gives a software-oriented treatment centered on the Bertini system, and HomotopyContinuation.jl [Breiding–Timme 2018], the tracking engine underlying HomotopyGetsReal.jl, provides the same functionality natively in Julia. Throughout the pipeline described in §3, path tracking is used in two distinct roles: to solve the zero-dimensional systems that produce critical, boundary, and singular vertices (§3.4), and, via *parameter homotopies* that vary a single coordinate or slab parameter continuously rather than an artificial homotopy parameter t, to track the positive-dimensional edges and faces of the decomposition itself (§3.6, §3.8).

### 2.2 Real solutions of complex systems

A polynomial system with real coefficients defines a variety V(F) ⊂ ℂⁿ whose real points V(F) ∩ ℝⁿ are, in general, a proper subset of the full complex solution set and are not distinguished by the tracking process itself: continuation is carried out over ℂ, and a solution is accepted as real only after the fact, by testing whether its imaginary part falls below a tolerance (`critical_point_tol` in Table 2). This is unavoidable in floating-point arithmetic, where an exactly real solution and one perturbed by rounding error are numerically indistinguishable without such a threshold; the same tolerance-based acceptance criterion is used, with the same justification, by Bertini_real.

### 2.3 Critical points and the topology of a real fiber

The decomposition strategy common to Bertini_real (§6) and to HomotopyGetsReal.jl rests on a classical fact about projections of real varieties: if π is a coordinate projection (or, more generally, a generic linear functional) and c₁ < c₂ are consecutive values among the critical values of π restricted to V(F) together with any values forced by a bounding region, then the fiber π⁻¹(c) ∩ V(F) ∩ ℝⁿ varies continuously — in particular, its cardinality is constant — for c ranging over the open interval (c₁, c₂). Consequently, a single witness fiber at the midpoint cₘ = (c₁+c₂)/2 determines the topology of the entire interval, and each point of that fiber can be continued (tracked) outward to c₁ and to c₂ without encountering a topological change along the way. This is the mathematical content underlying both the curve algorithm of Lu et al. [2007] and the surface algorithm of Besana et al. [2013], and it is exactly the justification for the MidSlice-First strategy of §3.6: tracking from a midpoint outward toward, rather than away from, the critical values is what guarantees that no topological change is missed between witness and target.

The critical values of π|_V(F) are computed as the solutions of an augmented system {F(x), det J(x)} = 0, where J stacks the Jacobian of F with the differential of π; this is the same augmentation used in §3.4 (e.g. {f, ∂ᵧf} for a plane curve projected onto x). A solution of this augmented system is a critical point of π|_V(F) for one of two distinct reasons: either x is a genuinely singular point of V(F) itself, where the Jacobian of F already drops rank, or x is a smooth point of V(F) at which the tangent space happens to be tangent to the fibers of π (a turning point of the projection, not of the variety). HomotopyGetsReal.jl distinguishes these two cases from a single solve by applying the rank/singular-value test of §3.4 to the Jacobian of F alone, independently of π: a solution is classified `Singular` only in the first case, while turning points that are smooth on V(F) are retained as ordinary `Critical` vertices, since both are needed to correctly bound the slicing intervals but only the former requires the special handling of §3.6–§3.8.

### 2.4 Scope relative to the general theory

The full numerical-algebraic-geometry framework treats a possibly reducible, possibly non-equidimensional solution set by first computing a *numerical irreducible decomposition* [Sommese–Verschelde–Wampler 2001]: a witness set for each irreducible component, obtained via generic slicing and monodromy, from which components can be decomposed independently and non-real conjugate pairs of components can be identified and set aside. Bertini_real is built directly on top of this machinery: its input is a witness set for a single irreducible component, already verified to be self-conjugate, produced by a separate run of Bertini (§6). HomotopyGetsReal.jl does not implement numerical irreducible decomposition or witness-set generation at all: it takes a single defining system F as input and applies the pipeline of §3 directly, implicitly assuming that V(F) is a single equidimensional real curve or surface rather than a union of several components of mixed dimension or conjugate pairs. This is a genuine restriction of scope rather than an implementation detail, and it is the reason the validation examples of §4–§5 are all irreducible hypersurfaces cut out by a single polynomial.

---

## 3. Architecture

HomotopyGetsReal.jl is organized as a single module of approximately 5800+ lines of Julia, built on top of HomotopyContinuation.jl for numerical path tracking and GLMakie for visualization. The library implements the same conceptual pipeline as Bertini_real — decomposition of a real algebraic curve or surface into a cell complex via homotopy continuation — but is organized around three design principles that depart from that implementation: type-stable parametric data structures, centralized numerical configuration, and a single reusable path-tracking engine shared by the one- and two-dimensional stages of the pipeline.

### 3.1 Design principles

**Type stability.** All geometric objects are parametrized by a floating-point type `T <: AbstractFloat`, allowing the same code path to run in `Float64` for production use or `BigFloat` for precision-sensitive validation (e.g. singularity classification near ill-conditioned Jacobians), without code duplication.

**Centralized configuration.** Every numerical tolerance used across the pipeline is a field of a single struct, `HomotopyConfig{T}`, passed explicitly to every function that requires it. No tolerance is hard-coded at a call site. This is a deliberate departure from the more common pattern of scattering magic constants throughout a numerical codebase, and it is what makes the precision–performance trade-offs of the library auditable from a single location (Table 2).

**Engine reuse.** The one-dimensional curve tracer and the patch-based surface tracer are both built on top of a single adaptive path-tracking primitive (§3.7), configured differently but sharing the same bisection-on-failure and singularity-detection logic. This avoids the duplication of tracking logic that would otherwise be required to support both the 1D and 3D stages of the pipeline.

### 3.2 Core data types

Table 1 summarizes the parametric structs used to represent the decomposition. `NativeVertex{T}` carries a classification (`Critical`, `Boundary`, `Singular`, or `Artificial`) together with an untyped metadata dictionary used to record diagnostic information — Jacobian rank, singular values, the tolerance under which the vertex was accepted, and provenance flags such as `:endpoint_fallback` for vertices synthesized when a tracked path fails to land cleanly. `Edge{T}` and `Face{T}` store, respectively, arc-length-resampled real point sequences and welded triangular mesh data.

**Table 1. Core parametric types.**

| Type | Role |
|---|---|
| `VertexType` | Enum: `Critical`, `Boundary`, `Singular`, `Artificial` |
| `NativeVertex{T}` | Complex coordinates, classification, diagnostic metadata |
| `Edge{T}` | Endpoint vertex ids, resampled real points, singularity flag |
| `Face{T}` | Seed z-slice, boundary edges, welded mesh vertices/topology |
| `HomotopyConfig{T}` | All numerical tolerances (Table 2) |

**Table 2. Fields of `HomotopyConfig{T}` and their role in the pipeline.**

| Field | Default | Role |
|---|---|---|
| `critical_point_tol` | 10⁻⁶ | Imaginary-part cutoff for solutions; Newton polish target |
| `boundary_point_tol` | 10⁻⁵ | Bounding-box containment slack |
| `vertex_match_tol` | 10⁻⁴ | Vertex identity / clustering; mesh weld threshold |
| `jacobian_rank_tol` | 10⁻⁸ | Singular-value cutoff for numerical rank |
| `singular_value_threshold` | 10⁻⁶ | Minimum-singular-value gate for `Singular` classification |
| `path_tracker_precision` | 10⁻¹⁰ | Minimum step size passed to the underlying tracker |
| `patch_transversality_cos_tol` | 0.9 | Re-anchoring gate for surface patches (≈26°) |
| `max_path_steps` | 1000 | Step budget per tracked path |
| `edge_sample_density` | 50 | Arc-length resampling count per edge |
| `midslice_sample_density` | 100 | z-samples per direction in the surface sweep |

### 3.3 Pipeline overview

The decomposition pipeline consists of six stages, summarized in Figure 1: (1) a solver core computing critical points, bounding-box intersections, and singularity classification; (2) a one-dimensional curve decomposer following a "MidSlice-First" tracking strategy; (3) a generalized path-tracking engine shared by stages 2 and 4; (4) a patch-based surface sweep operating on z-slices; (5) mesh welding with gradient-consistent triangle winding; and (6) GLMakie-based visualization.

> **[Figure 1 — pipeline diagram, TikZ, not reproduced here]** Two parallel lanes (curve on top, surface on bottom) share a single adaptive path-tracking engine at Step 4 and draw every numerical tolerance from one `HomotopyConfig{T}` object (outer frame). Dashed boxes are opt-in: the generic-projection preprocessing step before Step 1, and the diagnostic isosingular deflation branch off Step 5, which stamps metadata but does not feed back into the main flow.

We describe this as "the same six steps" as Bertini_real (§6) because the correspondence is genuinely useful for orientation, not because it is exact: the two pipelines diverge in specific, identifiable ways, most visibly in how the bounding region is constructed (an axis-aligned box here versus Bertini_real's bounding sphere), in how "merge" and "sample" are packaged relative to one another, and in how surface faces are constructed (a z-sweep engine here versus Bertini_real's midtrack approach). A reader comparing the two codebases step by step will find these divergences; we name them here rather than let a careful referee find them first.

### 3.4 Solver core: critical points, boundary intersection, and singularity classification

Given a polynomial system F, critical points are computed by augmenting F with its partial derivatives when necessary (a single equation in three variables is augmented to {f, ∂ₓf, ∂ᵧf}), solving in `Float64` via HomotopyContinuation.jl, discarding solutions whose imaginary part exceeds `critical_point_tol`, and, when T ≠ `Float64`, refining the surviving real solutions with a fixed-iteration Newton polish. Boundary vertices are obtained analogously by fixing each variable at each face of the bounding box and solving the resulting reduced system; this step is restricted to plane and space curves (nₑ = nᵥ − 1, nᵥ ∈ {2,3}), since intersecting a curve with a face of the bounding box is well posed while intersecting a surface with a one-dimensional edge of the box is not.

Every vertex is classified by evaluating the Jacobian of F at the candidate point and computing its singular values σ₁ ≥ ... ≥ σₖ via singular value decomposition (using GenericLinearAlgebra.jl when T is not a hardware float, e.g. `BigFloat`). Writing r = |{i : σᵢ > `jacobian_rank_tol`}| for the numerical rank and r_exp for the number of defining equations, a vertex is classified `Singular` if r < r_exp or σₖ < `singular_value_threshold`, and `Critical` or `Boundary` otherwise, according to its origin. The two gates are intentionally independent: a full-rank Jacobian with a small trailing singular value indicates near-singular but numerically distinguishable behavior, which the rank gate alone would miss. Duplicate vertices arising from redundant solves are merged by `cluster_vertices` under `vertex_match_tol`, with `Singular` classification taking precedence in ties.

### 3.5 Isosingular deflation

The rank/singular-value test of §3.4 tells us *that* a vertex is singular but not, on its own, anything about the local structure of V(F) there. HomotopyGetsReal.jl includes a second, diagnostic-only subsystem that answers a sharper question at a `Singular` vertex: what is its isosingular local dimension, in the sense of Hauenstein and Wampler?

We state the rule plainly here, rather than leave it to be inferred from the worked examples below. Isosingular local dimension is defined operationally through the *corank sequence* produced by repeated deflation: each round of deflation augments the current system with new equations capturing the vanishing of certain minors of its Jacobian at the point in question, and the corank — the rank deficiency of that Jacobian — is recorded after each round. If this sequence reaches corank 0 after finitely many rounds, the point is isolated under repeated deflation and its isosingular dimension is 0. If the sequence instead plateaus at a nonzero value k rather than continuing to decrease, its isosingular dimension is k, and the point lies on a genuinely k-dimensional piece of the singular locus that deflation does not reduce to an isolated point. This is Hauenstein and Wampler's own definition; the stabilization test described below is what turns "plateaus" from an informal description into a precise, checkable criterion.

Four primitives implement this. `estimate_corank` measures the rank deficiency of the Jacobian at the point, keeping the row-rank and column-rank conventions deliberately distinct rather than conflating them, a distinction found necessary in practice. `deflate_once` performs a single round of isosingular deflation directly on the bare defining system, using a minor-based construction that requires no witness-set or witness-slice machinery to set up. `verify_isosingular_dimension` is the stabilization test that decides whether a candidate isosingular dimension is authoritative rather than merely plausible. `resolve_isosingular_dimension` orchestrates repeated rounds of deflation until `verify_isosingular_dimension` accepts a dimension or a round budget is exhausted.

**Ground-truth validation.** We validated this subsystem exactly against Hauenstein and Wampler's own worked examples. On the node curve (y² − x² = 0) the corank sequence is [2,0] in one round; on the cusp curve (y² − x³ = 0) it is [2,1,0] in two rounds. On the Whitney umbrella (x² − y²z = 0), the tip (0,0,0) yields corank sequence [3,2,0] and isosingular dimension 0 (an isolated point, resolved in two rounds), while a point (0,0,1) on the umbrella's one-dimensional singular z-axis yields corank sequence [3,1,1,...] — a genuine plateau that never reaches 0 — and isosingular dimension 1. This is the correct behavior for a point on a positive-dimensional singular locus, not a failure to converge, and it matches Hauenstein and Wampler's own published examples exactly.

**Why the stabilization test is nontrivial.** A naive criterion — accept a dimension once two consecutive corank estimates agree — is shown in Hauenstein and Wampler's own analysis to be necessary but not sufficient for correct stabilization; their f_{k,l} counterexample family exhibits a plateau that is not yet stable. Two alternative, cheaper verification strategies were tried and rejected on evidence before settling on the current design: discarding equations from the accumulated deflated system produces a false positive on the Whitney umbrella's singular axis, and replacing the full system with a randomized reduction R(f) = A·f was measured to produce a false-positive rate of four cases out of six on inputs known not to have stabilized. `verify_isosingular_dimension` instead draws a fresh set of d generic hyperplanes through the point on every retry — rather than reusing the same draw across attempts — and checks the result against the full, undiscarded, original system rather than a reduced or randomized one. Acceptance requires both a successful tracking return code *and* a residual check against the true undeflated system; the residual check alone was shown necessary, since trusting the tracking return code alone would have silently accepted the false positives above.

**Surface-level validation.** On the Taubin heart (§4.2), every one of 22 real deflation firings encountered during a full surface decomposition with deflation enabled was inspected individually, including firings that resolve and are absorbed during mesh clustering rather than surviving to the final output. Seventeen resolved in a single round; the remaining five, all at the slice-level silhouette extremes of the heart's cusp, required four rounds (corank sequence [2,1,1,1,0]) — a genuine local-degeneracy signal at a geometrically distinguished point, not a red flag. Across all 22 firings, zero were classified ambiguous and zero exhausted their round budget without resolving.

**Scope.** Deflation in HomotopyGetsReal.jl is diagnostic only. A `Singular` vertex on which `deflate_once`/`resolve_isosingular_dimension` is run has its isosingular dimension, verdict, and corank sequence stamped into its metadata dictionary (§3.4), but its coordinates and the surrounding mesh are not altered as a consequence. We state this as a limitation rather than a design choice to be proud of: Bertini_real uses the analogous diagnosis to decompose a singular curve as its own geometric object (§6); this library does not yet do so.

### 3.6 One-dimensional decomposition: MidSlice-First

Rather than tracking outward from the singular and critical vertices themselves, the curve decomposer tracks outward from smooth witness points located at the midpoint between adjacent distinct x-coordinates of the vertex set (a "midslice"). This ordering is deliberate: a midslice positioned strictly between two consecutive critical or boundary x-values is smooth by construction, whereas step-size control in a predictor-corrector path tracker degrades sharply in the neighborhood of a singular fiber. Singular and critical vertices are therefore only ever approached as *targets* of a track, never used as its starting point. For each midslice, every real root in y seeds a bidirectional track toward the adjacent x-values on both sides (`connect_the_dots!`); each resulting edge is then resampled to a fixed arc-length density (`sample_edge`). Tracks that fail to land within tolerance of a known vertex synthesize an `Artificial` vertex, tagged with provenance metadata, rather than silently discarding the path segment.

### 3.7 Path-tracking engine

Both the curve decomposer and the surface sweep are built on a single adaptive tracking primitive operating on complex state vectors (of length 1 for plane curves, length 2 for surface patches). Given a parameter homotopy and a start/target parameter pair, the tracker advances the homotopy parameter and bisects the step whenever the result is non-finite, exceeds a poor-accuracy threshold, or is flagged as near a singular fiber by an auxiliary rank/singular-value test on the Jacobian — subject to a step budget and a minimum step width. Termination is not gated on HomotopyContinuation.jl's own success flag, since paths legitimately terminate near branch points without satisfying it. The two call sites differ chiefly in compilation mode: the curve tracer compiles the full homotopy system once per edge (`compile=:all`), while the surface tracer, which builds a fresh linear patch system at every anchor point, uses interpreted evaluation (`compile=:none`) to avoid recompilation overhead that would otherwise dominate runtime.

### 3.8 Surface sweeping: z-slicing and patch-based face tracking

The three-dimensional stage decomposes a surface F(x,y,z) = 0 by slicing along z. Critical z-values are computed by treating z as a third free variable in the critical-point solve of §3.4 and clustering the resulting values; together with the bounding-box limits, these partition the z-range into slabs. Within each slab, a "robust" midslice routine locates a curve at a representative z-value, retrying with a small perturbation (bounded by `max_z_mid_retries`) if the initial slice lands on a degenerate or singular configuration — a situation that arises, for instance, at z=0 for the Taubin heart surface (§4.2), where the defining polynomial factors non-trivially in that plane.

Each edge of the midslice curve seeds a surface patch: given the gradient ∇F at an anchor point, a local linear patch a(x − x₀) + b(y − y₀) = 0 with (a,b) = (F_y, −F_x) is constructed so that tracking along z follows the surface transversally to the gradient's projection onto the xy-plane. The patch is re-anchored whenever the cosine between the gradient direction at the current anchor and at the original anchor drops below `patch_transversality_cos_tol`, preventing the linear approximation from drifting too far from the true tangent direction as the track progresses in z.

### 3.9 Mesh welding and gradient-consistent winding

Face meshes produced independently by each patch track are merged into a single mesh by `weld_mesh`: vertex coordinates from all faces are pooled and clustered under `vertex_match_tol`, face triangles are remapped onto the resulting global vertex indices, and degenerate triangles are discarded. Triangle winding is then normalized so that the geometric normal (p₂ − p₁) × (p₃ − p₁) has non-negative inner product with ∇F evaluated at p₁, reversing the vertex order when it does not. This "gradient-consistent" convention guarantees a globally coherent orientation for the welded mesh even though individual face tracks are produced independently and without a shared winding convention.

**Watertightness.** We state this precisely rather than describe the welded mesh as unconditionally watertight, which it is not. On the Taubin heart, an instrumented count of naked (unmatched) mesh edges finds 188 before incidence-based stitching is enabled. Enabling incidence-based stitching closes the naked edges arising at the surface's fold-type or point-type singular features (its two tips) completely, while naked edges at the two boundaries where more than two faces meet along a shared, multi-face edge (the singular notch and a saddle pair) are reduced substantially but not to zero: the residual count there ranges from 31 to 35 across repeated decomposes, itself not exactly reproducible run to run, exhibiting the same kind of cross-run solver jitter documented elsewhere in this paper for the torus (§4.3). We record this as a genuine, current scope limitation of the welding strategy at multi-face boundaries, not as a solved problem understated by the word "watertight."

### 3.10 Visualization

GLMakie-based plotting is provided for both stages of the pipeline: `plot_curve_decomposition` renders the planar edge set with vertices colored and marker-coded by classification (critical, boundary, singular, artificial), and `plot_surface_decomposition` renders either a welded mesh or a raw per-face collection, with coloring by coordinate or by an arbitrary scalar function. Visualization of the raw per-face collection does not apply the winding correction of the previous section and is intended for diagnostic inspection rather than for presentation-quality output.

---

## 4. Usage examples

We illustrate the library on three surfaces chosen to exercise distinct aspects of the pipeline: the unit sphere, a smooth surface with no critical fibers beyond its two poles; the Taubin heart, a surface with a cusp-like indentation whose defining polynomial degenerates along a coordinate plane; and a genus-one torus, which requires the generic projection mode of §6 rather than the library's default fixed-coordinate slicing. All timings below were measured on a single machine in a persistent Julia session; first-call timings include just-in-time compilation of HomotopyContinuation.jl's tracker and are reported separately from steady-state ("warm") timings, consistent with common practice for Julia numerical software.

### 4.1 Unit sphere

```julia
using HomotopyContinuation, HomotopyGetsReal, CairoMakie

@var x y z
F = System([x^2 + y^2 + z^2 - 1], variables = [x, y, z])
cfg = HomotopyConfig{Float64}()
vertices, edges, faces, mesh = decompose_3d_surface(F, cfg)
fig = plot_surface_decomposition(mesh; color_by = :z, cfg = cfg)
CairoMakie.save("sphere.pdf", fig)
```

At the library's default sampling density (`edge_sample_density` = 50, `midslice_sample_density` = 100), this call identifies 2 critical vertices (the poles at z = ±1), 2 edges, and 2 faces, producing a welded mesh of 19,504 vertices and 39,004 triangles in 3.94 s (warm) after a first-call cost of 23.93 s. At the coarser sampling density used throughout the test suite (`edge_sample_density` = 6, `midslice_sample_density` = 8), the same topological output (2 vertices, 2 edges, 2 faces) is produced with a mesh of 152 vertices and 300 triangles, illustrating that sampling density controls mesh resolution without affecting the combinatorial structure of the decomposition.

*(Figure 2: decomposition of the unit sphere, colored by z-coordinate.)*

### 4.2 Taubin heart

```julia
using HomotopyContinuation, HomotopyGetsReal, CairoMakie

@var x y z
f = (x^2 + (1.2y)^2 + z^2 - 1)^3 - x^2*z^3 - 0.1*(1.2y)^2*z^3
F = System([f], variables = [x, y, z])
cfg = HomotopyConfig{Float64}(
    bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3),
    edge_sample_density = 8, midslice_sample_density = 8)
vertices, edges, faces, mesh = decompose_3d_surface(F, cfg)
fig = plot_surface_decomposition(mesh; color_by = :z, cfg = cfg)
CairoMakie.save("taubin.pdf", fig)
```

The Taubin heart is a substantially richer test case than the sphere: `compute_critical_z_slices` identifies four critical z-values, z_crit ≈ {−1.0, 1.0, 1.0648, 1.2367}, matching the values used as hard assertions in the test suite to within 10⁻⁷. The plain return value of `decompose_3d_surface` produces 14 vertices — 10 classified `Critical` and 4 synthesized as `Artificial`, none classified `Singular` — together with 14 edges and 14 faces, welding to a mesh of 1,638 vertices and 3,124 triangles. The four artificial vertices, absent in both the sphere and ellipsoid examples, are produced by the endpoint-fallback mechanism of §3.6 when a tracked path does not land within tolerance of a known vertex; their presence here and their absence on the two quadric examples is consistent with the greater topological complexity of the heart surface near its cusp.

This plain result, however, understates why the Taubin heart was chosen as a validation fixture in the first place: the heart's two real cusps (its top and bottom tips) are genuine singular points of the surface, and the plain decomposition tuple does not report them as such. Enabling the surface-level incidence diagnostic and taking the union of the plain vertex set with its reported critical vertices yields the paper's headline result for this example: 20 vertices in total — 14 `Critical`, 4 `Artificial`, and **2 genuinely `Singular`** vertices at the cusps, at identical mesh statistics (1,638 vertices, 3,124 triangles; residual statistics below are unaffected, since the overlay adds classified vertices without altering the mesh). We report this combined view, rather than the plain tuple, as the representative outcome for this fixture, since it is what actually substantiates a topologically nontrivial singular configuration; we also record the plain-tuple behavior explicitly as a usability note — a caller who inspects only the default return value will not see the real cusps as `Singular` unless the separate diagnostic is also invoked (§6 discusses this alongside Bertini_real's own, more integrated treatment of singular structure). Wall-clock time was 18.52 s for `compute_critical_z_slices` and 24.86 s for the full `decompose_3d_surface` call (first-call timings).

*(Figure 3: decomposition of the Taubin heart surface, colored by z-coordinate.)*

### 4.3 Torus

```julia
using HomotopyContinuation, HomotopyGetsReal, CairoMakie, Random

@var x y z
f = (x^2 + y^2 + z^2 + 3)^2 - 16*(x^2 + y^2)
F = System([f], variables = [x, y, z])
cfg = HomotopyConfig{Float64}(projection = :random, rng = Xoshiro(42))
vertices, edges, faces, mesh = decompose_3d_surface(F, cfg)
fig = plot_surface_decomposition(mesh; color_by = :z, cfg = cfg)
CairoMakie.save("torus.pdf", fig)
```

The torus (x²+y²+z²+3)² − 16(x²+y²) = 0 is a genus-one surface, included specifically to validate a different axis from the other three examples: the sphere and ellipsoid are genus-zero correctness controls and the Taubin heart exercises singular structure, but none of the three has a topologically nontrivial critical locus under z-slicing. The torus does: its fold circles at z = ±1 form a genuinely one-dimensional critical locus under the default fixed-coordinate slicing of §3.8, which the isolated-point critical-value solver cannot represent, producing unusable output (residuals as large as 10⁸) rather than a graceful failure. Switching to the generic (randomized) projection mode described in §6 resolves this outright: the decomposition produces 8 vertices, all correctly classified `Critical` (the torus has no real singularities), and 8 edges and faces, welding to a mesh of 78,340 vertices and 160,381 triangles.

Unlike the sphere, ellipsoid, and Taubin heart, which reproduce identical mesh statistics across repeated runs, the torus result is sensitive to run-to-run variability with a precisely identified cause, not an open question. HomotopyContinuation.jl's `solve` exposes no seeding path into the uniform lifting sampler that MixedSubdivisions.jl uses internally to construct its polyhedral start system; that sampler instead draws from Julia's global `default_rng()`, which is itself seeded from OS entropy independently at every process launch. `decompose_3d_surface`'s own `rng` keyword reaches only the projection-matrix draw, not this internal solve, so fixing that keyword has no effect on this source of variability. We confirmed the mechanism directly: two otherwise-independent runs produce bit-identical output once Julia's global RNG is seeded explicitly (`Random.seed!`) before the call, rather than relying on the pipeline's own `rng` keyword alone. Under the pipeline's current, unseeded default behavior, repeated verification runs at the same configuration and the same fixed projection seed produced mesh counts varying by a few tenths of a percent (e.g. 78,360/168,962 and 78,339/160,346 triangle/vertex pairs were also observed). We report 78,340 vertices and 160,086 triangles as the number tied to the specific run underlying Figure 4, not as an exactly reproducible invariant of the configuration, and we likewise report its residual statistics as characteristic rather than exact: median pointwise residual on the order of 1.4×10⁻⁶ and maximum on the order of 2.0×10⁻⁵, one to two orders of magnitude looser than the sphere and ellipsoid, consistent across the observed runs. The fix — seeding Julia's global RNG explicitly before calling `decompose_3d_surface` — is straightforward but not applied by default in this library today, so the numbers reported here reflect current default behavior rather than an inherent limit on reproducibility.

*(Figure 4: decomposition of a genus-one torus via generic projection, colored by z-coordinate.)*

---

## 5. Validation

Correctness is assessed in two complementary ways: pointwise residuals of the decomposed mesh against the defining polynomial, evaluated on the unit sphere, an asymmetric ellipsoid, the Taubin heart, and a genus-one torus; and a formal test suite exercising the library's internal invariants.

### 5.1 Residual distributions

For each example we evaluate |F(p)| at every vertex p of the final welded mesh and summarize the distribution in Table 3. Figure 5 shows the corresponding histograms on a logarithmic scale. Sphere, ellipsoid, and Taubin heart figures reproduced exactly across repeated runs; the torus figures did not (solver run-to-run nondeterminism; see §4.3 for the observed range) and are reported as characteristic order-of-magnitude values rather than exact statistics, marked with *.

**Table 3. Pointwise residual statistics of |F(p)| on the welded mesh.**

| Example (density) | n | mean | median | p₉₀ | p₉₉ | max |
|---|---|---|---|---|---|---|
| Sphere (production) | 19,504 | 2.85×10⁻⁸ | 2.62×10⁻⁸ | 5.32×10⁻⁸ | 7.05×10⁻⁸ | 7.58×10⁻⁷ |
| Sphere (coarse) | 152 | 2.31×10⁻⁸ | 2.14×10⁻⁸ | 4.92×10⁻⁸ | 5.61×10⁻⁸ | 5.61×10⁻⁸ |
| Ellipsoid (production) | 19,500 | 2.70×10⁻⁷ | 2.79×10⁻⁸ | 5.85×10⁻⁸ | 8.99×10⁻⁸ | 1.12×10⁻⁴ |
| Ellipsoid (coarse) | 160 | 5.26×10⁻⁷ | 2.57×10⁻⁸ | 5.96×10⁻⁸ | 1.000×10⁻⁵ | 1.000×10⁻⁵ |
| Taubin heart | 1,638 | 1.034×10⁻⁷ | 3.44×10⁻⁸ | 1.847×10⁻⁷ | 2.000×10⁻⁶ | 2.42×10⁻⁶ |
| Torus* | 78,340 | 2.211×10⁻⁶ | ~1.4×10⁻⁶ | 5.491×10⁻⁶ | 1.008×10⁻⁵ | ~2.0×10⁻⁵ |

*(Figure 5: distribution of pointwise residuals |F(p)| (log scale, axis clamped to [10⁻¹⁰, 10⁻³]) for the sphere, ellipsoid, and Taubin heart examples. A handful of mesh points (2 each on the sphere and ellipsoid) converge to a literal residual of 0 — genuine bit-exact landings on true critical points, such as the sphere's own (±1,0,0), not an error — and are binned at the axis's lower edge rather than plotted at their true value, which would otherwise stretch the axis to ~10⁻³⁰⁰.)*

Across the sphere, ellipsoid, and Taubin heart, the median residual clusters near 2×10⁻⁸, roughly 50 times tighter than the `critical_point_tol` setting of 10⁻⁶ used to accept tracked solutions. The maximum residual on the production-density ellipsoid mesh (1.12×10⁻⁴) exceeds, by a small margin, the 10⁻⁴ pointwise gate used as a pass/fail criterion in the test suite; that gate is applied only to the coarser mesh used in automated testing (maximum residual 1.000×10⁻⁵ at that density), and does not constitute a claim about arbitrary sampling densities. We report this explicitly rather than round it away: the effect is consistent with accumulated floating-point error at higher sample density in regions of higher curvature, remains four orders of magnitude below the geometric scale of the bounding box, and does not correspond to a topological error in the decomposition (vertex, edge, and face counts are identical to the coarse run). The torus's residuals are one to two orders of magnitude looser still (~2.0×10⁻⁵ maximum) and are, as noted above, themselves only characteristic rather than exactly reproducible; this is consistent with it being the one example that requires the generic projection mode rather than the library's simpler default path (§6).

**A regression caught by residual testing.** The residual figures above reflect a fix, made during this rebuild, to `sample_edge`, the arc-length resampling routine used by both the curve and surface stages (§3.6, §3.8). The routine previously approximated a tracked edge between consecutive arc-length samples by a straight chord rather than continuing to track the underlying curve, an error that is essentially invisible on sharply curved features such as the neighborhood of a cusp, where chord and arc coincide closely, but degrades pointwise accuracy substantially on smooth, gently curved edges, such as a circle, where consecutive samples can cut visibly inside the true curve. On the worst affected test case, the maximum pointwise residual was 0.4998 before the fix and 7.3×10⁻⁷ after it, a reduction of roughly six orders of magnitude. We record this here as a concrete illustration of what the residual checks in this section are for: the error was caught by exactly this kind of pointwise residual testing, not by visual inspection, and would have been easy to miss on test geometry dominated by sharp features.

*(Figure 6: the `sample_edge` regression on a circle — chord approximation before the fix (left, max residual 0.4998) versus correct arc-length tracking after it (right, max residual 7.3×10⁻⁷).)*

### 5.2 Robustness at a degenerate slice

The Taubin heart exercises the robust slicing mechanism of §3.8 directly: the naive midslice of the slab z ∈ [−1, 1] at z_mid = 0 is rejected, since at z = 0 the defining polynomial factors as (x² + 1.44y² − 1)³ and the resulting curve slice is topologically degenerate (4 vertices fall back to `Artificial`, 2 are classified `Singular`). The retry mechanism perturbs the candidate slice up to `max_z_mid_retries` times and accepts z_mid = 0.06 after 5 retries, after which decomposition proceeds without further incident. This is, to our knowledge, the only slab among the four validation examples where the retry mechanism is triggered, and it occurs precisely at the z-value where the surface's algebraic structure is degenerate — the case the mechanism was designed for.

### 5.3 Test suite

The formal test suite comprises 537 `@test`/`@inferred` assertions across thirteen test files. A live, full run (environment variable `HOMOTOPYGETSREAL_RUN_SLOW_TESTS=1` set, exercising every test file including the Taubin-heart- and deflation-specific suites) reports 537/537 passing, 0 failures. We observed a documented ±1 jitter in the total assertion count across otherwise-identical runs, traced to the astroid fixture's cusp-detection step firing a variable number of `@test` invocations inside a loop rather than to any regression; we cite 537/537 as the representative figure rather than a strictly invariant one.

**Table 4. Test suite results (full configuration, live run).**

| Configuration | Pass / Total | Failures |
|---|---|---|
| Full (`HOMOTOPYGETSREAL_RUN_SLOW_TESTS=1`) | 537 / 537 | 0 |

The environment variable continues to gate the slower, Taubin- and deflation-specific assertions behind an opt-in flag so that a default (unset) run remains fast; we do not cite a fast-only total here, since this paper's figures come from a live full run and we prefer not to restate an unverified fast-only count from an earlier snapshot of the codebase.

---

## 6. Comparison with Bertini_real

HomotopyGetsReal.jl reimplements the conceptual pipeline introduced by Bertini_real, which decomposes the real points of a one- or two-dimensional complex algebraic set into a cell complex, following the curve algorithm of Lu, Bates, Sommese, and Wampler and the surface algorithm of Besana, Di Rocco, Hauenstein, Sommese, and Wampler. Both the curve and surface stages of Bertini_real proceed through six named steps — critical points, sphere intersection, slicing, connect-the-dots, merge/singular handling, and refinement — and HomotopyGetsReal.jl deliberately preserves this terminology where the underlying operation is analogous. The two libraries diverge, however, in several load-bearing design decisions, summarized in Table 5 and discussed below.

**Table 5. Design comparison with Bertini_real.**

| | Bertini_real | HomotopyGetsReal.jl |
|---|---|---|
| Runtime dependencies | Compiled Bertini library, Python (sympy), Boost, GMP/MPFR; MPI-parallel | Julia + HomotopyContinuation.jl + GLMakie only |
| Decomposition coordinate | Random real projection(s) π₁ (curves), π₁,π₂ (surfaces), always | Fixed ambient coordinate by default; optional generic/random projection mode, load-bearing for genus-1 examples |
| Bounding region | Intersection with a sphere of chosen center/radius | Axis-aligned bounding box |
| Singular points | Isosingular deflation used to decompose singular curves as first-class output | Isosingular deflation diagnoses local dimension (validated exactly against Hauenstein–Wampler); diagnosis not yet used to decompose a singular curve or alter mesh geometry |
| Refinement | Separate `sampler` program; fixed-count or adaptive tolerance-based | Static per-run sampling density in `HomotopyConfig` |
| Demonstrated scale | Curves up to 14 variables, degree 630; surfaces up to ~6 variables | Trivariate/genus-1 surfaces of low degree (sphere, ellipsoid, Taubin heart, torus) |

**Dependency footprint.** Bertini_real is a compiled, MPI-parallel C program that requires a library-compiled build of Bertini itself, Python (sympy) for symbolic computation (deflation, symbolic Jacobians and determinants), and the Boost, GMP, and MPFR libraries. HomotopyGetsReal.jl has no external binary dependency: the entire pipeline runs inside a single Julia process on top of HomotopyContinuation.jl, with no separate Python/sympy environment or Bertini compilation step required. This trades Bertini_real's MPI parallelism and mature multiple-precision infrastructure for a substantially simpler installation and dependency surface, consistent with the library's design goal of a native Julia implementation with no Bertini binaries. It is worth being precise about the granularity of that traded-away parallelism rather than leaving it vague: reading Bertini_real's own source directly (the `serial_tracker_loop` routine governing an individual path) confirms that path tracking itself is carried out serially within a given process, with Bertini_real's MPI parallelism obtained by distributing distinct paths across processes rather than by parallelizing the tracking of any single path. HomotopyGetsReal.jl's shared adaptive tracker (§3.7) is, correspondingly, single-process and single-threaded; neither library parallelizes the tracking of an individual path.

**Coordinate choice.** Bertini_real always decomposes with respect to one or two *random* real linear projections π₁,π₂ of the ambient variables, chosen specifically to guarantee genericity (transversality of the projection to the variety's tangent spaces) independent of how the defining system happens to be coordinatized. HomotopyGetsReal.jl defaults to slicing directly along a fixed ambient coordinate — x for curves, z for surfaces — but also supports an optional generic (randomized) projection mode. The fixed-coordinate default is a simplification, not a generalization: correctness of the MidSlice-First strategy (§3.6) implicitly assumes the chosen coordinate is generic enough that critical fibers are isolated and transversality failures are rare, an assumption verified case-by-case for the sphere, ellipsoid, and Taubin heart (§5) but not guaranteed in general. A defining system pathologically aligned with the z-axis — for instance, one invariant under z ↦ −z with a degenerate critical fiber, as already observed at z=0 for the Taubin heart — requires the ad hoc perturbation-and-retry mechanism of `_robust_slice_at_z` rather than being avoided by construction. The generic projection mode exists precisely for the cases this default cannot handle at all rather than merely handles awkwardly: a torus sliced along its hole axis has fold circles at z = ±1 that form a genuinely one-dimensional critical locus, which the isolated-point critical-value solver of §3.4 cannot represent under the fixed coordinate, producing unusable output (residuals as large as 10⁸). Switching to a generic projection resolves this outright and is, to our knowledge, the only one of the four validation surfaces for which the fixed-coordinate default is not merely imprecise but fails completely.

**Singularity handling.** Both libraries now use isosingular deflation to reason about singular points, but the mechanisms differ in ways worth stating precisely rather than glossing as "the same idea implemented twice." Bertini_real's own stabilization test, `isosingularDimTest`, verifies a candidate isosingular dimension via a slice-moving homotopy that reuses the *same* draw of random hyperplanes across all retries (up to a fixed iteration cap) and checks the result against a reduced, randomized internal system. HomotopyGetsReal.jl's own `verify_isosingular_dimension` (§3.5) instead draws a *fresh* set of hyperplanes on every retry and checks against the full, original, undeflated system directly — a deliberate, evidence-driven departure motivated by measured false-positive rates in the reduced-system approach, not an oversight. Bertini_real also requires an externally generated witness set for the whole component, obtained from Bertini's own numerical irreducible decomposition, before any specific singular point is identified; HomotopyGetsReal.jl's `deflate_once` works directly on the bare defining system, with no witness-set stage, which we confirmed by exactly reproducing Hauenstein and Wampler's own published deflation sequences on the Whitney umbrella without constructing a witness set at all.

Despite this parity in diagnostic mechanism, the practical consequence still differs, and this remains the most significant scope gap between the two libraries. Bertini_real uses its deflation diagnosis to decompose the singular curve itself as a first-class geometric object, presented in Brake et al. [2014] as that paper's principal technical advance over the almost-smoothness restriction of the earlier surface algorithm. HomotopyGetsReal.jl's deflation subsystem is diagnostic only (§3.5): a singular vertex is correctly classified and its isosingular dimension correctly computed, but that computation does not yet feed back into decomposing the local branch structure at that point or into the resulting mesh. This is a real scope limitation relative to Bertini_real, not a stylistic difference, and we state it as such rather than implying that matching Bertini_real's diagnostic theory closes the practical gap.

**Incidence and projection-degeneracy diagnostics.** Two smaller findings from a direct comparison against Bertini_real's own source are worth recording. First, HomotopyGetsReal.jl's cross-slab identity mechanism includes a `continuity_ok` incidence-consistency diagnostic with no direct analogue we could find anywhere in Bertini_real's source; however, it is currently write-only, recorded as metadata rather than enforced as a hard gate, so this is a real capability difference but not yet a behavioral one, since neither tool currently rejects output on this basis. Second, both libraries check for projection degeneracy before decomposing, but inconsistently in Bertini_real's case: its own `verify_projection_ok` hard-aborts on curves but only warns on surfaces. HomotopyGetsReal.jl's `_verify_projection_ok` throws in both cases, a uniformly stricter enforcement policy — but it is testing a different condition, not the same one applied more carefully: Bertini_real's check is a determinant condition on the projection, while HomotopyGetsReal.jl's checks directly for vanishing of the augmenting partials. We are careful not to claim these are the same test enforced more consistently; they are related but distinct tests, and only the enforcement policy, not the underlying criterion, is directly comparable.

**A concrete instance: the Taubin heart's cusps.** The Taubin heart (§4) makes the practical consequence of the previous paragraph concrete rather than abstract. Its two real cusps are genuine singular points, but the plain return value of `decompose_3d_surface` does not classify either of them as `Singular`: doing so requires separately enabling the surface-level incidence diagnostic and taking the union of the plain vertex set with the vertices it reports as critical, which then correctly surfaces both cusps. Bertini_real does not have this two-tier structure: because singular-curve decomposition is integrated into its main decomposition pipeline rather than an optional diagnostic invoked separately, a user of Bertini_real does not need to know to ask for singular structure in a second pass to see it in the output. This is a real usability difference in Bertini_real's favor, not merely a cosmetic one, and it is the concrete face of the scope gap described above: HomotopyGetsReal.jl can be made to report the Taubin heart's cusps correctly, but does not do so by default, whereas Bertini_real's integrated treatment does.

**A worked comparison: the astroid.** The astroid curve (x²+y²−1)³ + 27x²y² = 0 gives a fully coordinate-verified, citable comparison. Bertini_real's raw output for this curve reports 32 vertices and 20 edges, which overstates the curve's real topology substantially: inspecting Bertini_real's own internal vertex-type bitflags shows that 18 of the 32 vertices are `Midpoint` witness points (one per edge, an artifact of Bertini_real's own midslice construction rather than a distinct topological feature), 8 more are `New`/`Removed` bookkeeping from sphere compactification, and 14 of the 20 reported edges are degenerate placeholder self-edges. Filtered to genuine topology, Bertini_real finds exactly the same 4 singular cusps as HomotopyGetsReal.jl, at (±1, 0) and (0, ±1). The remaining discrepancy — Bertini_real's 6 non-degenerate edges against HomotopyGetsReal.jl's 4 — is explained, not a sign that either tool is missing an arc: Bertini_real's own random projection happens to place two additional smooth projection-critical points that split two of the four cusp-to-cusp arcs into halves. Direct coordinate matching confirms that HomotopyGetsReal.jl's simpler four-edge decomposition is a complete, correct cycle through the four cusps: the same underlying one-manifold, partitioned differently as a consequence of each tool's own choice of projection, not a completeness gap in either. This exact curve and its four-cusp, four-edge topology is also given as a worked example, independently of both Bertini_real and this library, by Amethyst, Hauenstein, and Wampler, who report "four singular points ... connected by four edges" for the same defining polynomial, providing a third, independent confirmation of the same topological count.

A related citation caution is worth recording here rather than silently avoiding: Bertini_real's own packaged Griffis–Duffy test fixture publishes a reference deflation sequence of 4,1,1 with 588 candidate minors and 8 kept. A live, independent re-run against a real Bertini_real installation reproduced only the leading term of that sequence; the rest could not be reproduced, traced to a type error in Bertini_real's own stock deflation script under the currently installed version of its symbolic-algebra dependency, unrelated to HomotopyGetsReal.jl. We do not cite Bertini_real's published 588/8 figure as an independently verified number in this paper, and note it here only as Bertini_real's own published claim, not independently reproducible with the toolchain available to us.

**Bounding region and refinement.** Bertini_real bounds unbounded components by intersecting with a sphere of user-chosen center and radius, treating the resulting intersection curve as part of the critical curve to be decomposed alongside genuine critical and singular loci. HomotopyGetsReal.jl instead uses an axis-aligned bounding box, with boundary vertices computed by fixing each coordinate at each face of the box. Refinement in Bertini_real is a distinct post-processing stage, performed by a separate `sampler` program that supports both fixed-count and adaptive tolerance-driven resampling of an already-computed coarse decomposition. HomotopyGetsReal.jl instead fixes the sampling density (`edge_sample_density`, `midslice_sample_density`) once, as part of the initial configuration; there is currently no adaptive, error-driven refinement pass analogous to Bertini_real's `sampler`.

**Demonstrated scale.** Bertini_real has been applied to substantially larger and higher-degree problems than those presented in §4: the authors report decomposing a 3-3 Burmester curve in 14 variables of degree 630, and note that their surface decomposer is tractable for systems of roughly six variables without randomization, with an eight-polynomial, ten-variable Burmester surface identified as an open challenge due to the size of the symbolic determinant required by the critical-curve computation. The validation presented in this paper is restricted to trivariate surfaces of low degree. This difference in demonstrated scale reflects the current stage of HomotopyGetsReal.jl's development rather than an architectural ceiling, but it should not be understated: no claim is made here that HomotopyGetsReal.jl matches Bertini_real's demonstrated range of applicability.

---

## 7. Availability

HomotopyGetsReal.jl is released under the MIT license. The source code, test suite, and documentation are hosted at

**https://github.com/JuanCaGC/HomotopyGetsReal.jl**

The package's direct dependencies — HomotopyContinuation.jl, GLMakie.jl (and the underlying Makie.jl), GeometryBasics.jl, GenericLinearAlgebra.jl, Parameters.jl, and Combinatorics.jl (used by the isosingular deflation subsystem, §3.5), together with the Julia standard library (including `Random`, also used by deflation) — are all released under the MIT license (Table 6), and are therefore free software under the Free Software Foundation's definition, as required for JSAG submission. The transitive dependency closure includes one package under the GPLv3 (FLINT, via HomotopyContinuation.jl's arbitrary-precision arithmetic), which is likewise free software and does not conflict with this requirement, and one package, TreeViews.jl, whose repository carries no explicit license file; the latter is a display-only utility with no bearing on the correctness of the numerical pipeline, but its licensing status is technically unresolved and is noted here rather than silently assumed to be unproblematic.

**Table 6. Direct dependencies and their licenses (audited against each package's own `LICENSE` file).**

| Package | License |
|---|---|
| HomotopyContinuation.jl | MIT |
| GLMakie.jl / Makie.jl | MIT |
| GeometryBasics.jl | MIT |
| GenericLinearAlgebra.jl | MIT |
| Parameters.jl | MIT |
| Combinatorics.jl | MIT |
| LinearAlgebra, Statistics, Random (Julia stdlib) | MIT |

### 7.1 Requirements and installation

HomotopyGetsReal.jl requires Julia 1.12 or later; this is enforced by the package's own `[compat]` bounds and by those of its dependencies, and was verified by confirming that dependency resolution fails on Julia 1.10 and 1.11 under the current `Project.toml`.

HomotopyGetsReal.jl is registered in Julia's General registry as `v0.2.1`, confirmed by a live installation in a clean environment:

```julia
using Pkg
Pkg.add("HomotopyGetsReal")
```

This registered release includes every feature described in this paper, including the isosingular deflation subsystem (§3.5), the generic-projection mode (§6), and incidence-based mesh welding, and its `Project.toml` declares the full dependency set required by the current codebase, including Combinatorics.jl and `Random`. A reader who runs `Pkg.add("HomotopyGetsReal")` obtains exactly the feature set this paper describes.

The test suite is run with `Pkg.test()` for the default (fast) configuration and with the environment variable `HOMOTOPYGETSREAL_RUN_SLOW_TESTS=1` set for the full configuration including the Taubin heart examples (§5).

---

## Acknowledgments

Substantial portions of this manuscript's prose, including its exposition and the argumentative structure of its comparison with Bertini_real, were drafted with the assistance of AI tools, primarily Claude (Anthropic) and Cursor. These tools were also used to verify technical claims and numerical results against the author's own codebase and test suite, and to check bibliographic citations against their original sources, including a correction to the title of a cited reference (Amethyst, Hauenstein, and Wampler). The pipeline diagram (Figure 1) was likewise produced through an AI-assisted, iterative drafting process. All technical content, numerical results, and citations in this paper originate from the author's own software and its verified behavior; the author has reviewed, independently verified, and takes full responsibility for every claim, number, and reference appearing here.

---

## References

1. Silviana Amethyst, Jonathan D. Hauenstein, and Charles W. Wampler, *Kinematic synthesis over curves using cellular decompositions and Chebyshev interpolants*, 3rd IMA Conference on Mathematics of Robotics (IMA 2025), Springer, 2026, pp. 27–36.
2. Dennis S. Arnon, George E. Collins, and Scott McCallum, *Cylindrical algebraic decomposition I: The basic algorithm*, SIAM Journal on Computing **13** (1984), no. 4, 865–877.
3. Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler, *Numerically solving polynomial systems with Bertini*, Software, Environments, and Tools, vol. 25, SIAM, 2013.
4. Gian Mario Besana, Sandra Di Rocco, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler, *Cell decomposition of almost smooth real algebraic surfaces*, Numerical Algorithms **63** (2013), no. 4, 645–678.
5. Daniel A. Brake, Daniel J. Bates, Wenrui Hao, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler, *Bertini_real: software for one- and two-dimensional real algebraic sets*, Mathematical Software — ICMS 2014, LNCS vol. 8592, Springer, 2014, pp. 175–182.
6. Daniel A. Brake, Daniel J. Bates, Wenrui Hao, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler, *Algorithm 976: Bertini_real: numerical decomposition of real algebraic curves and surfaces*, ACM Transactions on Mathematical Software **44** (2017), no. 1, 1–30.
7. Paul Breiding and Sascha Timme, *HomotopyContinuation.jl: a package for homotopy continuation in Julia*, Mathematical Software — ICMS 2018, LNCS vol. 10931, Springer, 2018, pp. 458–465.
8. Jonathan D. Hauenstein and Charles W. Wampler, *Isosingular sets and deflation*, Foundations of Computational Mathematics **13** (2013), no. 3, 371–403.
9. Ye Lu, Daniel J. Bates, Andrew J. Sommese, and Charles W. Wampler, *Finding all real points of a complex curve*, Contemporary Mathematics **448** (2007), 183–205.
10. Andrew J. Sommese, Jan Verschelde, and Charles W. Wampler, *Numerical decomposition of the solution sets of polynomial systems into irreducible components*, SIAM Journal on Numerical Analysis **38** (2001), no. 6, 2022–2046.
11. Andrew J. Sommese and Charles W. Wampler, *The numerical solution of systems of polynomials arising in engineering and science*, World Scientific, Hackensack, NJ, 2005.

---

*This markdown is a faithful transcription of the current LaTeX sources (`main.tex` and its `\input` files) as of the last recompile (22 pages, clean). Figures are noted by placeholder since they are binary image files, not reproduced here; the pipeline diagram (Figure 1) is TikZ source, also not reproduced. Section numbers above follow the compiled paper's actual numbering.*
