# HomotopyGetsReal.jl — Technical Repository Summary

**Prepared for:** study-advisor briefing / academic presentation prep
**Package version:** `0.2.1` (`Project.toml:4`), MIT, registered in Julia's General registry
**Julia requirement:** `1.12+` (`Project.toml:27`, `README.md:27`)
**Scope of the software:** real plane curves in ℝ² (1 equation, 2 variables) and real **hypersurfaces** in ℝ³ (1 equation, 3 variables). Nothing else is supported.
**Upstream engine:** HomotopyContinuation.jl `2.17.2` (`Project.toml:22`)

Everything below was read directly out of the working tree. Where a repository document contradicts the live code or another document, the contradiction is flagged in §5 rather than silently resolved.

---

## 0. Orientation: the six-step framework and where it lives

The project is organized around the six algorithmic steps of the original *Homotopy gets real* / Bertini_real pipeline. The mapping from step to code is explicit in the source headers (`src/Topology.jl:3-11`, `src/Config.jl:4-7`):

| # | Step | 1D curve implementation | 3D surface implementation |
|---|------|-------------------------|---------------------------|
| 1 | `compute_critical_points` | `Solver.jl:733` | `SurfaceDecomposition.jl:48` (`compute_critical_z_slices`) |
| 2 | `intersect_bounding_object` | `Solver.jl:825` | (AABB clip via `bbox_z` in `_slab_bounds`, `SurfaceDecomposition.jl:411`) |
| 3 | `interslice` / `MidSlice!` | `Topology.jl:61` (`compute_midslice`) | `SurfaceDecomposition.jl:286` (`_robust_slice_at_z`) |
| 4 | `ConnectTheDots!` | `Topology.jl:177` | `FaceTracking.jl:574` (`track_face`) |
| 5 | `Merge` / `GetMergeCandidates` | `Clustering.jl:40` (`cluster_vertices`) | `SurfaceDecomposition.jl:1581` (`weld_mesh`) |
| 6 | `sample_edge` | `Topology.jl:290` | (per-face mesh grid from `track_face`) |

The package is a **single flat module** — every file is `include`d into one namespace with no per-file `module` blocks, and include order is load-bearing and documented at `src/HomotopyGetsReal.jl:1-38`. Order: `Config → Types → Clustering → Solver → PathTracking → Topology → FaceTracking → Projection → SurfaceDecomposition → Visuals` (`src/HomotopyGetsReal.jl:50-59`).

**Central architectural idea — the "MidSlice-First Strategy"** (`src/Topology.jl:13-19`): path tracking *never* starts from a critical or singular vertex. It always starts from a provably-smooth interior witness point produced by `compute_midslice`, and tracks *outward* toward the vertices. This is why singular points are reached as limits rather than as solve targets.

**Precision architecture.** The package is precision-parametric (`HomotopyConfig{T}`, `NativeVertex{T}`, `Edge{T}`, `Face{T}`) but there is a hard, deliberate precision boundary: *all path tracking happens in Float64 regardless of `T`*, because HomotopyContinuation.jl's `NewtonCache`/`NewtonResult` are hardcoded to `ComplexF64` (`src/Solver.jl:7-12`). BigFloat only enters afterward, through `Solver._newton_polish` (`src/Solver.jl:665`) applied to final landing points. Consequences of this boundary appear throughout §4.

---

## 1. Architecture & Core Modules

### 1.1 File inventory

| File | Lines | Role |
|------|-------|------|
| `src/Config.jl` | 139 | `HomotopyConfig{T}` — all tolerances and knobs |
| `src/Types.jl` | 196 | `VertexType`, `NativeVertex`, `Edge`, `Face` |
| `src/Clustering.jl` | 482 | Vertex/scalar/point clustering, `VertexRegistry` |
| `src/Solver.jl` | 922 | Steps 1–2, Jacobian rank, isosingular deflation |
| `src/PathTracking.jl` | 259 | Adaptive bidirectional tracking engine |
| `src/Topology.jl` | 401 | Steps 3–6, 1D curve pipeline |
| `src/FaceTracking.jl` | 645 | Patch systems, z-sweep face tracking |
| `src/Projection.jl` | 265 | Generic SO(3) projection wrapper |
| `src/SurfaceDecomposition.jl` | 2042 | 3D orchestrator, incidence, mesh welding |
| `src/Visuals.jl` | 416 | Makie plotting |
| **Total src** | **5767** | |

---

### 1.2 `src/Types.jl` — the decomposition data model

```
@enum VertexType  Critical  Boundary  Singular  Artificial     # Types.jl:26-31
```

Semantics (`Types.jl:19-25`):
- `Critical` — critical point of the projection map π
- `Boundary` — intersection with the bounding object (AABB)
- `Singular` — rank-deficient Jacobian at the point
- `Artificial` — synthetic vertex from merging or slicing

**`NativeVertex{T}`** (`Types.jl:48-55`): `id::Int`, `coordinates::Vector{Complex{T}}`, `v_type::VertexType`, `metadata::Dict{Symbol,Any}`. Coordinates are **complex** because homotopy paths are tracked in ℂⁿ; the real decomposition is recovered by filtering on imaginary part.

The config-aware constructor (`Types.jl:84-105`) records *which* tolerance admitted the vertex into `metadata[:tolerance_used]`, dispatching on type: `Critical → critical_point_tol`, `Boundary → boundary_point_tol`, `Singular → singular_value_threshold`, `Artificial → vertex_match_tol` (`Types.jl:93-101`). This is a small but important provenance mechanism — every vertex carries the tolerance that created it.

**`Edge{T}`** (`Types.jl:121-129`): `id`, `left_vertex_id`, `right_vertex_id`, `sampled_points::Vector{Vector{T}}` (real-space), `is_singular::Bool`.

**`Face{T}`** (`Types.jl:167-173`): `id`, `mid_slice_z::T`, `boundary_edges::Vector{Int}`, `mesh_vertices::Matrix{T}`, `mesh_topology::Matrix{Int}`. Note `mid_slice_z` is the *chart* sweep coordinate, not literal world z, when a projection is active (`Types.jl:160-162`).

---

### 1.3 `src/Solver.jl` — critical points, boundary points, isosingular deflation

This is the mathematically densest file. It contains 14 functions, 2 enums, 2 structs.

#### Rank and classification primitives

**`jacobian_rank_info(F::System, point, cfg::HomotopyConfig{T}) → (rank::Int, singular_values::Vector{T})`** — `Solver.jl:47`

Forms the complexified Jacobian J = ∂F/∂x by symbolic differentiation, evaluates at T-precision (`bits = precision(T)`), and computes σ₁ ≥ σ₂ ≥ … by SVD. Numerical rank = `count(>(cfg.jacobian_rank_tol), svals)` (`Solver.jl:55`). This is the standard thresholded-SVD rank estimate. BigFloat SVD requires `GenericLinearAlgebra.jl` because Base has no generic `svdvals` fallback (`Solver.jl:22-27`).

**`_classify_vertex_type(info, cfg, expected_rank::Int, base_type::VertexType) → VertexType`** — `Solver.jl:71`

Two independent, either-sufficient singularity conditions:
1. rank deficiency: `info.rank < expected_rank` (`Solver.jl:77`)
2. near-zero smallest singular value: `minimum(info.singular_values) < cfg.singular_value_threshold` (`Solver.jl:78`)

If either fires → `Singular`; otherwise the point keeps `base_type`.

#### Isosingular deflation (Hauenstein–Wampler)

This is the block with the strongest external validation. It is built directly against Hauenstein–Wampler's `D_det` construction, Definition 5.18, Lemma 3.4, and Algorithm 6.3, with explicit cross-references to Bertini1/BertiniReal source line numbers.

**`estimate_corank(F, point, cfg; expected_rank::Int) → Int`** — `Solver.jl:119`

Returns `expected_rank − rank(J(x))`. The keyword is **mandatory, no default** (`Solver.jl:95-96`) because two incompatible conventions are in use:
- **row-rank convention**, `expected_rank = length(F.expressions)` — "smooth means full row rank of the defining equations"; used by `_classify_vertex_type` and `intersect_bounding_object`
- **column-rank convention**, `expected_rank = length(F.variables)` — "isolated point in ambient space"; used by `deflate_once` and the H–W construction

They coincide only for square systems. Documented limitation (`Solver.jl:112-117`): this is purely first-order and *cannot distinguish singularities of equal first-order corank* — a node and a cusp both report corank 1 under the row convention. Explicitly labelled "expected behavior, not a bug."

**`_corank_plateau_hint(corank_sequence) → Bool`** — `Solver.jl:162`

Enforces H–W Lemma 3.4 monotonicity `N ≥ d_k ≥ d_{k+1} ≥ 0` (throws `ArgumentError` on violation) and returns `true` if the sequence ends in 0 **or** the last two entries are equal. Explicitly a **necessary-but-not-sufficient** pre-filter (`Solver.jl:133-134`), quoting H–W §6: *"a necessary condition for stabilization is that two consecutive terms in the deflation sequence must be equal, but this is not sufficient"* (`Solver.jl:138-141`).

**`deflate_once(F, x0, cfg; expected_rank = length(F.variables)) → (F_new::System, corank_new::Int)`** — `Solver.jl:209`

One H–W isosingular deflation step. Computes corank `c` under the column convention, sets `minorSize = expected_rank − c + 1 = rank(J) + 1` (formula cross-referenced to `isosingular.cpp:195` at `Solver.jl:193`), and appends every `minorSize × minorSize` minor determinant of `jacobian(F)` that is not *structurally* the zero polynomial. The zero test is `iszero(expand(det(...)))` — a CAS structural test, **not** numerical evaluation at `x0`, because all such minors vanish numerically at a corank-`c` point by definition of rank (`Solver.jl:186-190`).

**`_deflation_applicable(F_original, x0, cfg; expected_rank) → Bool`** — `Solver.jl:293`

Stage-4c gate: would `deflate_once` throw here? Necessary because `corank > 0` alone is the *wrong* test — on a bare single-equation curve, even a smooth point reports corank 1 under the ambient convention (`Solver.jl:255-259`). Minimal witness given in the source: `f = x − y³` at the origin (`Solver.jl:260-264`).

**`IsosingularVerdict`** — `Solver.jl:327` — `Verified | NotTerminal | Inconclusive`

**`VerifyResult{T}`** — `Solver.jl:344` — `verdict`, `endpoint`, `attempts`, `inconsistent_count`

**`verify_isosingular_dimension(F_current, F_original, x0, d, cfg; rng) → VerifyResult{T}`** — `Solver.jl:419`

H–W Algorithm 6.3 / Bertini1 `isosingularDimTest` analogue (`isosingular.c:328-450`, cited at `Solver.jl:358-359`). Builds `d` random real linear slices Lᵢ(x) = Σⱼ aᵢⱼ(xⱼ − x0ⱼ) − sᵢ through the real witness `x0`, augments `F_current`, and tracks a parameter homotopy from `s = 0` (hyperplanes through x0) to `s = 1`. Verdicts:

| Verdict | Condition | Line |
|---|---|---|
| `Verified` | `d == 0` short-circuit (no tracking), **or** an attempt with `return_code == :success` **and** endpoint residual against `F_original` ≤ `critical_point_tol` | 427-428, 455-460 |
| `Inconclusive` | **First** attempt that tracks successfully but fails the `F_original` residual check — fires immediately, does not consume remaining budget | 461-463 |
| `NotTerminal` | All `cfg.isosingular_verify_retries` attempts were clean tracking failures with no inconsistent success | 466-469 |

Design decision documented at `Solver.jl:389-391`: checking only tracker success "would produce real false positives" — the residual check against the *original, undeflated* system is what makes the test sound.

**`ResolveVerdict`** — `Solver.jl:500` — `Resolved | Ambiguous | Exhausted`
**`ResolveResult{T}`** — `Solver.jl:518` — `verdict`, `isosingular_dimension`, `endpoint`, `corank_sequence::Vector{Int}`, `rounds::Int`, `F_final::System`

**`resolve_isosingular_dimension(F_original, x0, cfg; rng) → ResolveResult{T}`** — `Solver.jl:598`

Orchestrates H–W Definition 5.18 (isosingular local dimension as the limit of the deflation sequence). Loop: if `d == 0` or the plateau hint fires, verify; `Verified → Resolved`, `Inconclusive → Ambiguous`, `NotTerminal → deflate and continue`. Bounded by `cfg.max_deflations` → `Exhausted` (`Solver.jl:625-628`). A genuine nonzero-limit plateau always costs one extra round beyond the true limit, because the hint needs a *repeated* value (`Solver.jl:563-566`).

**Important scope caveat:** `ResolveResult.F_final` is retained but "Not consumed anywhere yet — Stage 4a is diagnostic-only" (`Solver.jl:513-515`). See §4.6.

#### Newton refinement

**`_newton_polish(F, x0::Vector{Complex{T}}, cfg) → Vector{Complex{T}}`** — `Solver.jl:665`

Hand-rolled `x ← x − J(x)⁻¹F(x)` at T-precision, stopping at `‖F(x)‖ ≤ cfg.critical_point_tol` or 50 iterations. No-op when `T === Float64`. Documented convergence caveat (`Solver.jl:647-662`): near a multiple root — i.e. exactly the points classified `Singular` — the Jacobian is rank-deficient and convergence degrades from quadratic to **linear**, empirically ~0.25× residual shrink per iteration, needing 41 iterations to reach 1e-40 where a simple root needs 4. A `Singular` vertex can therefore sit at or just above `critical_point_tol` after the 50-iteration cap, and this is documented as expected.

#### Steps 1 and 2

**`compute_critical_points(F, cfg; deflate = false, F_original = nothing) → Vector{NativeVertex{T}}`** — `Solver.jl:733`

Two input shapes:
- square (`ne == nv`) — `F` is taken as an already-augmented 0-dimensional system
- raw surface (`ne == 1, nv == 3`) — auto-builds `F_aug = {f, ∂f/∂x, ∂f/∂y}`, i.e. the locus where the gradient is parallel to the z-axis and the surface fails to be a local graph z = g(x,y)

For plane curves the augmentation is done by the *caller*: `decompose_1d_curve` builds `F_aug = System([f, ∂f/∂y], F.variables)` (`Topology.jl:369`) — two equations in two unknowns.

Pipeline: Float64 solve → filter `maximum(abs, imag.(s)) ≤ crit_tol64` (`Solver.jl:772`) → optional Newton polish → classify via Jacobian rank **on `F_aug`** → optional isosingular resolution → `cluster_vertices(candidates, cfg.vertex_match_tol)` (`Solver.jl:790`).

**`intersect_bounding_object(F, cfg; deflate = false, F_original = nothing) → Vector{NativeVertex{T}}`** — `Solver.jl:825`

For each coordinate and each of its two `bbox_*` bounds, fixes that variable, solves the reduced 0-dimensional subsystem, keeps nearly-real solutions, checks the remaining coordinates lie inside the box with `boundary_point_tol` slack (`Solver.jl:884`), embeds back into ambient coordinates, classifies against the **bare** `F` under the row-rank convention, and clusters. Raw surfaces are explicitly out of scope here: "a face intersection is a curve, not isolated points" (`Solver.jl:836-837`).

---

### 1.4 `src/PathTracking.jl` — the adaptive bidirectional tracking engine

Float64/ComplexF64-only by construction; `HomotopyConfig{T}` appears only in `build_tracker`'s signature (`PathTracking.jl:32-38`).

**`is_near_singular(F, point::Vector{ComplexF64}, expected_rank, rank_tol, sv_thresh) → Bool`** — `PathTracking.jl:51`
Mirrors `_classify_vertex_type`'s two-part rule during tracking: `rank_deficient || near_zero_sv` (`PathTracking.jl:68-70`).

**`build_tracker(H_sys, x_start, cfg; compile = :all) → (ph::ParameterHomotopy, tracker::Tracker)`** — `PathTracking.jl:84`
The only place `cfg.path_tracker_precision` reaches HomotopyContinuation, mapped to `TrackerOptions(min_step_size)` (`PathTracking.jl:89`). Explicitly noted as "closest HC lever, not exact accuracy guarantee" (`PathTracking.jl:87`). `compile = :none` is for many short-lived per-anchor systems; `:all` amortizes one reused system.

**`_track_path_segment!(...)`** — `PathTracking.jl:139` — the recursive core.

This function encodes one of the most important design decisions in the package. The landing value is taken from `res.solution` whenever finite and is **deliberately not gated on `HomotopyContinuation.is_success(res)`** (`PathTracking.jl:101-110`): right at a branch point (e.g. approaching a `Critical` vertex where ∂f/∂y = 0) the tracked system has a genuine multiple root, the adaptive internal step size legitimately collapses, and the tracker terminates with `:terminated_step_size_too_small` *after already landing extremely close to the true root*. Treating that as failure would discard a good answer.

Bisection triggers when the landing is non-finite, or `res.accuracy > poor_acc_tol`, or `is_near_singular` fires — each independently sufficient (`PathTracking.jl:120-124`) — provided the step budget is unexhausted and the interval is wider than `min_width` (`PathTracking.jl:166`).

Critical honesty note in the docstring (`PathTracking.jl:126-134`): once the budget or the width floor is reached, "the *last attempted* landing value is accepted as a fallback… there is no 'discard' outcome; every call always produces some point, good or bad." The only downstream defense is the caller's vertex matching against `vertex_match_tol`, which surfaces a bad landing as an unexpected new `Artificial` vertex.

**`track_path(...) → (y_final, path)`** — `PathTracking.jl:193` — sets `expected_rank = length(F.expressions)` (`PathTracking.jl:209`).
**`track_bidirectional(...) → (full_path, y_land_left, y_land_right)`** — `PathTracking.jl:228` — two independent `track_path` calls, assembled as `reverse(path_left) ++ midpoint ++ path_right`.

---

### 1.5 `src/Topology.jl` — the 1D curve pipeline

**`compute_midslice(F, x_left, x_right, cfg) → Vector{Complex{T}}`** — `Topology.jl:61`
Fixes x at the midpoint, solves for y, returns roots with `maximum(abs, imag.(s)) ≤ critical_point_tol` (`Topology.jl:85`). Substitutes `Float64(x_mid)` before solving — required because HC's polyhedral start-system construction has no `Complex{BigFloat}` method (`Topology.jl:71-76`).

**`_resolve_endpoint(F, x_var, y_var, x_val, y_guess, vertices, cfg) → (vertex_id, y_final)`** — `Topology.jl:127`
Newton-polishes the guess against the reduced square system, then searches the **full** vertex list by actual Euclidean distance (never by sorted-list position — a Phase 3 counterexample is cited at `Topology.jl:104-105`), matching within `cfg.vertex_match_tol` (`Topology.jl:151`). On no match it appends a new `Artificial` vertex tagged `metadata[:origin] = :endpoint_fallback` (`Topology.jl:155`). That tag distinguishes "a tracked path failed to close onto any known vertex" (a warning sign) from the benign `Artificial` produced by clustering a coincident `Critical`/`Boundary` pair (`Topology.jl:112-125`). `_robust_slice_at_z` consumes exactly this tag.

**`connect_the_dots!(F, x_left, x_mid, x_right, y_mid, edge_id, vertices, cfg) → Edge{T}`** — `Topology.jl:177`
Builds `H_sys = System([f], variables = [y_var], parameters = [x_var])`, tracks bidirectionally from the midslice witness, resolves both endpoints (mutating `vertices`), and marks `is_singular` if either endpoint is `Singular` (`Topology.jl:223`). Tolerance wiring is explicit at `Topology.jl:200-203`: `rank_tol ← jacobian_rank_tol`, `sv_thresh ← singular_value_threshold`, `poor_acc_tol ← critical_point_tol`, `min_width ← vertex_match_tol`.

**`_project_to_curve(f, fx, fy, x_var, y_var, x0, y0, cfg) → (x, y)`** — `Topology.jl:248`
Single-constraint Gauss–Newton correction: `step = f/(fx² + fy²)`, `(x,y) −= step·(fx,fy)` — the minimal-norm Newton step for one equation in two unknowns, iterating until `|f| ≤ critical_point_tol`, capped at 50 iterations.

**`sample_edge(F, edge, cfg) → Edge{T}`** — `Topology.jl:290`
Resamples to `cfg.edge_sample_density` arc-length-equidistant points. The docstring documents a real 2026-07 bug fix (`Topology.jl:271-289`): this used to be "pure geometry," but a smooth arc that never needed bisection can have a raw path as sparse as **2 points**, and interpolating along the straight chord between them lands every interior point measurably off the curve — **residuals up to 0.28 measured on the Taubin heart fixture**, orders of magnitude above `critical_point_tol`. Every interpolated point is now re-projected via `_project_to_curve`, which is why `F` became a required argument. The paper-artifact measurement of this fix is in `results.json:418-422`: unit circle max residual **0.4998 before → 7.34e-7 after**.

**`decompose_1d_curve(F, cfg; deflate = false) → (vertices, edges)`** — `Topology.jl:359`
Full 1D pipeline: `F_aug = {f, ∂f/∂y}` → `compute_critical_points` → `intersect_bounding_object` → id renumbering with offset → `cluster_vertices(…, vertex_match_tol)` → `cluster_scalars(xs, vertex_match_tol)` to get distinct x-slots → per adjacent interval, `compute_midslice` then one `connect_the_dots!` per y-root → `sample_edge` on every edge.

---

### 1.6 `src/FaceTracking.jl` — sweeping a slice curve into a face

The header (`FaceTracking.jl:1-79`) states the two structural gaps that made plain `PathTracking` insufficient:

1. `f(x,y,z) = 0` at fixed z is **one equation in two unknowns** — underdetermined. Tracking a *specific* point as z sweeps needs a second "patch" equation pinning it, using the local surface gradient `(f_y, −f_x)` as the patch direction. The gradient-based patch was chosen over the old prototype's origin-radial line, which degenerates for curves not centered on the origin (`FaceTracking.jl:9-15`).
2. `track_path`/`track_bidirectional` are **endpoint** trackers — they only record a point when they accept a leaf, so a smooth run may produce 1–3 raw points. Mesh building needs a dense uniform z-grid, hence `track_dense_path` (`FaceTracking.jl:16-25`).

| Function | Line | Role |
|---|---|---|
| `build_patch_system(F)` | 95 | One-time symbolic prep; returns NamedTuple `(f, fx, fy, fz, x_var, y_var, z_var, F_for_tracking)`. `F_for_tracking` reorders variables to `[z, x, y]` so z is the swept parameter |
| `_gradient_at(patch, x0, y0, z0, cfg)` | 158 | T-generic gradient evaluation at `bits = precision(T)` |
| `_residual_at(patch, x0, y0, z0, cfg)` | 183 | Ground-truth `\|f\|` surface-membership check |
| `_project_to_slice(patch, x0, y0, z_val, cfg)` | 222 | Gauss–Newton re-projection onto the slice curve (3D analogue of `_project_to_curve`) |
| `patch_direction(patch, x0, y0, z0, cfg)` | 245 | The `(f_y, −f_x)` in-slice patch direction |
| `build_face_tracker(patch, x0, y0, z0, cfg)` | 263 | Per-anchor tracker with **literal** `x0,y0,a,b` coefficients, `compile = :none` (Option B, benchmarked — `FaceTracking.jl:27-35`) |
| `track_dense_path(...)` | 290 | Walks a caller-supplied sequence of z-targets, reusing `_track_path_segment!`'s retry logic between each consecutive pair, guaranteeing ≥1 accepted point per hop |
| `_sweep_hop!(...)` | 394 | Single hop with both quality gates |
| `_sweep_direction(...)` | 479 | Re-anchoring decision |
| `sweep_face_bidirectional(...)` | 524 | Sweeps up and down from `z_mid`, each direction building its own tracker |
| `track_face(...)` | 574 | Top-level: one mid-slice edge → one `Face` mesh |

**Adaptive re-anchoring** (`FaceTracking.jl:37-64`) is the key robustness mechanism and uses **two complementary gates, both required**:

1. A dimensionless **cosine-similarity** check between the anchor's fixed gradient direction and the current local gradient, threshold `cfg.patch_transversality_cos_tol`. A fixed patch line is only guaranteed to keep intersecting the level curve if the gradient is radially symmetric about the shrink axis — true for a sphere, false for an ellipsoid. When the fixed line loses transversal intersection entirely, *no amount of bisection can fix it* because the system genuinely has no nearby solution.
2. A **residual-based bisection** gate reusing `cfg.critical_point_tol`, which recursively halves a hop whose landing fails `f ≈ 0`. Needed because right at a z-critical target the level curve degenerates to a point and **no** patch stays transversal — only shrinking the step converges to the true limit.

Each of the two sweep directions builds its own tracker from the same anchor and re-anchors independently; re-anchoring in one direction must never affect the other (`FaceTracking.jl:66-72`).

---

### 1.7 `src/Projection.jl` — generic SO(3) projections

**Why it exists** (`Projection.jl:10-17`): the auto-augmented critical system `{f, ∂f/∂x, ∂f/∂y}` contains an **identically-zero polynomial** whenever the surface is independent of an augmenting axis (e.g. `z − x²`, or any y-independent surface), and HomotopyContinuation's polyhedral start system then throws `OverflowError: Cannot compute a start system` from deep inside `solve`. A generic rotation makes that configuration measure-zero.

**Why SO(3) and not O(3)** (`Projection.jl:19-24`): `weld_mesh` aligns triangle windings with `+∇f` in the chart. For an orthogonal map, `cross(Qa, Qb) = det(Q)·Q·cross(a,b)` while gradients map as `∇f_world = Q·∇f_chart`, so a reflection (det = −1) would silently flip every mapped-back normal against the winding convention. Reflections are **rejected, not auto-fixed** (`Projection.jl:94-98`).

| Function | Line | Notes |
|---|---|---|
| `random_orthogonal_matrix(::Type{T}, n; rng)` | 41 | Haar-uniform SO(n) via QR of a Gaussian with `diag(R)` sign correction (Mezzadri, *Notices AMS* 54, 2007 — cited at `Projection.jl:44`), then column negation if det < 0. Distributionally equivalent to Bertini's Stewart Householder construction `make_matrix_random_real_mp` (`Projection.jl:37-38`) |
| `_resolve_projection(projection, rng, cfg)` | 80 | `:random` draws fresh; a matrix must be 3×3, satisfy `norm(Q'Q − I) ≤ cfg.projection_orthonormality_tol`, and have det > 0 |
| `_rotate_system(F, Q)` | 115 | Builds `F'(x') := F(Q x')` by simultaneous substitution |
| `_chart_config(cfg, Q)` | 131 | Replaces `bbox_*` with the enclosing AABB of the rotated world box (`Q'·corners`). **The chart working domain is therefore strictly larger than the world box at the same cfg** (`Projection.jl:126-129`) |
| `_verify_projection_ok(F_chart, cfg)` | 172 | BertiniReal's `verify_projection_ok` analogue |
| `_map_to_world(vertices, edges, faces, mesh, Q)` | 212 | Maps everything back via `p = Q p'`; triangles unchanged since det = +1 |

`_verify_projection_ok` deserves emphasis for the presentation because its scope is narrower than the name suggests. It checks **exactly the crash class**: whether either *augmenting* partial (∂f'/∂x, ∂f'/∂y) vanishes identically, tested by evaluation at 3 fixed deterministic complex probe points (`_PROJECTION_PROBES`, `Projection.jl:147-151`) against a relative threshold of `cfg.jacobian_rank_tol`. Explicitly **not** checked (`Projection.jl:165-168`): ∂f'/∂z ≡ 0 (a chart-z-independent surface like a cylinder has an empty critical system and decomposes fine), and positive-dimensional critical loci. It is a crash guard, **not a transversality test** — see §4.7.

---

### 1.8 `src/SurfaceDecomposition.jl` — the 3D orchestrator

`decompose_3d_surface(F, cfg; projection = nothing, rng = Random.default_rng(), incidence = false, deflate = false)` — `SurfaceDecomposition.jl:1798`.
Returns `(vertices, edges, faces, mesh)`, or a 5-tuple additionally returning a `SurfaceIncidence` when `incidence = true` (`SurfaceDecomposition.jl:2041`).

**Hard input contract** (`SurfaceDecomposition.jl:1806-1810`, mirrored at `:49-50`): exactly 1 equation in exactly 3 variables, ordered `[x, y, z]` with **z last**.

#### Execution order

| Stage | Lines | Mathematics |
|---|---|---|
| Projection wrapper (Phase 8) | 1812–1831 | Resolve Q, build `F_chart = F(Qx')`, verify, build chart cfg, **recurse with `projection = nothing`**, map results back |
| Slab construction | 1833–1835 | `_slab_bounds`: sorted `bbox_z` endpoints + in-range critical z, merged below `min_slab_width` |
| Patch prep | 1839 | `build_patch_system(F)` |
| Registry + incidence init | 1841–1883 | `VertexRegistry{T}(cfg.vertex_match_tol)` when `incidence = true`; `_surface_critical_vertices` for fold anchors |
| Per-slab loop | 1923–2017 | mid-slice → renumber → append → crit-slice prep → `track_face` per edge → landing assignment → continuity check |
| Weld | 2019–2041 | `weld_mesh`, with or without incidence-based stitching |

#### Key functions

**`compute_critical_z_slices(F, cfg) → Vector{T}`** — `:48`
Calls `compute_critical_points` on the raw surface (which auto-builds `{f, ∂f/∂x, ∂f/∂y}`), extracts `Re(z)`, clusters with `cluster_scalars(…, vertex_match_tol)`. These are the z-values where the surface fails to be a local graph over the (x,y) plane.

**`slice_at_z(F, z_val, cfg; deflate = false) → (vertices_3d, edges_3d)`** — `:92`
Substitutes `Float64(z_val)`, runs `decompose_1d_curve` on the resulting 2D system, appends z. **Scoping note** (`:78-82`): a vertex flagged `Singular` here was flagged against the **2D slice curve**, not the 3-variable surface — deflating against `F` instead would ask a different, higher-dimensional question.

**`_robust_slice_at_z(F, patch, z_bottom, z_top, cfg; deflate) → (vertices_3d, edges_3d, z_mid)`** — `:286`
Choosing z_mid naively as the slab midpoint can land exactly on a degenerate configuration. The Taubin heart at `z_mid = 0` is the canonical case: the equation becomes an exact cube (`:137-138`). What then happens is documented very precisely (`:139-146`):

> "decompose_1d_curve has no way to know this ahead of time; what actually happens is that its true critical points get correctly classified Singular … but connect_the_dots! then cannot track paths back onto those Singular vertices within vertex_match_tol, so the Topology._resolve_endpoint fallback fabricates new Artificial vertices wherever the paths actually (wrongly) landed -- silently producing edges that do not follow the true curve at all."

Two gates guard against this:
1. **Topology gate** — reject a candidate z_mid producing an `Artificial :endpoint_fallback` vertex co-occurring with a `Singular` vertex (unless reference slices also carry `Singular`, which indicates a genuine transversal singular curve rather than a bad slice).
2. **Gradient gate** — compare the minimum per-edge `patch_direction` magnitude at the candidate against a cross-z reference at the slab quarter-points (0.25 and 0.75 fractions, `:218-219`); suspect if below `z_mid_gradient_ratio_tol × ref_max`. This catches slices that are topologically clean but numerically ill-conditioned. Raw residual and `|∇_xy f|/|f_z|` thresholds were tried and rejected (`Config.jl:93-94`).

Retries perturb z_mid by `z_mid_retry_frac × (z_top − z_bottom)` with alternating offsets, capped at 45% of the slab (`:175-176`) and `max_z_mid_retries` attempts. If every attempt remains suspect it **throws** — "never a silent fallback to a slice already known to be wrong" (`:167-169`).

**`_slab_bounds(F, cfg) → Vector{T}`** — `:411` — bbox endpoints + filtered critical z, merged via `cluster_scalars(…, min_slab_width)`, outer bounds clamped to exact bbox.

#### Incidence machinery (Phases 9a/9b/9c)

**`CritSlice{T}`** — `:447` — `boundary_index`, `z`, `vertices`, `edges`, `is_degenerate`. The 1D decomposition at an interior critical-z boundary; `is_degenerate = true` when it has no edges (a fold or point boundary — e.g. sphere/ellipsoid poles).

**`ColumnLanding{T}`** — `:472` — `kind::Symbol` ∈ `{:edge, :crit_slice_vertex, :critical_point, :bbox, :none}`, `id`, `dist`. Distances are **recorded, never thresholded** in Phase 9a: assignment is nearest-wins, "deliberately not guessed ahead of the data" (`:466-470`).

**`SurfaceIncidence{T}`** — `:526` — 12 fields spanning crit-slices, critical vertices, per-face bottom/top edge sets and anchors, the edge→faces map, per-side column landings, and continuity diagnostics.

**`_landing_confidence(point, landing, edge_spacing, scale_ref, patch, cfg) → Bool`** — `:722`
A landing is snap-eligible iff `dist ≤ incidence_snap_tol_ratio × local_scale` **and**, for `:edge` kind, its own `|f|` residual ≤ `critical_point_tol`. `local_scale` is the assigned edge's median chord spacing or the face's own generating-column spacing. `:bbox` and `:none` are never confident.

**`_monotone_snap_targets(points, targets) → Vector{Int}`** — `:951` — dynamic program minimizing total distance subject to monotone target indices.

**`_snap_boundary_points!(...)`** — `:1032` — two passes: snap confident vertex/anchor landings to exact cell coordinates and pool edge columns; then a joint monotone snap per crit-slice edge across all faces, so the subsequent clustering merges them and closes the cracks.

**`_split_t_junctions(triangles, edge_polylines)`** — `:1141` — fan-splits triangles whose edge spans intermediate polyline vertices, up to 4 fixed-point passes.

**`_bridge_ribbon(P, Q)`** — `:1271` — O(n·m) DP minimum-cost triangulation between two polylines (Phase 9c lofting).

**`weld_mesh(faces, patch, cfg; incidence = nothing) → GeometryBasics.Mesh`** — `:1581`
Collect points → optional snap/loft → `cluster_points_indexed(…, vertex_match_tol)` → remap triangles → append loft triangles → winding correction against `+∇f` via `_gradient_at` → optional T-junction splitting.

**It returns a bare `GeometryBasics.Mesh` only — no watertightness flag, no naked-edge count.** The audit is a separate utility, `_naked_mesh_edges(mesh)` (`:1697`), which returns the undirected edges belonging to exactly one triangle.

---

### 1.9 `src/Clustering.jl` and `src/Visuals.jl`

| Function | Line | Role |
|---|---|---|
| `cluster_vertices(vertices, tol)` | `Clustering.jl:40` | Merge coincident `NativeVertex`es within `tol` |
| `_merge_cluster(vertices, member_idxs)` | `Clustering.jl:78` | Representative selection within a cluster |
| `merge_metadata(metadatas, ids)` | `Clustering.jl:124` | Combine per-vertex metadata across a merge |
| `cluster_scalars(xs, tol)` | `Clustering.jl:196` | 1D clustering — used for distinct x-slots and critical z-values |
| `cluster_points_indexed(points, tol)` | `Clustering.jl:232` | Point clustering returning an index map — the welding primitive |
| `VertexRegistry{T}` | `Clustering.jl:376` | Mutable global vertex-id unifier |
| `register!(reg, candidate)` | `Clustering.jl:434` | Register-or-match into the global namespace |

`Visuals.jl` provides `plot_curve_decomposition` (`:131`), two `plot_surface_decomposition` methods — one on `GeometryBasics.Mesh` (`:249`), one on `Vector{Face{T}}` (`:338`) — and `interactive_3d_viewer` (`:394`), with per-`VertexType` color/marker tables at `:59` and `:66`.

---

## 2. Configuration & Tolerances (`src/Config.jl`)

`HomotopyConfig{T<:AbstractFloat}` is a `@with_kw` struct (`Config.jl:50`). Every numeric default is written `T(...)` so `HomotopyConfig{BigFloat}()` gets full-precision literals rather than Float64-rounded ones (`Config.jl:20-21`).

The file opens with a warning that is worth quoting in a talk, because it explains why there are three separate "closeness" tolerances rather than one (`Config.jl:51-55`):

> `critical_point_tol`: solution quality at compute_critical_points ("did the solver converge?").
> `boundary_point_tol`: geometric containment at intersect_bounding_object ("is this on the box?").
> `vertex_match_tol`: deduplication at ConnectTheDots!/Merge ("are these the same vertex?").
> Mixing these three is a common source of subtle bugs — each bounds "close enough" at a different pipeline stage (residual vs containment vs identity).

### 2.1 Tolerances

| Field | Default | Line | Where applied |
|---|---|---|---|
| `critical_point_tol` | `1e-6` | 56 | Nearly-real filter in `compute_critical_points` (`Solver.jl:766,772`) and `compute_midslice` (`Topology.jl:82,85`); Newton stopping criterion in `_newton_polish` (`Solver.jl:668,673`) and `_project_to_curve` (`Topology.jl:253`); `poor_acc_tol` for tracker bisection (`Topology.jl:202`, `PathTracking.jl:162`); residual acceptance against `F_original` in `verify_isosingular_dimension` (`Solver.jl:432,459`); `:edge`-landing residual gate in `_landing_confidence` (`SurfaceDecomposition.jl:700,703,711,736`); the residual bisection gate in `FaceTracking`'s re-anchoring (`FaceTracking.jl:52-53`) |
| `boundary_point_tol` | `1e-5` | 57 | `intersect_bounding_object` only — nearly-real filter (`Solver.jl:845,878`) and the inside-box slack `lo − btol ≤ val ≤ hi + btol` (`Solver.jl:884`) |
| `vertex_match_tol` | `1e-4` | 58 | Vertex dedup in `cluster_vertices` (`Solver.jl:790,921`); distinct-x clustering (`Topology.jl:382`); endpoint match radius in `_resolve_endpoint` (`Topology.jl:151`); `min_width` bisection floor for tracking (`Topology.jl:203`, `PathTracking.jl:166`); critical-z clustering (`SurfaceDecomposition.jl:59`); coincidence test in `_cells_adjacent` (`:746`); mesh welding in `cluster_points_indexed` (`:1615`); `VertexRegistry` tolerance (`:1856,1866`) |
| `jacobian_rank_tol` | `1e-8` | 59 | Singular-value cutoff for numerical rank: `count(>(cfg.jacobian_rank_tol), svals)` (`Solver.jl:55`); propagates into `estimate_corank`, `deflate_once`, `_deflation_applicable`; `rank_tol` during tracking (`Topology.jl:200`, `PathTracking.jl:68`); relative probe threshold in `_verify_projection_ok` (`Projection.jl:185,187`) |
| `singular_value_threshold` | `1e-6` | 60 | Second, independent singularity test `min(σ) < threshold` in `_classify_vertex_type` (`Solver.jl:78`) and `is_near_singular` (`PathTracking.jl:69`) |
| `path_tracker_precision` | `1e-10` | 61 | `TrackerOptions(min_step_size)` in `build_tracker` (`PathTracking.jl:89`) — the only consumer |
| `patch_transversality_cos_tol` | `0.9` | 67 | `FaceTracking._sweep_direction` re-anchoring. Dimensionless cosine bound ≈ 26° rotation. Comment at `Config.jl:62-66` records that reusing `singular_value_threshold` for this "fired one hop too late on asymmetric-ellipsoid regression because that scale tracks equation magnitude, not patch skew" |
| `projection_orthonormality_tol` | `1e-8` | 77 | Acceptance threshold for `norm(Q'Q − I)` in `_resolve_projection` (`Projection.jl:89-90`). Deliberately its own field, not a reuse of `jacobian_rank_tol` despite the same default magnitude — a matrix-identity defect is a different physical quantity from a Jacobian singular-value cutoff (`Config.jl:68-76`) |

### 2.2 Knobs

| Field | Default | Line | Where applied |
|---|---|---|---|
| `max_path_steps` | `1000` | 78 | Per-direction step budget for `track_path` (`PathTracking.jl:208`) |
| `bbox_x`, `bbox_y`, `bbox_z` | `(-4.0, 4.0)` each | 79–81 | Bound substitution in `intersect_bounding_object` (`Solver.jl:844`); slab outer bounds (`SurfaceDecomposition.jl:412-413`); rotated to a chart AABB by `_chart_config` (`Projection.jl:131`) |
| `edge_sample_density` | `50` | 82 | Output point count in `sample_edge` (`Topology.jl:302`) |
| `midslice_sample_density` | `100` | 83 | Face grid rows: `n_z = 2·density + 1` (`SurfaceDecomposition.jl:1064,1221,1392,1626,1999`) |
| `z_mid_retry_frac` | `0.01` | 88 | Perturbation as a **fraction of slab width**, not an absolute offset (`SurfaceDecomposition.jl:378,388`) |
| `max_z_mid_retries` | `8` | 90 | Attempt cap before `_robust_slice_at_z` throws (`:375,387`) |
| `z_mid_gradient_ratio_tol` | `0.01` | 95 | Gradient-degeneracy gate (`:354,358,389`) |
| `min_slab_width` | `1e-3` | 105 | Resolution floor merging critical z-values in `_slab_bounds`; also the fold-anchor z-proximity window in `_assign_landings` (`:611,625`) |
| `incidence_snap_tol_ratio` | `1.5` | 120 | Phase 9b confidence gate (`:699,733,829`) |
| `isosingular_verify_retries` | `20` | 131 | Fresh-hyperplane retry budget in `verify_isosingular_dimension` (`Solver.jl:436,469`) |
| `max_deflations` | `10` | 138 | Round budget in `resolve_isosingular_dimension` (`Solver.jl:625`) |

### 2.3 Calibration provenance (the reason these numbers are what they are)

Several defaults carry unusually detailed empirical justification in-source. These are good talk material because they show the tolerances were measured, not guessed.

- **`min_slab_width = 1e-3`** (`Config.jl:96-105`): path endpoints landing on a point singularity of multiplicity *m* carry `~accuracy^(1/m)` scatter in z — **~2e-4 observed**, beyond `vertex_match_tol = 1e-4`. That scatter mints ~2e-4-wide slabs centered *on* a singular point that no z_mid choice can slice. The default sits **5× above the observed scatter and 65× below the smallest genuine fixture gap** (fixed-axis Taubin: 0.065), so all fixed-axis fixture bounds are unchanged. Explicitly a resolution limit: genuinely distinct critical values closer than this are **not** separated.
- **`incidence_snap_tol_ratio = 1.5`** (`Config.jl:106-120`): calibrated on the fixed-axis Taubin fixture. Honest `dist/spacing` ratios topped out at **0.828** across two structurally different edge-type boundaries (z = 1.0: median 0.094 / max 0.828; z = 1.0648: median 0.389 / max 0.825). The default carries ~1.8× margin above that ceiling. It is a single dimensionless constant rather than a per-boundary absolute cutoff because raw distances scale with local sample spacing — **12× spread within one boundary alone**.
- **`isosingular_verify_retries = 20`** (`Config.jl:121-131`): matches Bertini1's `isosingularDimTest maxIts` (`isosingular.c:336`). A safety margin well above the observed range, **not** a minimum requirement — and it explicitly does *not* gate `Inconclusive`, which fires on the first residual-failing success regardless of remaining budget.
- **`max_deflations = 10`** (`Config.jl:132-138`) and **`patch_transversality_cos_tol = 0.9`** are likewise documented as safety margins rather than empirically-required minima.

---

## 3. Fixtures & Test Suite

### 3.1 Suite structure

`test/runtests.jl` (41 lines) always includes 11 files (`runtests.jl:12-21,40`): `test_docstring_rendering`, `test_types`, `test_vertex_registry`, `test_solver`, `test_topology`, `test_pathtracking`, `test_surfacedecomposition`, `test_projection`, `test_incidence`, `test_visuals`, `test_isosingular_deflation`.

The slow gate is a single env var (`runtests.jl:23-27`):

```julia
if get(ENV, "HOMOTOPYGETSREAL_RUN_SLOW_TESTS", "0") == "1"
    include("test_taubin.jl")
else
    @info "Skipping test_taubin.jl (set HOMOTOPYGETSREAL_RUN_SLOW_TESTS=1 to run)"
end
```

Three further files gate *portions* internally on the same variable: `test_isosingular_deflation.jl:95-214` (historical curves + Taubin), `test_projection.jl:116-183` (rotated Taubin), `test_visuals.jl:83-105` (Taubin visuals).

**There are zero `@test_broken`, zero `@test_skip`, and zero commented-out `@test` calls in the suite.** Known-weak areas are instead handled by (a) slow-gating, (b) loose bounds with an explanatory comment, or (c) soft `@info` rather than a hard assertion.

### 3.2 Fixture-by-fixture

#### Nodal cubic — `f = y² − x³ − x²`
Defined at `test_topology.jl:8`, `test_pathtracking.jl:8`, `test_solver.jl:9`. The workhorse 1D fixture: it has one genuine node at the origin and one smooth x-critical point.

Topology assertions (`test_topology.jl:33-42`) encode the expected decomposition exactly:

```julia
@test length(vertices) == 4
@test count(v -> v.v_type == Critical, vertices) == 1
@test count(v -> v.v_type == Singular, vertices) == 1
@test count(v -> v.v_type == Boundary, vertices) == 2
```
plus coordinate assertions: `Critical` at `(-1, 0)` and `Singular` at `(0, 0)`, both `atol = 1e-3`, and `Boundary` at `y = ±4`.

Edge structure (`test_topology.jl:54-64`): `@test length(edges) == 4`, connectivity broken down as 2 critical–singular, 1 singular–boundary⁺, 1 singular–boundary⁻, 0 boundary–boundary, and `@test all(e -> e.is_singular, edges)`. Distinct x-slots: `@test length(distinct_xs) == 3` (`:76`). Sampling decoupling is tested by re-running at density 200 (`:93`). BigFloat T-genericity reproduces the same vertex-type counts (`:118-120`).

Solver-level (`test_solver.jl:48-51,80-82`): 2 critical points with the singular one at the origin and the critical one at x = −1; 2 boundary points at `|y| = 4`. Tightening the bbox moves them (`:97-98`) — a wiring test proving `bbox_*` is actually read.

#### Unit sphere — `f = x² + y² + z² − 1`
`test_surfacedecomposition.jl:49`, `test_incidence.jl:10`, `test_projection.jl:31`, `test_visuals.jl:14`. The reference smooth-surface fixture.

```julia
@test length(z_crits) == 2                                    # :64
@test isapprox(sort(z_crits), [-1.0, 1.0]; atol = 1e-6)       # :65
```
Equator slice: 2 vertices at `|x| = 1, y = 0`, 2 edges, exact `z == 0.0` on all sampled points (`:75-79, :91`). Patch variable reordering is asserted literally: `@test string.(patch.F_for_tracking.variables) == string.([z, x, y])` (`:101`). Dense sweep points stay on the sphere to `atol = 1e-4` (`:137-138`), and — notably — the sphere needs **zero** re-anchoring: `@test n_rebuilds_down == 0`, `@test n_rebuilds_up == 0` (`:152-153`). This is the fixture that would hide the ellipsoid bug, which is exactly why the ellipsoid exists. Face mesh residuals `atol = 1e-3` (`:166`), welded mesh `atol = 1e-2` (`:189`), outward winding (`:203`), and `@test isapprox(z_used, z_naive; atol = 1e-12)` confirming the naive midpoint is accepted unmodified (`:218`).

#### Ellipsoid — `f = x² + 4y² + 9z² − 1`
`test_surfacedecomposition.jl:232`, `test_incidence.jl:55`, `test_visuals.jl:76`. Asymmetric by construction so an x/y/z swap cannot pass silently, and it is the fixture that forced adaptive re-anchoring into existence.

```julia
@test length(z_crits_ell) == 2                                     # :238
@test isapprox(sort(z_crits_ell), [-1/3, 1/3]; atol = 1e-6)        # :239
@test all(<=(1e-4), ell_residuals)                                 # :263
```
where `ell_residuals = [abs(p[1]^2 + 4p[2]^2 + 9p[3]^2 - 1.0) for p in ell_mesh_pts]` (`:252`) — a direct pointwise residual check against the defining polynomial at every welded mesh vertex. The bound was tightened to 1e-4 after re-anchoring landed empirically below 1e-6 (`:245-251`). A cosine-tolerance sweep asserts monotone behavior: `@test issorted(rebuild_counts)` and `@test results[1].max_resid >= results[end].max_resid - 1e-8` (`:340,342`) — i.e. a looser transversality tolerance rebuilds less and is never more accurate.

#### Taubin heart — `f = (x² + (1.2y)² + z² − 1)³ − x²z³ − 0.1(1.2y)²z³`
`test_taubin.jl:150`, `test_isosingular_deflation.jl:145`, `test_projection.jl:125`, `test_visuals.jl:86`. Degree 9, slow-gated, and the fixture that drives most of the robustness machinery. It has a positive-dimensional singular curve `{x² + 1.44y² = 1, z = 0}` (`SurfaceDecomposition.jl:251`) plus cusp and tip point singularities.

```julia
@test length(z_crits) == 4                                                                   # :195
@test isapprox(sort(z_crits), [-1.0, 1.0, 1.0647678179140714, 1.2366591700121616]; atol=1e-6) # :196
@test all(s -> s.max_resid < 1e-4, slab_stats)                                               # :266
@test all(<=(1e-4), heart_residuals)                                                         # :305
@test n_degenerate == 0                                                                      # :318
```

The z_mid fix is regression-tested from both sides. Before the fix, `z_mid = 0` produced max `|f|` up to **1.64** (`:256-264`); a perturbation to 0.02 is still insufficient — `@test resid_002 > 1e-2` (`:426`) — while the accepted fix gives `@test n_fallback_fixed == 0`, 2 vertices, 2 edges, `@test resid_fixed < 1e-4` (`:442-445`).

Crit-slice census under `incidence = true`: 4 slices with edge counts 0/4/12/0 and matching degeneracy flags (`:485-492`).

#### Whitney umbrella — `f = x² − y²z`
`test_isosingular_deflation.jl:23,63`, `test_solver.jl:327`. This is the **published-ground-truth** fixture, reproducing Hauenstein–Wampler's deflation sequences on the bare, unsliced equation:

```julia
@test seq_tip == [3, 2, 0]                              # test_isosingular_deflation.jl:82
@test HomotopyGetsReal._corank_plateau_hint(seq_tip)    # :83
@test seq_handle == [3, 1, 1]                           # :92
```
The tip `(0,0,0)` is an isolated point that resolves in two rounds; the handle `(0,0,1)` sits on the 1-dimensional singular locus and produces a genuine plateau. The handle sequence is hard-capped at 2 rounds because **round 5 of this exact case needs ~5.3e9 symbolic minors and never terminates** — it had to be killed (`test_isosingular_deflation.jl:86-90`).

Verdict-level assertions (`test_solver.jl`): handle verify returns `Verified` with `1 ≤ attempts ≤ isosingular_verify_retries` and `inconsistent_count == 0` (`:344-356`); a deliberately wrong `d` returns `Inconclusive` (`:366`); round-1 returns `NotTerminal` and round-2 `Verified` with `attempts == 0` (`:371-378`); full resolve gives `@test r_handle_full.corank_sequence == [3, 1, 1]` and `@test r_handle_full.rounds == 2` (`:459-463`).

Both failure verdicts are exercised with constructed cases: `Exhausted` via `max_deflations = 1` (`:550-552`), and `Ambiguous` via a shifted umbrella `f = x² − y²z − 1e-3` (`test_solver.jl:563`) whose Jacobian is identical to the true umbrella but whose residual against the original system fails — `@test r_ambiguous.verdict == Ambiguous`, `@test r_ambiguous.corank_sequence == [3, 1, 1]` (`:574-577`).

#### Node and cusp toys — `f = y² − x²`, `f = y² − x³`
`test_solver.jl:208-209`. Corank at origin vs smooth point (`:223-226`); deflated equation counts (`:303-308`); resolve round counts — node resolves in 1 round with equation counts `[3]`, cusp in 2 rounds with `[3, 6]` (`:440-443`).

#### Astroid — `f = (x² + y² − 1)³ + 27x²y²`
`test_isosingular_deflation.jl:109`, slow-gated. End-to-end `decompose_1d_curve(...; deflate = true)`, asserting every flagged `Singular` vertex resolves:

```julia
@test verdict == Resolved                               # test_isosingular_deflation.jl:131
```

This loop is the source of the documented ±1 test-count variance — it fires **one assertion per resulting `Singular` vertex**, and that count depends on what a given live run happens to flag. Paper artifacts record the astroid as 4 vertices / 4 edges, all 4 `Singular` (`results.json:359-369`).

Two companion historical curves share the same testset: `f2 = (x³ − xy² + y + 1)²(x² + y² − 1) + y² − 5` (`:110`) and a degree-8 curve `f3 = x⁸ − 28x⁶y² + 70x⁴y⁴ − 28x²y⁶ + y⁸ + 15x⁴y² − 15x²y⁴` (`:111`), the latter run twice, at default bbox and `bbox = (-8, 8)`. This testset alone was measured at **6m24.8s**, which is why it was moved behind the slow gate (`:103-107`).

#### Parabolic sheet — `f = z − x²`
`test_projection.jl:52`, `test_incidence.jl:68`. The crash-class fixture: bare decomposition **throws** (`test_projection.jl:81`), works under `projection = :random` with world residual `@test res_p < 1e-3` (`:90`), and has no interior crit-slices with all landings `:bbox` (`test_incidence.jl:71,74`).

#### Cylinder — `f = x² + y² − 1`
`test_projection.jl:70`. A scoping test only: `@test HomotopyGetsReal._verify_projection_ok(F_cyl, cfg) === nothing` (`:71`) — proving the guard does *not* reject chart-z-independent surfaces.

### 3.3 Torus — a scope correction

**The torus is not a test fixture.** `rg -n 'torus|Torus' test/` returns no matches. It exists only as a paper artifact, generated by `paper_artifacts/torus_example.jl` and recorded in `results.json:501-536`. Nothing in `test/` asserts anything about it, so its numbers are reproducible measurements rather than enforced invariants — an important distinction if the torus appears in a talk alongside the tested fixtures.

### 3.4 Measured results (paper artifacts, `paper_artifacts/data/results.json`)

These are recorded measurements, not test assertions.

| Fixture | Config | V (by type) | E | F | Mesh V / △ | mean \|f\| | p99 \|f\| | max \|f\| | wall time |
|---|---|---|---|---|---|---|---|---|---|
| Sphere | default (50/100) | 2 Critical | 2 | 2 | 19504 / 39004 | 2.85e-8 | 7.05e-8 | 7.58e-7 | 0.96 s |
| Ellipsoid | default (50/100) | 2 Critical | 2 | 2 | 19500 / 39120 | 2.70e-7 | 8.99e-8 | 1.12e-4 | 1.25 s |
| Taubin heart | paper (8/8) | 10 Critical + 4 Artificial | 14 | 14 | 1638 / 3124 | 1.03e-7 | 2.00e-6 | 2.42e-6 | 34.2 s |
| Torus | `:random`, seed 42 (50/100) | 8 Critical | 8 | 8 | 78340 / 160381 | 2.20e-6 | 1.00e-5 | 1.78e-5 | 54.6 s |

Source lines: sphere `results.json:21-42`, ellipsoid `:109-130`, Taubin `:201-236`, torus `:510-535`.

Two details worth carrying into a talk:

- **The Taubin heart's default surface decomposition reports zero `Singular` vertices** (`results.json:205-210`) — 10 `Critical` and 4 `Artificial`. The cusps only appear when the incidence diagnostic is enabled, which yields 14 `Critical` + 2 `Singular` in the overlay (`results.json:439-444`). The reason is the scoping note at `SurfaceDecomposition.jl:78-82`: `Singular` in the default path is classified against the 2D slice curve, not the 3-variable surface.
- **The Taubin `[-1, 1]` slab genuinely exercises the retry ladder** (`results.json:222-233`): naive `z_mid = 0` failed both gates with `endpoint_fallback = 4` and `Singular = 2`; after **5 retries** it accepted `z_mid = 0.06`, taking 6.30 s.

Curve artifacts (`results.json:358-409`): astroid 4V/4E all `Singular`; cusp 3V/2E (2 `Boundary`, 1 `Singular`); nodal cubic 4V/4E (2 `Boundary`, 1 `Critical`, 1 `Singular`); and — flagged in the artifact itself — the `node` figure is **illustrative, hand-placed**, because `decompose_1d_curve` finds 0 critical vertices there (`results.json:408`).

---

## 4. Weak Spots & True Code State

This section is the honest inventory. Most of it comes from `docs/DESIGN_NOTES.md`, which is unusually candid and contains an explicit backlog and several retractions.

### 4.1 Multiplicity ≥ 2 — the single biggest correctness risk

`docs/DESIGN_NOTES.md:75-178` consolidates a recurring failure pattern. In every case, the critical-point-finding augmented system has a solution of **multiplicity ≥ 2 at exactly the point being searched for**, and HomotopyContinuation.jl's default polyhedral `solve` finds it unreliably.

| Fixture | Failure mode | Source |
|---|---|---|
| Node curve `y² − x²` | `compute_critical_points` on `[y² − x², 2y]` returns **zero** solutions; both tracked paths report `return_code = excess_solution`, neither converging to the genuine double root at the origin | `DESIGN_NOTES.md:75-178` |
| Astroid's 4 cusps | Across 5 independent runs, **2 came back with 1 or 2 cusps missing entirely**, replaced by `Artificial` fallback vertices — not coordinate jitter, actual feature loss | same |
| Whitney umbrella apex | `compute_critical_z_slices` returned `Float64[]` — no critical z found at all | same |
| Cone apex | `compute_critical_z_slices` returned `Float64[]`; `decompose_3d_surface` **silently returns a completely empty decomposition** (0 vertices, 0 faces) with no crash — "arguably worse to notice" | same |
| Horn torus pinch | Same empty decomposition bare; **succeeds** under `projection = :random, rng = Xoshiro(42)` | same |
| Three concurrent lines `x(x−y)(x+y)` | Triple point found but **misclassified `Artificial`, not `Singular`**; separately, the vertical line `x = 0` is structurally invisible to x-parametrized slicing and never appears in the edge graph | `DESIGN_NOTES.md:490-503, 492-495` |

The working hypothesis is explicitly labelled unverified (`DESIGN_NOTES.md:171-178`): "not yet verified against HC.jl internals: `solve`'s default settings struggle specifically when the critical-point-finding augmented system has a solution of multiplicity ≥ 2 at the point being searched for."

**A rejected workaround worth mentioning as a methodological example.** An asymmetric `bbox_z = (-4.0, 4.3)` makes the cone produce a *non-empty* decomposition, but the result is geometrically corrupted: median residual is fine at 1.2e-7 while **p99 = 12.2 and max = 16.0**. `DESIGN_NOTES.md:118-129` records the conclusion: "This is more dangerous than the original empty-mesh failure, not a partial fix. Do not present this as a usable workaround under any framing." Net: **no known safe workaround exists for the cone.**

A related, falsified prediction is flagged rather than resolved (`DESIGN_NOTES.md:160-169`): the geometric argument that "the singular point sits at the origin, fixed under any rotation, so `projection = :random` shouldn't help" correctly predicted the cone's failure but **not** the horn torus's success under identical treatment. The divergence was not investigated.

**Downstream consequence.** There is an uncaught-exception path when `compute_critical_z_slices` finds zero critical z-values *and* the naive bbox midpoint is itself degenerate: `slice_at_z(F, 0.0, cfg)` throws `OverflowError: Cannot compute a start system`, and `_robust_slice_at_z`'s retry ladder has **no `try`/`catch` around the tracking calls**, so the exception propagates uncaught (`DESIGN_NOTES.md:54-73`).

**A separate class, do not conflate** (`DESIGN_NOTES.md:180-191`): the multiplicity issue above is a *detection-reliability bug*. The absence of first-class singular-curve decomposition (§4.6) is a *missing capability*.

**Also documented:** smooth critical points correctly located but **misclassified `Singular`**, with downstream edge/mesh construction bypassing them entirely and wiring to off-curve `Artificial :endpoint_fallback` vertices instead — residual up to ~1.0 on the surface case (`squircle_quartic`, `quartic_superellipsoid`; `DESIGN_NOTES.md:461-477`).

### 4.2 Polyhedral solver and BigFloat boundaries

- **Path tracking is Float64-only, unconditionally.** HC's `NewtonCache`/`NewtonResult` are hardcoded to `ComplexF64`; "there is no way to path-track directly in BigFloat" (`Solver.jl:7-12`). BigFloat is a *post-processing* precision, applied by `_newton_polish` to landing points only.
- **`Complex{BigFloat}` has no polyhedral start system.** `ToricHomotopy` has no method for it, so a system built from literal BigFloat coefficients throws a `MethodError` (`Solver.jl:859-865`). The uniform workaround is Float64 substitution before `solve`, applied identically at `Solver.jl:868-872`, `Topology.jl:71-77`, and `SurfaceDecomposition.jl:99`.
- **HC's convenience wrappers silently default to 53 bits.** `jacobian(F, x)` and `F(x)` ignore `x`'s element type, so the codebase uses low-level `evaluate(...; bits = precision(T))` throughout (`Solver.jl:13-21`).
- **BigFloat SVD needs `GenericLinearAlgebra.jl`** — Base's `svdvals` has no generic fallback and throws `MethodError` (`Solver.jl:22-27`).
- **The polyhedral lifting is unseedable from HGR's API.** `HomotopyContinuation.solve()` exposes no `rng`/seed parameter; `MixedSubdivisions.jl`'s `uniform_lifting_sampler` draws from Julia's global `default_rng()` with no way to seed it through this codebase's own `rng` kwarg — that kwarg only ever reaches the projection-matrix draw (`docs/ORCHESTRATOR_BRIEFING.md:95-108`). This is the mechanism behind run-to-run mesh-count variation, and it means **decompositions are not bit-reproducible across processes even at a fixed seed**.

### 4.3 Mesh watertightness — measured, partially closed, explicitly deferred

`weld_mesh` produces a mesh but makes no watertightness guarantee. The audit utility `_naked_mesh_edges` (`SurfaceDecomposition.jl:1697`) counts undirected edges belonging to exactly one triangle; a closed surface within the bbox should have none.

The fixed-axis Taubin heart is the canonical measurement, and the numbers are **hard-asserted in the test suite**, not just documented (`test_taubin.jl:469-473`):

```julia
@test length(naked) == 188
@test count_naked_near(-1.0)    == 10
@test count_naked_near(1.0)     == 70
@test count_naked_near(1.0648)  == 84
@test count_naked_near(1.2367)  == 24
```

Progression across the three stitching mechanisms (`DESIGN_NOTES.md:1146-1162`):

| Stage | Naked edges | Per-boundary (z = −1 / 1.0 / 1.0648 / 1.2367) |
|---|---|---|
| Tolerance-only weld (`incidence = nothing`) | **188** | 10 / 70 / 84 / 24 |
| + Phase 9b incidence snap-stitching | **58** | 0 / 26 / 32 / 0 |
| + Phase 9c coordinated loft | **31–35** | fold/point boundaries fully closed |

Assertions on the stitched result are deliberately **loose bounds**, with the reason stated in the test file (`test_taubin.jl:598-601`): `@test count_stitched_near(-1.0) == 0`, `@test count_stitched_near(1.2367) == 0`, `@test length(naked_stitched) < 58`, `@test length(naked_stitched) <= 40` — "regression guard around the measured 31-35 range," with the spread attributed to cross-process HC.jl jitter.

The residual naked edges break down as **roughly one third** genuine cross-edge-junction cases, **one third** seam artifacts, and **one third undiagnosed after three separate investigation attempts** (`DESIGN_NOTES.md:1113-1120`, echoed in `test_taubin.jl:586-591`). Full closure is **DEFERRED as a future, separately-scoped phase** (`SurfaceDecomposition.jl:1567-1576`, `:1690-1691`).

Three further honest details:

- The `sample_edge` re-projection fix (§1.5) **did not reduce the naked-edge count at all** — re-measured at 188 unstitched / 32–34 stitched, within the same jitter range, not a reduction (`DESIGN_NOTES.md:1122-1137`). Two real defects with independent causes.
- Counts are **density-dependent**: at survey density (`edge = 6, midslice = 8`) Taubin gives 29 from 132 bare; at paper density (8/8) ~30–31; at full production defaults (50/100) **~56–63** (`DESIGN_NOTES.md:1189-1195`). Bigger meshes have more cracks.
- A methodological caveat: a bbox-clipped rim on an *unbounded* surface is naked **by construction** and is not a stitching failure — `hyperboloid_one_sheet` shows 196 naked edges purely from its clipped boundary (`DESIGN_NOTES.md:1235-1245`).

The torus at seed 42 records 18 naked edges, with the artifact itself flagging the number as "characteristic, not deterministic -- observed to vary run to run at this same seed (21 vs 144 across two runs)" (`results.json:515-516`).

Framing for a talk: naked edges are a defect in the **mesh realization**, not evidence of a wrong decomposition (`notes/TALK_PREP.md:645-666`).

### 4.4 No formal certification

**There is none.** HomotopyContinuation.jl ships `certify`, based on interval-arithmetic Krawczyk / Smale α-theory, and it is not used anywhere: `grep -rn "certify" src/` returns nothing (`notes/TALK_PREP.md:154-159`). Real-root filtering is purely tolerance-based against `critical_point_tol = 1e-6`.

What is and is not guaranteed (`notes/TALK_PREP.md:141-152`):

| Guarantee | Gröbner / CAD / msolve | HGR |
|---|---|---|
| Solution completeness | Yes | **No proof.** Path count is a bound; missed paths are silent |
| Certification that a point is near a true root | Not needed (exact) | **Not performed** |
| Topological correctness | CAD: yes | Inferred; checked *a posteriori* by residuals and test invariants |

Even if `certify` were wired in, it would only handle **nonsingular isolated** solutions of square systems. It could not certify (i) the *completeness* of the critical-point list, (ii) the path-tracked edges and faces, or (iii) the singular points — the very things this pipeline is about. Certifying a *decomposition* is a substantially harder open problem than certifying a root (`notes/TALK_PREP.md:483-489`).

The four mechanisms actually used as correctness evidence (per `CLAUDE.md`, all with live evidence in-repo) are: pointwise residuals against the defining polynomial; cross-validation against Hauenstein–Wampler's published deflation sequences; coordinate-level comparison against Bertini_real's real astroid output (this comparison lives in an external, gitignored auditor workspace by design); and the `Test.jl` suite as a topology-invariant encoding.

**The load-bearing conditional:** every correctness claim is conditional on the completeness of the critical-point solve. If one critical value is missed, two slabs silently merge and a topological feature disappears with **no error raised** — the cone fixture returns 0 vertices, 0 faces, and no exception (`notes/TALK_PREP.md:234-241`).

### 4.5 Test-suite fragility

- **±1 variance is expected, not a bug** (`DESIGN_NOTES.md:759-769`): the historical-curves testset fires one `@test verdict == Resolved` per resulting `Singular` vertex (`test_isosingular_deflation.jl:126-131`), and that count depends on what a given live run flags. Counts of 536, 537, and 538 all appear in repository documents.
- **The full suite has occasionally errored outright**: 1 of 3 live attempts in one session, root cause unconfirmed (`DESIGN_NOTES.md:805-809`, `CLAUDE.md:78-79`). A plausible but unconfirmed partial explanation is a native segfault (signal 11) inside GLMakie's window-creation path in `test_visuals.jl` when the machine's display has gone to sleep (`DESIGN_NOTES.md:824-839`).
- **CI `test-slow` hung twice** pre-fix, last output at `test_isosingular_deflation.jl:134`, then ~29–30 minutes of zero output until the 90-minute cap (`DESIGN_NOTES.md:341-354`). After a `compile=:none` fix the job completed, but the note is explicit that this is **n = 1** and "still not permanently fixed" (`DESIGN_NOTES.md:416-449`), and that the fix was never observed to reproduce the hang locally in the first place (`:407-414`).
- **Runtime is dominated by one file**: `test_isosingular_deflation.jl` alone is ~76% of full-suite runtime — 18m15.2s of ~24m (`DESIGN_NOTES.md:774-790`).
- **Some behaviors are deliberately not asserted**: rotated-sphere face count (`test_projection.jl:114`), and on-manifold fidelity of chord-interpolated `sample_edge` points at low density (`test_surfacedecomposition.jl:80-90`).
- **Cross-process nondeterminism reaches correctness**: 1 of 3 identical `decompose_1d_curve` runs on `folium_descartes` silently orphaned a genuine on-curve `Boundary` vertex, **with no exception** (`DESIGN_NOTES.md:479-488`).

### 4.6 Deflation is diagnostic-only

The isosingular machinery works and reproduces published ground truth, but it **does not feed back into the geometry**. `ResolveResult.F_final` is retained and "not consumed anywhere yet" (`Solver.jl:513-515`); `compute_critical_points` writes verdicts into vertex metadata and leaves coordinates untouched (`Solver.jl:693-698`). Under a generic projection, `deflate = true` "does not itself change mesh quality in this band — it is diagnostic-only… actually USING that diagnosis to improve mesh quality here… remains future work" (`SurfaceDecomposition.jl:1777-1790`).

This is described as the largest real scope gap in the project (`notes/TALK_PREP.md:330-335`): Bertini_real's whole point was to *use* the diagnosis for singular-curve decomposition. HGR has no first-class representation of a singular curve as its own decomposed object.

Two hard walls in the deflation itself: round *k* appends C(rows, m)·C(cols, m) symbolic determinants, and round 5 of the Whitney handle needs **~5.3e9 symbolic minors and never terminates** (`test_isosingular_deflation.jl:86-90`); and `_deflation_applicable` exists because a real crash was hit — `"deflate_once: minorSize=2 exceeds available rows"` on Taubin heart crit-slices (`Solver.jl:708-709`).

### 4.7 Projection and fixed-axis limitations

- **`_verify_projection_ok` is a crash guard, not a transversality test.** It checks only that the augmenting partials are not identically zero. By contrast, Bertini_real computes the determinant of the stacked Jacobian + projection rows, which *is* a genuine transversality test (`notes/TALK_PREP.md:746-780`, `Projection.jl:161-168`). The almost-sure genericity argument for a random rotation is a **mathematical** argument; the code does not verify it (`notes/TALK_PREP.md:401-408`).
- **Positive-dimensional critical loci cannot be represented at all.** At the fixed-axis torus fold `z = ±1`, `∂f/∂x = ∂f/∂y = 0` identically for every point on a circle — a genuine 1-dimensional critical locus that `compute_critical_points` "cannot represent at all." This is not logged as a backlog item because the pipeline was never designed to detect them (`DESIGN_NOTES.md:270-279`). Fixed-axis torus produced ~1400 s and a catastrophically wrong mesh; `projection = :random` is the fix (`DESIGN_NOTES.md:303-304`).
- **Mesh quality degrades near a transversal singular curve.** Under a generic projection of the Taubin heart, median world residual is 6e-9 but **~2.5% of points exceed 1e-4** (max 0.26), all confined to `|z_world| ≤ 0.14` (`DESIGN_NOTES.md:1248-1253`), asserted in `test_projection.jl:173-179` as `@test sorted_h[cld(end,2)] < 1e-6`, `@test n_bad / length(res_h) < 0.05`, `@test band <= 0.2`. Labelled a known limitation at `SurfaceDecomposition.jl:1777-1790`.
- **The chart working domain is larger than the world bbox** at the same cfg, since `_chart_config` takes the enclosing AABB of the rotated box (`Projection.jl:126-129`).
- **`continuity_ok` is a write-only diagnostic** — populated by `_check_continuity!` but read by nothing; "it does not affect any mesh actually produced today" (`SurfaceDecomposition.jl:504-508`). Promotion to a hard error is deliberately deferred until the flags are trustworthy (`:512-515`, `:822-823`).

### 4.8 Scalability and dimension limits

- **Ambient dimension is architecturally capped.** Surfaces are hypersurfaces in ℝ³, hard-enforced — `compute_critical_z_slices` throws unless exactly 3 variables and 1 equation (`SurfaceDecomposition.jl:49-53`). There is **no complete-intersection or higher-codimension path at all** (`notes/TALK_PREP.md:547-551`, `docs/BERTINIREAL_AUDIT.md:322-323`).
- **Demonstrated scale is modest**: low-degree trivariate surfaces (sphere, ellipsoid, Taubin heart degree 9, torus degree 4), against Bertini_real's reported degree-630 curve in 14 variables (`notes/TALK_PREP.md:541-543`).
- **The AABB choice scales badly**: an axis-aligned box requires a separate boundary solve per face — 2N faces in N variables. The Bertini_real manual is quoted on exactly this: "for higher dimensions, it meant doing more and more curve decompositions, with each of the 2N planes. It got messy" (`notes/TALK_PREP.md:552-558`).
- **`incidence = true` costs ~+9 s on a ~16 s Taubin decompose** (`DESIGN_NOTES.md:1164-1167`), and a duplicate `compute_critical_points` solve under `incidence = true` is a known, deliberately deferred optimization (`DESIGN_NOTES.md:250-253`).
- Capability survey across 24 fixtures: **14 clean, 7 caveats, 1 exception, 2 empty** (`notes/TALK_PREP.md:943`).

### 4.9 Documented retractions and corrections

`docs/DESIGN_NOTES.md` records several claims that were made and then withdrawn. These are worth knowing before quoting any older document.

| Claim | Correction | Source |
|---|---|---|
| `deflation_stabilized` should return `false` for `[1,1,1]` (only "ends in 0" counts) | **Retracted** — contradicted H–W Def. 5.18 and the Whitney handle data `[3,1,1,1,1]`. Now returns `true`; function renamed `_corank_plateau_hint` | `DESIGN_NOTES.md:563-572`, `Solver.jl:155-160` |
| Witness-slice construction is required before deflation | **Retracted** — `deflate_once` on the bare Whitney umbrella reproduces the published sequences exactly, no slice needed | `DESIGN_NOTES.md:583-600`, `Solver.jl:202-207` |
| "Every real case resolves in exactly 1 verify attempt" | **Retracted** after measuring `1,2,1,4,1,1,4,1,1,1` across 10 trials. The invariant is only `1 ≤ attempts ≤ isosingular_verify_retries` | `DESIGN_NOTES.md:620-630` |
| Cone asymmetric-bbox workaround | **Rejected as dangerous** (p99 = 12.2, max = 16.0) | `DESIGN_NOTES.md:118-129` |
| Rotation can't help the cone/horn torus | **Falsified** for horn torus; divergence uninvestigated | `DESIGN_NOTES.md:160-169` |
| `compile=:mixed` fix resolved the CI hang | **Caveated** — never reproduced locally; CI success is n = 1 | `DESIGN_NOTES.md:407-414, 416-449` |
| GLMakie confirmed headless | **Corrected** — conditional on display being awake; display-sleep segfault documented | `DESIGN_NOTES.md:1302-1314` |
| Griffis-Duffy 588/8 corank sequence | **Do not cite as re-verified** — only the leading term (4) confirmed live; a BertiniReal Python deflation bug was found, and patched reconstruction gave 2347 nonzero minors, not 588 | `DESIGN_NOTES.md:211-227` |
| Historical-curves timing "~3× documented baseline" | **Retracted** — invalid comparison (whole file vs single testset) | `DESIGN_NOTES.md:792-799` |
| Paper §2.3 gives the plane-curve critical system as `{f, ∂ₓf, ∂ᵧf}` | **The paper is wrong** (3 equations, 2 unknowns). The code is correct: `{f, f_y}` at `Topology.jl:369` | `notes/TALK_PREP.md:718-737` |

---

## 5. Flagged inconsistencies

Per the repository's own ground rules, these are reported rather than resolved.

1. **`paper_artifacts/data/results.json` reports test counts an order of magnitude below every other source.** It records fast `162/162` in 1m40.4s and full `178/178` in 3m15.1s (`results.json:260-305`), against `CLAUDE.md:75-77`'s fast `477/477` (~3.5 min) and full `537/537` (~30–34 min). The `results.json` figures date from a July 14 artifact run and appear simply stale, but nothing in the file marks them as such. Do not quote them.

2. **Full-suite pass count is reported as 536, 537, and 538 in different places.** `README.md:71` and `CLAUDE.md:75` say 537 ±1; `DESIGN_NOTES.md:428-429` records a CI run at 538/538; `paper_artifacts/figures/src/talk_script.md:7,91` says 538/538. The ±1 mechanism is genuinely documented (§4.5), but 536-vs-538 is a 3-wide spread, not ±1. Worth pinning with a live run before it goes on a slide.

3. **Full-suite runtime is reported as both ~30–34 min and 7m51s.** `README.md:71` / `CLAUDE.md:77` say 30–34 min; `DESIGN_NOTES.md:402-404` records 537/537 in **7m51s** locally after the `compile=:none` fix. If the 7m51s figure is post-fix and the 30–34 min figure is pre-fix, the README and CLAUDE.md are stale.

4. **`.cursorrules` states "Julia (v1.10+)"** while `Project.toml:27` and `README.md:27` both require **1.12**.

5. **`docs/BERTINIREAL_AUDIT.md:306` says isosingular deflation is "NOT IMPLEMENTED — no formal deflation."** That is contradicted by all of `src/Solver.jl:119-663` and by `CHANGELOG.md:24-30`, which adds it in v0.2.0. The audit is dated 2026-07-15 and `docs/ORCHESTRATOR_BRIEFING.md:116-128` marks it a historical snapshot, but the file itself carries no in-place staleness warning at that line.

6. **The torus is presented alongside tested fixtures in the paper artifacts but has no test coverage** (§3.3). Its residuals and counts are measurements, not enforced invariants, and its naked-edge count is explicitly non-deterministic.

---

## 6. Suggested talk framing

If the goal is a defensible academic presentation, the strongest honest structure is:

**What is genuinely solid.** The six-step pipeline is fully implemented for plane curves and ℝ³ hypersurfaces, registered, and covered by a topology-encoding test suite that asserts vertex-type counts and coordinates rather than "it ran." Pointwise residuals against the defining polynomial are excellent on smooth fixtures: sphere mean 2.85e-8, ellipsoid asserted `all(<=(1e-4))` at every welded mesh vertex, Taubin heart mean 1.03e-7. The isosingular deflation reproduces Hauenstein–Wampler's published Whitney-umbrella sequences `[3,2,0]` and `[3,1,1]` exactly, on the bare unsliced equation. Tolerances are calibrated against measured data with documented margins, not guessed.

**What is genuinely engineered.** Three pieces are non-obvious and defensible design work: the MidSlice-First strategy (never start tracking from a singularity), the deliberate refusal to trust HC's `is_success` at branch points, and the two-gate adaptive re-anchoring that makes non-radially-symmetric surfaces tractable.

**What is honestly open.** Multiplicity ≥ 2 detection is unreliable and fails *silently* — the cone returns an empty decomposition with no exception. There is no formal certification, and the standard tool (`certify`) would not certify the things that matter here anyway. Meshes are not watertight: 188 → 58 → 31–35 naked edges on the Taubin heart, with a third of the remainder undiagnosed. Deflation is diagnostic-only and does not yet feed singular-curve decomposition. Ambient dimension is capped at hypersurfaces in ℝ³.

The sharpest one-line framing available in the repository's own notes: **HGR trades a proof for a decomposition you can actually compute and look at** — its evidence of correctness is empirical, not deductive (`notes/TALK_PREP.md:141-152`).
