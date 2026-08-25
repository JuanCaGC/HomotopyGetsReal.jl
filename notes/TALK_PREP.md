# Talk prep — "When Homotopy Gets Real" (Albatross 2026)

**Purpose:** personal exam-style prep, not slide content. This is the material to
read *before* the talk so that questions about the mathematics and the code can be
answered from understanding rather than from the slide deck.

**Audience model:** heavily symbolic/exact-methods researchers — Gröbner bases,
`AlgebraicSolving.jl`, `msolve`, OSCAR-adjacent. Assume they know commutative
algebra and exact computation *better than you do*, and know numerical algebraic
geometry *less* well than you do. The failure mode to prepare for is not "they
don't understand homotopy continuation"; it is "they understand certification and
exactness precisely, and will ask what your numbers actually guarantee."

**Grounding rule used throughout:** every claim below is cited to a live file,
line, or measured number in this repo. Where a source document and the code
disagree, that is flagged rather than smoothed over (see §4).

---

## Scope caveat — read this once, then forget it

This document was written in the **Developer** environment (`CLAUDE.md`, "Roles"),
not in the external Correctness Auditor workspace. That matters for §4 in one
specific way, stated plainly rather than papered over:

- **Available and used here:** the full Bertini_real C++ source tree
  (`codigo_cplusplus/`, gitignored), the Bertini_real manual and Bertini user
  manual (`referencias_academicas/*.pdf`, gitignored), a Bertini1 binary at
  `/usr/local/bin/bertini`, the complete HGR source and test suite, and the
  current paper draft (`paper/HomotopyGetsReal_paper_current.md`).
- **Not available here:** a working `bertini_real` binary (`which bertini_real`
  → not found). The coordinate-level astroid cross-check described in
  `CLAUDE.md` mechanism 3 lives in the external Auditor workspace; this repo
  holds only the rendered figure (`paper_artifacts/astroid_bertini_real_sampled.pdf`),
  not the raw comparison data.

So §4 is an independent audit **against Bertini_real's source and published
manuals**, plus an independent read of HGR's own mathematics against its own
docstrings — not a re-run of Bertini_real. Everything in §4 was verified directly
in this session. Nothing in §4 is a restatement of `docs/DESIGN_NOTES.md`.

---

# Section 1 — Pipeline stage → math → code → validated example

The one table to have memorized. `Xy` in "Math" means ∂f/∂y throughout.

| Stage | Mathematical concept | Implementation | Best fixture |
|---|---|---|---|
| **Critical points** | Critical points of the projection π restricted to V(f). For a plane curve projected to x: the tangent is vertical, i.e. `{f, f_y}` — 2 equations, 2 unknowns. For a surface projected to z: ∇f ∥ (0,0,1), i.e. `{f, f_x, f_y}` — 3 and 3. Generically isolated, so a 0-dimensional solve is well posed. | `compute_critical_points`, `src/Solver.jl:732`. Curve augmentation is done by the *caller*: `src/Topology.jl:369`. Surface augmentation is automatic: `src/Solver.jl:749-752`. | **Sphere** — clean, 2 critical z at exactly ±1, asserted `atol=1e-6` (`test/test_surfacedecomposition.jl:64-65`). **Torus** is the counterexample: critical locus is 1-dimensional, so this stage fails outright. |
| **Boundary intersection** | Compactify an unbounded real variety by intersecting with a bounding region, so every edge has two endpoints. HGR fixes each variable at each face of an **axis-aligned box** and solves the reduced system. Bertini_real intersects a **sphere** instead. | `intersect_bounding_object`, `src/Solver.jl:824`; faces looped at `:843-917`. Explicitly **rejects** raw surfaces (`:835-841`) — a surface meets a box *edge* in a non-well-posed way. | **Nodal cubic** `y²−x³−x²` — 2 `Boundary` vertices at y=±4, asserted at `test/test_topology.jl:41-42`. |
| **Singularity classification** | Two independent numerical gates on the Jacobian's singular values σ₁≥…≥σₖ: rank deficiency (`rank < expected_rank`, rank = #{σᵢ > `jacobian_rank_tol`}) **or** a small trailing value (`σₖ < singular_value_threshold`). Distinguishes a genuine singular point of V(f) from a mere turning point of π. | `jacobian_rank_info` `src/Solver.jl:47`; the rule itself `_classify_vertex_type` `src/Solver.jl:71-80`. Genuinely `T`-generic — `BigFloat` SVD via GenericLinearAlgebra.jl. | **Nodal cubic** — exactly 1 `Critical` at (−1,0) and 1 `Singular` at (0,0), coordinates asserted (`test/test_topology.jl:33-40`). |
| **1D curve decomposition** | The load-bearing theorem: between consecutive critical values c₁<c₂ the real fiber π⁻¹(c) is topologically constant, so a *single* witness fiber at the midpoint determines the whole open interval. | `decompose_1d_curve` `src/Topology.jl:359`; `compute_midslice` `:61`; `connect_the_dots!` `:177`; `sample_edge` `:290`. | **Astroid** `(x²+y²−1)³+27x²y²` — 4 cusps, 4 edges; the paper's worked Bertini_real comparison (§6). |
| **Path tracking** | Parameter homotopy in a *real geometric* parameter (the slice coordinate), not an artificial t. Predictor–corrector with bisection on failure. Deliberately **not** gated on HC.jl's own success flag — paths legitimately die near branch points. | `build_tracker` `src/PathTracking.jl:84`; bisection core `_track_path_segment!` `:139`; `track_bidirectional` `:228`. | **Nodal cubic** (`test/test_pathtracking.jl`) — lands on `Critical`, `Singular`, and `Boundary` targets. |
| **3D surface sweeping** | Slice along z at critical values plus box bounds → slabs. In each slab decompose one 2D curve, then sweep each of its edges in z, keeping the track on the surface with a local linear **patch** `a(x−x₀)+b(y−y₀)=0`, `(a,b)=(f_y,−f_x)`, re-anchored when the gradient direction rotates too far. | `decompose_3d_surface` `src/SurfaceDecomposition.jl:1798`; `compute_critical_z_slices` `:48`; `slice_at_z` `:92`; `_robust_slice_at_z` `:286`; `build_patch_system` `src/FaceTracking.jl:95`; `track_face` `:574`. | **Ellipsoid** `x²+4y²+9z²−1` — asymmetric so it can't hide a swapped variable role; critical z at ±1/3; **every** mesh vertex residual ≤1e-4 asserted (`test/test_surfacedecomposition.jl:263`). |
| **Mesh welding** | Faces are tracked independently, so they carry no shared orientation. Pool and cluster vertices under `vertex_match_tol`, remap triangles, then fix winding so the geometric normal (p₂−p₁)×(p₃−p₁) agrees with ∇f. | `weld_mesh` `src/SurfaceDecomposition.jl:1581`; winding at `:1658-1664`. | **Sphere** — every welded normal outward, asserted `dot(n, centroid) ≥ 0` (`test/test_surfacedecomposition.jl:196-203`). **Taubin heart** is the honest one: 188 naked edges bare, 32–37 after stitching (measured across 28 trials, median/mode 34). |
| **Isosingular deflation** | At a singular point, append minors of the Jacobian to force the point to become regular in a larger system. `minor_size = expected_rank − corank + 1`. Iterate; the **corank sequence** either reaches 0 (isolated under deflation, isosingular dimension 0) or plateaus at k>0 (the point lies on a k-dimensional isosingular set). | `estimate_corank` `src/Solver.jl:119`; `deflate_once` `:209`; `_corank_plateau_hint` `:162`; `verify_isosingular_dimension` `:419`; `resolve_isosingular_dimension` `:597`. | **Whitney umbrella** `x²−y²z` — reproduces Hauenstein–Wampler exactly: tip (0,0,0) → `[3,2,0]`; handle (0,0,1) → `[3,1,1]`. Raw asserts at `test/test_isosingular_deflation.jl:82` and `:92`. |
| **Generic projection** | A random rotation Q ∈ SO(3) makes the *coordinate-alignment* degeneracies measure-zero. Implemented as pure change of coordinates: rotate the system, run the unchanged pipeline in the chart, map back. SO(3) not O(3), because a reflection would flip the welding convention. | `random_orthogonal_matrix` `src/Projection.jl:41`; `_rotate_system` `:115`; `_verify_projection_ok` `:172`; `_map_to_world` `:212`. Called from `src/SurfaceDecomposition.jl:1818-1830`. | **Torus** `(x²+y²+z²+3)²−16(x²+y²)` — the fixture that *only* works this way. Fixed-axis gives residuals ~1e⁸; `projection=:random` gives max ~2e-5. |

### Fixture cheat-sheet — what each one is *for*

| Fixture | Why it exists | Number to remember |
|---|---|---|
| **Sphere** | Smoothness control; tightest residuals | mean 2.85e-8 over 19,504 pts |
| **Ellipsoid** | Asymmetric — catches swapped variable roles | max 1.12e-4 at production density (**above** the 1e-4 test gate; see §3 Q7) |
| **Taubin heart** | Singular structure + the degenerate-slice mechanism | crit z ≈ {−1.0, 1.0, 1.0648, 1.2367}; naive z_mid=0 rejected, z=0.06 accepted after 5 retries |
| **Torus** | Positive-dimensional critical locus | needs `projection=:random`; mean 2.2e-6, p99 1.0e-5 |
| **Whitney umbrella** | Deflation ground truth vs published literature | `[3,2,0]` tip, `[3,1,1]` handle |
| **Astroid** | Coordinate-verified Bertini_real comparison | HGR 4 edges vs BR 6; reconciled, not a gap |

### Pipeline steps in full — curves, exact I/O types

Driver: `decompose_1d_curve(F::System, cfg::HomotopyConfig{T}; deflate::Bool=false) where {T<:AbstractFloat}` → `(Vector{NativeVertex{T}}, Vector{Edge{T}})`, `src/Topology.jl:359-401`. Its own inline comment states the real execution order — and it does **not** match the 1-2-3-4-5-6 numbering below: both merge operations (step 5) run *before* midslice/tracking (steps 3-4), front-loading dedup rather than cleaning up after. Numbering below is the conventional "six-step framework" naming (`src/Topology.jl:1-19`), not literal call order.

| Step | Function(s) — file:line | Input → Output | Mechanism |
|---|---|---|---|
| **1. Critical points** | `compute_critical_points`, `src/Solver.jl:733-791` | `(F::System, cfg::HomotopyConfig{T}; deflate::Bool=false, F_original=nothing)` → `Vector{NativeVertex{T}}` | Auto-augments a raw 3-var surface with `∂f/∂x, ∂f/∂y`; caller pre-augments a curve. HC `solve`, keeps near-real roots (`critical_point_tol`), Newton-polishes to `T`, classifies `Critical`/`Singular`, dedups via `cluster_vertices`. Finer-grained version of the table above's "Critical points" row — same function, now with the exact type signature. |
| **2. Boundary intersection** | `intersect_bounding_object`, `src/Solver.jl:825-922` | `(F::System, cfg::HomotopyConfig{T}; deflate=false, F_original=nothing)` → `Vector{NativeVertex{T}}` | Fixes each variable at each `bbox_x`/`bbox_y`/`bbox_z` face; Float64-casts before `solve` (HC's polyhedral start system can't take `Complex{BigFloat}` coefficients — `Solver.jl:856-872`); filters/polishes/classifies `Boundary`/`Singular`; dedups. |
| **3. Midslice** | `compute_midslice`, `src/Topology.jl:61-89` | `(F::System, x_left::T, x_right::T, cfg::HomotopyConfig{T})` → `Vector{Complex{T}}` | `F` must be a raw plane curve. Solves the 1-var fiber at `x_mid = (x_left+x_right)/2`, keeps near-real roots. This is the witness point everything else tracks from — never a critical/singular vertex directly. |
| **4. Connect the dots** | `connect_the_dots!`, `src/Topology.jl:177-232`; engine in `src/PathTracking.jl` — `build_tracker` `:84-92`, `track_bidirectional` `:228-259`, bisection core `_track_path_segment!` `:139-178` | `(F, x_left, x_mid, x_right, y_mid, edge_id, vertices [mutated], cfg)` → `Edge{T}` | Builds a `ParameterHomotopy`/`Tracker` with `x` as the parameter, tracks bidirectionally from the witness point. Poor accuracy or near-singular (rank/smallest-SV test) triggers adaptive bisection of the `x`-interval, not a discard — a bad landing only surfaces later as a vertex-match-tolerance miss. Landings resolve against the existing vertex list or append a new `Artificial` fallback vertex (`vertices` is mutated in place). |
| **5. Merge** | **Not one function.** `cluster_vertices`, `src/Clustering.jl:40-76`; `cluster_scalars`, `src/Clustering.jl:196-220` | `cluster_vertices(vertices::Vector{NativeVertex{T}}, tol::T)` → `Vector{NativeVertex{T}}`. `cluster_scalars(xs::Vector{T}, tol::T)` → `Vector{T}` | `cluster_vertices`: union-find over Euclidean distance ≤ `tol`, `Singular`-wins type merge, integer fields merged by minimum, floats by mean. `cluster_scalars`: sorts and merges consecutive x-values within `tol` before the per-interval midslice/tracking loop runs. The file's own "six-step framework" header (`Topology.jl:8-10`) names this step "Merge / GetMergeCandidates" — the two-function split is the code's actual shape, not an omission. |
| **6. Sample** | `sample_edge`, `src/Topology.jl:290-341` | `(F::System, edge::Edge{T}, cfg::HomotopyConfig{T})` → `Edge{T}` | Arc-length-equidistant resampling of `edge.sampled_points`, then re-projects every interpolated point back onto `f=0` via a capped-50-iteration Gauss-Newton correction (`_project_to_curve`, `:248-263`) — a 2026-07 fix, since raw chord interpolation was measured to land off-curve by residuals up to 0.28 on the Taubin heart. |

Supporting types: `NativeVertex{T}` — `id::Int, coordinates::Vector{Complex{T}}, v_type::VertexType, metadata::Dict{Symbol,Any}` (`src/Types.jl:48-55`). `VertexType` — `@enum Critical | Boundary | Singular | Artificial` (`src/Types.jl:26-31`). `Edge{T}` — `id::Int, left_vertex_id::Int, right_vertex_id::Int, sampled_points::Vector{Vector{T}}, is_singular::Bool` (`src/Types.jl:121-129`).

### Pipeline steps in full — surfaces, exact I/O types

Driver: `decompose_3d_surface(F::System, cfg::HomotopyConfig{T}; projection=nothing, rng=Random.default_rng(), incidence::Bool=false, deflate::Bool=false) where {T<:AbstractFloat}`, `src/SurfaceDecomposition.jl:1798-2042` → `(vertices, edges, faces, mesh)`, or with a 5th `incidence::SurfaceIncidence{T}` element when `incidence=true`.

| Step | Function(s) — file:line | Input → Output | Mechanism |
|---|---|---|---|
| **1. Critical z-values / slab boundaries** | `compute_critical_z_slices`, `:48-60`; `_slab_bounds`, `:411-423` | `(F::System [1 eq, 3 vars], cfg)` → sorted `Vector{T}` (both) | Reuses `compute_critical_points`'s 3-var auto-augmentation (`Solver.jl:748-753`) — no separate critical-point machinery in this file. `_slab_bounds` merges `bbox_z` endpoints with in-range critical z's, collapsing any closer than `min_slab_width` via `cluster_scalars` so no sliver slabs form. |
| **2. Robust midslice + retry** | `_robust_slice_at_z`, `:286-392`; underlying `slice_at_z`, `:92-126` | `_robust_slice_at_z(F, patch::NamedTuple, z_bottom::T, z_top::T, cfg; deflate=false)` → `(vertices_3d, edges_3d, z_mid::T)`. `slice_at_z(F, z_val::T, cfg; deflate=false)` → `(vertices_3d, edges_3d)` | `slice_at_z` substitutes `z`, runs `Topology.decompose_1d_curve` on the resulting 2-var curve, lifts back to 3D. Retry fires on two independent gates: a **topology gate** (an `Artificial`/`:endpoint_fallback` vertex co-occurs with a `Singular` one, but the slab's own quarter-point reference slices are *not* themselves singular) or a **gradient gate** (min anchor gradient falls below `z_mid_gradient_ratio_tol` × reference max). Retry perturbs the z-value only — never resampling density — stepping `+1,-1,+2,-2,...` frac-multiples of slab width off the midpoint, capped within 45% of the bounds; throws `ErrorException` if every attempt stays suspect rather than silently returning a bad slice. |
| **3. Critical-slice decomposition** | `_decompose_crit_slice`, `:559-577`; `_surface_critical_vertices`, `:597-598` | `(F, z_c::T, j::Int, cfg, vertex_registry::VertexRegistry{T}, e_offset::Int; deflate=false)` → `CritSlice{T}` (`boundary_index, z, vertices, edges, is_degenerate`) | Calls `slice_at_z` directly at the exact critical z — no perturbation, since the value is pinned by definition, not chosen. `is_degenerate = isempty(edges)` — flags fold/point-type boundaries (sphere/ellipsoid extremes: empty; cusps/tips: isolated `Singular` vertex, no edges). Only computed for interior slab boundaries, and only when `incidence=true`. |
| **4. Face sweep** | `build_patch_system` `FaceTracking.jl:95-116`; `patch_direction` `:245-250`; `build_face_tracker` `:263-272`; `sweep_face_bidirectional` `:524-552`; `track_face` `:574-645` | `track_face(F, patch, edge::Edge{T}, z_mid::T, z_bottom::T, z_top::T, face_id::Int, cfg)` → `Face{T}` (`id, mid_slice_z, boundary_edges, mesh_vertices::Matrix{T}, mesh_topology::Matrix{Int}`) | `patch_direction` gives a gradient-based transversal line `(a,b)=(f_y,-f_x)` — not radial. For each curve-arclength column of `edge`, sweeps bidirectionally in z and assembles a row-major `(n_z × n_curve)` grid, triangulating per quad. Winding is explicitly **not** normalized here — deferred to `weld_mesh` (docstring note, `FaceTracking.jl:569`). |
| **5. Mesh welding** | `weld_mesh`, `:1581-1676` | `(faces::Vector{Face{T}}, patch::NamedTuple, cfg; incidence=nothing)` → `GeometryBasics.Mesh` | Flattens all faces' `mesh_vertices`, clusters coincident points across faces at full `T`-precision via `cluster_points_indexed` (`Clustering.jl:232`), remaps triangles (dropping pinched ones), then per-triangle winding fix: flips to `(g1,g3,g2)` whenever `dot(normal, ∇f) < 0`, forcing every triangle's outward normal to agree with `+∇f` — since `track_face` emits no orientation guarantee of its own. |
| **6. Incidence-based stitching** | keyword confirmed: `incidence::Bool = false` on `decompose_3d_surface` (`:1798-1805`), forwarded to `weld_mesh`'s own `incidence::Union{Nothing,SurfaceIncidence{T}}=nothing` — **off by default in both places** | helpers: `_snap_boundary_points!` `:1032`, `_identify_edge_runs` `:1210`, `_reprojected_edge_targets` `:1000`, `_append_loft_triangles!` `:1373`, `_chained_edge_polylines` `:1464`, `_split_t_junctions` `:1141` | Snaps confident column landings onto shared crit-slice targets before clustering, builds ribbon/loft triangulation bridging face-boundary runs to crit-slice edges, then repairs T-junctions the loft introduces. Result per `_naked_mesh_edges`'s own docstring (`:1679-1707`): closes fold/point-type boundaries **completely**, reduces but does not fully close multi-face edge-type boundaries — exactly the Taubin heart's 188→32–37 naked-edge figure already in the table above (row "Mesh welding"), not a new number, just its mechanism spelled out. |
| **7. Visualization** | `plot_surface_decomposition(mesh::GeometryBasics.Mesh; ...)` `Visuals.jl:249-318`; `plot_surface_decomposition(faces::Vector{Face{T}}; ...)` `:338-384`; `interactive_3d_viewer` `:394-416` | mesh-method: welded `GeometryBasics.Mesh` → `Makie.Figure`. faces-method: pre-weld `Vector{Face{T}}` → `Makie.Figure` | Mesh-method: smooth per-vertex `:x`/`:y`/`:z` shading, no-ops gracefully on an empty mesh. Faces-method: flat per-face categorical coloring, explicitly does **not** apply `weld_mesh`'s winding correction (one-time `@warn`) — this is the exact function/code-path whose `n_faces=1` categorical-bin defect was root-caused and fixed earlier this session; not new information, just documented here for completeness. |

### How this codebase actually calls HomotopyContinuation.jl

One `using HomotopyContinuation` (`src/HomotopyGetsReal.jl:41`, bare, brings every export into scope). One deliberate fully-qualified exception: `HomotopyContinuation.evaluate(...)` at `Solver.jl:456` (elsewhere `evaluate` is called unqualified — stylistic disambiguation at that one site, not a structural difference). No `import` form anywhere.

| HC.jl function | Call sites | kwargs | Purpose | Role |
|---|---|---|---|---|
| `solve` | `Topology.jl:79` (midslice fiber); `Solver.jl:763` (`compute_critical_points`); `Solver.jl:874` (`intersect_bounding_object`) | `show_progress = false` at all three | Zero-dimensional solves feeding vertex classification | **(a) vertex solve** |
| `solve` | `Solver.jl:448`, inside `verify_isosingular_dimension` | `start_parameters=`, `target_parameters=`, `show_progress=false`, `compile=:none` | Walks one known point `x0_c` along a random-hyperplane-offset parameter homotopy to verify a single flagged `Singular` vertex's isosingular dimension | **hybrid — see Appendix flag** |
| `track` (via `Tracker`/`ParameterHomotopy`) | `PathTracking.jl:156`, built at `:88-90` | `compile` defaults `:all`, caller-overridable `:none`; `min_step_size = path_tracker_precision` | Per-step edge/face tracking, driven by `start_parameters!`/`target_parameters!` | **(b) edge/face tracking** |
| `certify` | **zero call sites**, confirmed by direct grep including comments | — | — | — |
| `@var` | **zero call sites** | — | Library never mints its own top-level variables; `System.variables` is always threaded through from caller-supplied systems | — |
| `Variable.(...)` | `Solver.jl:438` only | broadcast form | Names the `d` slack/hyperplane-offset parameters for the isosingular-verification homotopy above — the only site in `src/` that mints new symbolic variables | — |
| `solutions` | `Topology.jl:80`, `Solver.jl:764`, `Solver.jl:875` | `only_nonsingular = false` at all three, deliberately | Keeps singular solutions in-band; `Singular` vs. `Critical`/`Boundary` is decided downstream by this codebase's own rank test, not HC's flag | — |
| `path_results` | `Solver.jl:453` | — | `only(path_results(result))` — asserts the isosingular-verification homotopy produces exactly one result | — |
| `real_solutions`, `nonsingular` | **zero call sites** | — | — | — |
| `jacobian` | `Solver.jl:52, 230, 669` (symbolic); `PathTracking.jl:65` (numeric) | — | Feeds rank/SV-based vertex classification and the tracking bisection trigger | — |
| `differentiate` | 7 sites: `Projection.jl:175`; `FaceTracking.jl:107-109`; `Topology.jl:297,298,369`; `Solver.jl:753` | — | Builds Jacobian-augmentation equations for critical-point systems | — |
| `subs` | 7 sites: `Projection.jl:118`; `Topology.jl:77,137`; `Solver.jl:866,871`; `SurfaceDecomposition.jl:100` | — | The repeated "Float64-cast before `solve`" pattern — HC's polyhedral start system has no `Complex{BigFloat}` method, so any higher-precision `T` must be down-cast, then Newton-polished back up afterward | — |

**Discrepancy flagged against the paper** (`section3_background.tex`): the paper's clean two-role framing — zero-dim solve for vertices vs. parameter homotopy for positive-dimensional edges/faces — holds for 3 of the 4 `solve`/`track` call sites, and its augmented-system/rank-test description matches the code exactly down to the `{f, ∂_y f}` example. But `Solver.jl:448` doesn't fit either stated role cleanly: it's a `solve()` call carrying `start_parameters`/`target_parameters` (mechanically role (b)'s machinery) whose actual purpose is verifying a single already-flagged **vertex** (role (a)'s domain) — 0-dimensional in outcome (`only(path_results(...))`), never tracking a positive-dimensional object. The paper's §3 framing gives no reader a way to anticipate this call site exists. Not a contradiction of anything the code does — a real gap in the paper's description, flagged rather than smoothed over. See Appendix for the matching entry; not fixed here, same boundary as the naked-edges flag below.

---

# Section 2 — Mathematical background

Written for fielding questions. Each subsection ends with **where the subtlety
lives** — the thing a sharp question can expose.

## 2.1 Numerical algebraic geometry basics

**What a homotopy solve actually computes.** To solve a target system F, build a
start system G whose solutions are known, and connect them:
H(x,t) = (1−t)·γ·G(x) + t·F(x), with γ a generic complex constant (the "gamma
trick"). Each known solution of G is numerically tracked as t: 0→1 by a
predictor–corrector scheme. The γ is not decoration: it rotates the path off the
real/discriminant locus so that, for t<1, no two paths cross and none diverges —
this holds for generic γ, which is why "generic" is doing real work and not just
hedging.

What you get at t=1 is a **finite list of numerical approximations to the isolated
solutions**, each with an accuracy estimate. Three things you do *not* get for
free: (i) a proof that any given point is near a true solution, (ii) a proof that
the list is complete, (iii) anything at all about positive-dimensional components,
which the paths simply fail on.

**Witness sets and witness points.** For a positive-dimensional component of
dimension d, the standard object is a **witness set**: the component, a generic
linear slice of codimension d, and the finite set of intersection points. The
witness points are a finite handle on an infinite object; degree of the component
= number of witness points. Numerical irreducible decomposition (NID) computes
witness sets for every component of every dimension, using monodromy loops and the
trace test to group witness points by component.

**HGR does not implement any of this** (paper §2.4). Its input is one polynomial,
assumed to cut out an equidimensional real curve or surface. That is why every
validation fixture is an irreducible hypersurface. If asked "how do you handle a
reducible input" — the answer is "currently I don't, and the capability survey
shows exactly how it fails": on `three_concurrent_lines_reducible` the vertical
line x=0 is structurally invisible to x-parametrized slicing, and the triple point
is misclassified `Artificial` rather than `Singular`.

**Polyhedral vs. total-degree start systems.** Two ways to build G:

- *Total degree* (Bézout): for degrees d₁,…,dₙ, use G = {x₁^{d₁}−1, …, xₙ^{dₙ}−1},
  giving ∏dᵢ paths. Trivial to construct, but ∏dᵢ badly overcounts for sparse
  systems, so most paths diverge — wasted work.
- *Polyhedral* (Bernstein/BKK): the number of isolated solutions in the torus
  (ℂ*)ⁿ is bounded by the **mixed volume** of the Newton polytopes, which is
  ≤ ∏dᵢ and often far smaller. The start system comes from a regular mixed
  subdivision induced by a random lifting; you then track the "polyhedral
  homotopy" from it. Fewer paths, but a much more delicate construction.

HC.jl uses polyhedral by default, via `MixedSubdivisions.jl`. **This is the single
most load-bearing dependency risk in the project**, and it is where an
adversarial question will land hardest — see §3 Q3. Two concrete symptoms, both
live-confirmed:

1. The lifting sampler is **unseedable from HGR's API**. It draws from Julia's
   global `default_rng()`, which is seeded from OS entropy at process launch. So
   `decompose_3d_surface`'s own `rng` keyword reaches only the projection draw,
   not the solve. This is the mechanism behind the torus's run-to-run mesh-count
   variation (paper §4.3) — and the fix is known and one line
   (`Random.seed!` before the call), just not applied by default.
2. When the augmented critical system has an **identically-zero equation**, the
   polyhedral start-system construction throws
   `OverflowError: Cannot compute a start system`. `src/Projection.jl:12-17`
   documents exactly this crash class.

**What "numerical" costs you — be specific, this audience will press.** Do *not*
say "it's approximate but good enough." Say this:

| Guarantee | Gröbner / CAD / msolve | HGR as it stands |
|---|---|---|
| Solution completeness | Yes, by construction over the algebraic closure | No proof. Path count is a bound; missed paths are silent |
| Exactness of coordinates | Exact (rational, or real-algebraic via RUR/isolating intervals) | Float64, Newton-polished; optionally BigFloat |
| Certification that a point is near a true root | Not needed — exact | **Not performed** |
| Real-root counting | Yes (Sturm–Habicht, RUR, Descartes) | Imaginary part below a tolerance (`critical_point_tol = 1e-6`) |
| Topological correctness of the output complex | CAD: yes, by construction | Inferred from the critical-value theorem, checked *a posteriori* by residuals and test invariants |
| Cost | Doubly exponential worst case (CAD); Gröbner very sensitive to degree/variables | Roughly path-count × tracking cost; degrades gracefully |

The honest headline: **HGR trades a proof for a decomposition you can actually
compute and look at.** Its evidence of correctness is empirical (residuals,
literature cross-validation, topology-invariant tests), not deductive.

**Where the subtlety lives.** HC.jl *ships* a certifier — `certify`, based on
interval-arithmetic Krawczyk / Smale α-theory — which produces a rigorous proof
that a nonsingular solution of a square 0-dimensional system lies in a given box,
and can certify reality. `grep -rn "certify" src/` returns **nothing**. So the
right answer to "why no certification?" is not "certification isn't available in
this ecosystem" (false, and they may know it) but the far better answer in §3 Q2.

## 2.2 Critical-point-based decomposition

**The Morse-theoretic intuition.** Project the variety to one coordinate, π. As
the value c sweeps the line, look at the real fiber π⁻¹(c) ∩ V(f) ∩ ℝⁿ. Almost
always nothing interesting happens: the fiber's cardinality is locally constant and
the points move continuously. Topology changes only at **critical values** of
π restricted to V(f) — where the fiber is tangent to the projection, or where V(f)
itself is singular. This is Morse theory in flavour: critical points of a function
on a manifold are where sublevel-set topology changes.

So the algorithm writes itself. Compute all critical values c₁<…<c_m, add the
bounding values, and on each open interval (cᵢ, cᵢ₊₁) the fiber is topologically
constant. Take **one** witness fiber at the midpoint, and connect each of its
points outward to both ends by path tracking. This is the content of Lu–Bates–
Sommese–Wampler (curves) and Besana–Di Rocco–Hauenstein–Sommese–Wampler
(surfaces), and it is exactly why the implementation is "MidSlice-First"
(`src/Topology.jl:359`).

**Why MidSlice-First, and not the obvious thing.** You could track *from* the
critical vertices outward. You must not: a predictor–corrector tracker's step-size
control degrades sharply near a singular fiber, precisely because the Jacobian is
ill-conditioned there. Midslice points are smooth *by construction* — they sit
strictly between consecutive critical values. So singular points are only ever
**targets** of a track, never **sources**. This is the single cleanest design
answer in the whole project; have it ready.

**Why real is genuinely harder than complex — the actual pun.** This is worth
getting exactly right, because "when homotopy gets real" is a real mathematical
statement, not wordplay.

1. **ℝ is not algebraically closed.** All the machinery — Bézout, BKK, Bertini's
   theorem, genericity of γ — lives over ℂ. Tracking is *necessarily* complex.
   Real points are extracted only afterwards, by testing |Im(x)| against a
   tolerance. There is no homotopy that stays real.
2. **Real dimension is not determined by complex dimension.** A complex curve can
   have real points forming a 1-manifold, a finite set, or nothing at all. The
   fixtures show all three: `empty_curve` and `empty_surface` have empty real loci
   and decompose to nothing, gracefully.
3. **Genericity arguments weaken over ℝ.** Over ℂ, "generic" means "outside a
   proper Zariski-closed set", which is measure-zero *and* dense-complement. Over
   ℝ, a real random choice avoids a proper real algebraic subset almost surely,
   but the *real* structure — how many real points, whether a component is
   self-conjugate — is not a Zariski-generic property. Hence Bertini_real needs a
   self-conjugacy test (`checkSelfConjugate`) with no analogue in HGR at all.
4. **Real singularities are where the interesting geometry is, and where the
   numerics is worst.** Cusps, nodes, the Whitney umbrella's handle — exactly the
   features that make the picture worth drawing are exactly where the Jacobian
   drops rank and tracking degrades.

**What can go wrong.** Three named failure modes, all live-confirmed in this repo:

- **Degenerate critical values.** The critical fiber can itself be degenerate.
  Taubin heart at z=0: the polynomial factors as (x²+1.44y²−1)³, so the slice
  curve is a triple circle. The naive midslice of the slab [−1,1] is *exactly*
  z=0. Handled by `_robust_slice_at_z` (`src/SurfaceDecomposition.jl:286`), which
  perturbs and retries — accepting z=0.06 after 5 tries. Note this is a
  perturb-and-retry *workaround*, not a structural fix, and the gate threshold
  (`z_mid_gradient_ratio_tol = 0.01`) was calibrated on this one fixture:
  measured 0.013 for a good slab vs 0.0014–0.0057 for bad ones.
- **Multiplicity ≥ 2.** At a cusp, the critical point is a double root of the
  augmented system, and HC.jl's polyhedral solve finds such roots unreliably.
  Measured, not theorized: across 5 independent astroid runs, **2 came back with
  1 or 2 cusps missing entirely**, replaced by fallback vertices. On the node curve
  `y²−x²`, `compute_critical_points` on `[y²−x², 2y]` returns **zero solutions** —
  both paths report `excess_solution`.
- **Positive-dimensional critical locus.** The killer, because the whole method
  assumes critical points are isolated. Worth deriving live if asked: for the
  torus, f_x = 4x[(x²+y²+z²+3)−8] and f_y = 4y[(…)−8], so f_x=f_y=0 on
  x²+y²+z²=5; with f=0 this gives x²+y²=4, z=±1 — two whole **circles**. A
  0-dimensional solve on a positive-dimensional solution set is ill-posed, and the
  output is garbage (residuals ~1e⁸), not an error. This is precisely what generic
  projection fixes.

**Where the subtlety lives.** The critical-value theorem guarantees the fiber is
constant on the *open* interval. It says nothing about whether you found **all**
the critical values. If one is missed, two slabs silently merge and a topological
feature disappears with no error raised — the cone fixture returns 0 vertices, 0
faces, and no exception. Every correctness claim in the project is conditional on
the completeness of the critical-point solve, which is exactly the step with the
weakest guarantee. That conditional is the honest shape of the whole method, and
saying so out loud is much stronger than being caught not knowing it.

## 2.3 Isosingular deflation

**What corank means geometrically.** At a point x on V(f), the Jacobian J(x) is
the linearization. corank = (expected rank) − rank J(x) measures how far the
linearization is from telling you the local dimension. At a smooth point, rank is
maximal and the implicit function theorem applies: V(f) is locally a manifold of
the expected dimension. At a singular point the rank drops, and the first-order
data is no longer enough to pin down local structure — a node and a cusp both
report corank 1 under the row-rank convention (`src/Solver.jl:112-117` says this
explicitly).

**What deflation actually does.** If rank J(x) = r and the point is singular, then
*every* (r+1)×(r+1) minor of J vanishes at x — that is what rank deficiency
means. So those minors are **new equations satisfied at x**, and adjoining them
cuts down the solution set without losing x. Formally: append every
`minor_size = expected_rank − corank + 1` sized minor of the Jacobian
(`src/Solver.jl:222`), then recompute corank on the enlarged system. Iterate.

The **corank sequence** is the diagnostic. Hauenstein–Wampler Lemma 3.4: it is
nonincreasing — asserted as a hard `ArgumentError` at `src/Solver.jl:162-168`. Two
outcomes:

- Sequence reaches 0. The point has become **regular in the deflated system** —
  isolated under deflation, isosingular local dimension 0. Whitney tip:
  `[3,2,0]`.
- Sequence **plateaus** at k>0. The point lies on a genuinely k-dimensional
  isosingular set that deflation does not reduce. Whitney handle (0,0,1) on the
  singular z-axis: `[3,1,1,…]`, isosingular dimension 1. **This is the correct
  answer, not a failure to converge** — a point worth saying explicitly, because
  a plateau looks like non-termination.

**Why a plateau is not enough — the trap.** The obvious test ("two consecutive
coranks agree ⇒ done") is **necessary but not sufficient**, and Hauenstein–Wampler
say so themselves, with the counterexample family f_{k,l} = [x^k, y^l] exhibiting a
plateau that later drops. HGR names this honestly: the cheap test is called
`_corank_plateau_hint` (`src/Solver.jl:162`), documented as
"NECESSARY-BUT-NOT-SUFFICIENT", used only to decide whether it is worth paying for
a real verification. The authoritative test is
`verify_isosingular_dimension` (`:419`), which draws d generic hyperplanes through
the point, tracks, and checks the residual **against the original undeflated
system**.

Two cheaper alternatives were tried and rejected on measured evidence, which is a
good story to have: discarding equations from the accumulated system gives a false
positive on the Whitney handle, and replacing the system by a randomized reduction
R(f)=A·f gave a **4-of-6 false-positive rate** with residuals 0.08–1.79. Hence
verification insists on the full original system and a fresh hyperplane draw per
retry.

**How this relates to the actual Hauenstein–Wampler theory, cross-checked against
Bertini_real.** Three things I verified directly in the source and manual:

1. **The minor-size formula matches, modulo homogenization.** Bertini_real,
   `codigo_cplusplus/src/symbolics/isosingular.cpp:195`:
   `minorSize = numVars - declarations[1] - nullSpaceDim + 1`, where `nullSpaceDim`
   is the corank and `declarations[1]` counts homogeneous variable-group
   declarations. HGR's `expected_rank − corank + 1` is exactly the
   `declarations[1] = 0` case. HC.jl works affinely, so this is correct — and the
   docstring at `src/Solver.jl:193-194` already claims precisely this. **Claim
   verified.**
2. **The filter on minors is structural, not numerical — and this is the
   subtle part.** HGR keeps a minor iff `!iszero(expand(det(...)))`, i.e. iff it
   is not the zero *polynomial*. It must be structural: by definition of rank,
   every minor of that size vanishes *numerically* at x, so filtering on value at
   x would discard all of them. Bertini_real's Matlab path uses
   `simplify(det(...)) ~= 0` — the same structural test through a different CAS.
   **Claim verified.** This one is quietly a good answer for a symbolic audience:
   it is genuinely a symbolic computation sitting inside a numerical pipeline.
3. **Where the two genuinely differ — the witness point.** Bertini_real's own
   manual describes deflation as: "sub-determinants of the system are recursively
   added until the rank of the Jacobian stabilizes, where **this rank is computed
   using a witness point for the component being deflated**." HGR computes corank
   at *the specific singular point x₀*, with no witness point and no witness set —
   and the paper (§6) presents this as a strength, correctly noting it reproduces
   the published Whitney sequences with no witness-set stage. But see §4.3: the
   witness point in Bertini_real is not bureaucratic overhead, and this is the
   sharpest available question about the deflation subsystem.

**Where the subtlety lives.** Two places.

- **Combinatorial blowup — the thing this audience will smell immediately.** Round
  k appends C(rows, m)·C(cols, m) symbolic determinants. `test/test_isosingular_deflation.jl:86-90`
  records that round 5 of the Whitney handle would need **~5.3×10⁹ symbolic
  minors and never terminates** — it had to be killed, and the test is hard-capped
  at 2 rounds because of it. Deflation is the one part of HGR with genuinely
  symbolic, genuinely exponential cost. A Gröbner-basis audience will recognise
  this blowup instantly; own it rather than let them find it.
- **Scope.** Deflation is **diagnostic only** (paper §3.5, §6). It stamps
  isosingular dimension, verdict, and corank sequence into vertex metadata. It does
  not decompose a singular curve as a first-class object, and it does not move a
  single mesh vertex. Bertini_real's whole point in Brake et al. 2014 was to *use*
  the diagnosis for singular-curve decomposition. This is the largest real scope
  gap in the project, and the paper already says so.

## 2.4 Generic projection — the transversality argument

The claim to defend: *a random rotation avoids degenerate configurations almost
surely*. Here is the argument stated precisely enough to hold up.

**Setup.** Pick Q ∈ SO(3) and work in the chart x = Qx′, so the chart system is
f′(x′) = f(Qx′). The pipeline slices along the chart's z, which is the world
direction given by Q's third column. Nothing else in the pipeline reads world
coordinates — every frame-dependent decision goes through `cfg.bbox_*` or a
positional variable index (`src/Projection.jl:4-8`). So choosing a projection is
*exactly* choosing a rotation, and the pipeline itself needs no modification.

**The argument.** Consider the set of "bad" rotations B ⊂ SO(3) — those for which
the sweep direction fails to be generic with respect to V(f): the critical locus
of π_z is positive-dimensional, or a critical fiber is degenerate, or an augmenting
partial vanishes identically. Each of these is a **Zariski-closed condition on the
entries of Q** — the critical locus being positive-dimensional forces the
vanishing of appropriate resultants/minors in the coefficients of f′, which are
polynomial in Q's entries. For a fixed f that is not itself degenerate, at least
one Q is good, so B is a **proper** algebraic subset of SO(3). A proper algebraic
subset of a smooth irreducible real variety has measure zero for Haar measure.
Since `random_orthogonal_matrix` (`src/Projection.jl:41`) draws **Haar-uniformly**
on SO(3), P(Q ∈ B) = 0. Hence: a random rotation is generic almost surely.

Two implementation details that make the argument real rather than notional:

- **Haar-uniformity is not free.** A plain LAPACK `qr()` of a Gaussian matrix is
  *not* Haar-distributed — the sign convention in the R factor biases it
  (Mezzadri, *Notices AMS* 54 (2007)). Multiplying columns by `sign(diag(R))`
  restores Haar on O(n); negating one column when det<0 is a measure-preserving
  right translation onto SO(3). Verified empirically: 20,000 samples,
  E[(col)_z²] = 0.3323 against the theoretical 1/3. If someone asks "is your
  random rotation actually uniform," this is the answer.
- **SO(3), not O(3), and for a real reason.** For orthogonal Q,
  cross(Qa, Qb) = det(Q)·Q·cross(a,b), while gradients map as ∇f_world = Q·∇f_chart.
  So a reflection (det=−1) would silently flip every mapped-back normal against
  `weld_mesh`'s ∇f-alignment convention. `_resolve_projection` therefore
  **rejects** det<0 rather than auto-correcting it (`src/Projection.jl:94-98`).

**Answering "how do you know a random projection won't just move the problem
somewhere else?"** Give it in three parts.

1. **The degeneracies being avoided are properties of the *pair* (f, direction),
   not of f alone.** The torus's fold circles are a *coordinate-alignment*
   artifact: the torus has no real singularities at all (all 8 vertices classify
   `Critical`), and its critical locus is 1-dimensional *only* for sweep
   directions aligned with the symmetry axis. That is a measure-zero set of
   directions. So there is no "somewhere else" to move it to — genericity is not
   trading one bad direction for another, it is leaving a measure-zero set.
2. **Genuine singularities are *not* removed, and shouldn't be.** A cusp of the
   Taubin heart is a property of the variety. Rotating changes which coordinates
   express it, not whether it exists. Generic projection is *only* claimed to fix
   alignment artifacts. Be crisp about that boundary — it is where an overclaim
   would get caught.
3. **Then concede the honest part.** Empirically, `projection=:random` fixed the
   torus and `horn_torus` cleanly but **did not** fix the `cone`. And the
   geometric argument that looked like it explained the cone ("the singular point
   is at the origin, fixed under any rotation, so rotation can't help") **also
   predicts failure for horn_torus, which in fact succeeds.** That prediction is
   recorded as falsified and the divergence is *not yet understood*. Saying "my
   own explanation for one failure was falsified by another fixture and I haven't
   resolved it" is far stronger than offering a tidy story that a good questioner
   could break.

**Where the subtlety lives — and this is the sharpest single item in §2.** The
almost-sure argument above is a **mathematical** argument. The code does **not**
verify it. `_verify_projection_ok` (`src/Projection.jl:172`) checks exactly one
narrow thing: whether ∂f′/∂x or ∂f′/∂y is the **identically-zero polynomial**,
tested at 3 fixed complex probe points. It deliberately does not check ∂f′/∂z, and
explicitly leaves positive-dimensional critical loci "to the downstream machinery"
(`:166-169`). That is a crash guard, not a transversality test. See §4.2 — I
verified Bertini_real's corresponding check is a genuinely *different and stronger*
condition, and this is worth knowing before someone else notices.

---

# Section 3 — Anticipated questions, with real answers

Ordered by likelihood for *this* audience. Answers are written to be defensible,
not comfortable.

### Q1. "Why not just use a Gröbner basis / CAD / msolve for this?"

Because the output is a different object, and the exact tools don't produce it at
the scale where the pictures get interesting.

The goal is not "solve a system" — it is a **cell complex plus a triangulated mesh
of the real points**, suitable for visualisation and for downstream topological
questions. CAD does produce a genuine cell decomposition, and where it terminates
it is strictly better than this: it's exact and complete. But it is doubly
exponential in the number of variables, and for a degree-9 surface like the Taubin
heart in 3 variables a full CAD is already impractical. Gröbner bases and `msolve`
are superb at the zero-dimensional subproblem — and in fact the *critical-point*
step is exactly a zero-dimensional system where they'd be excellent — but they do
not by themselves give you the connectivity: which arcs join which critical points,
how the sheets glue across a critical slice. That's what path tracking supplies.

Then concede properly: for small enough problems, exact methods dominate, and this
audience's tools would give a *better* answer than mine. The regime where this
approach wins is moderate degree, 3 variables, where you want a mesh — and the
honest current demonstrated scale is trivariate low-degree surfaces (paper §6,
Table 5), against Bertini_real's degree-630 14-variable curve.

**The move that makes this a good answer rather than a defensive one:** propose the
hybrid. The critical-point solve is the least reliable stage (Q3) and is *exactly*
a zero-dimensional polynomial system. Computing it with `msolve`/`AlgebraicSolving.jl`
would give a **certified, complete** set of critical points, then hand off to
homotopy tracking for the connectivity. That fixes the project's weakest link with
this audience's strongest tool. And the plumbing already exists: `examples/oscar_integration.jl`
takes an OSCAR `MPolyIdeal` straight to an HC.jl `System`, verified end-to-end on
the ellipsoid (mean 7.4e-7, median 2.93e-8 against the baseline's 2.70e-7 /
2.79e-8). This is the most promising collaboration to offer from the stage.

**Note (no longer on a public slide, keep in mind for Q&A):** Rémi Prébet —
plenary speaker and tutorial instructor at this workshop — is a listed
maintainer of `AlgebraicSolving.jl`. He is not a RAGlib author (verified
against RAGlib's own source). Don't misspeak about which project he's
connected to if this comes up live.

### Q2. "How do you know your numerical results are actually correct, with no certification step?"

Concede the premise immediately: **there is no certification step, and "correct" is
the wrong word for what I have.** What I have is four independent kinds of
evidence, which is not the same as a proof.

1. **Pointwise residuals against the defining polynomial.** Every mesh vertex is
   plugged back into f. Sphere: mean 2.85e-8 over 19,504 points. Ellipsoid: every
   vertex asserted ≤1e-4 (`test/test_surfacedecomposition.jl:263`). Taubin: mean
   1.03e-7. This is a *necessary* condition — it proves points are near the
   variety. It says nothing about completeness.
2. **Cross-validation against published literature.** The deflation subsystem
   reproduces Hauenstein–Wampler's Whitney umbrella sequences exactly: `[3,2,0]`
   and `[3,1,1]` (`test/test_isosingular_deflation.jl:82`, `:92`). Independent of
   my own code's opinion of itself.
3. **Coordinate-level comparison against Bertini_real's actual output** on the
   astroid: BR's 6 real edges vs HGR's 4, reconciled exactly — the 2 extra come
   from BR's random projection subdividing 2 of the 4 cusp-to-cusp arcs at smooth
   projection-critical points S±, confirmed on-curve (|f(S±)| ~ 3.2e-10). Same
   1-manifold, different partition. And a *third* independent confirmation:
   Amethyst–Hauenstein–Wampler report "four singular points connected by four
   edges" for this curve.
4. **A test suite that encodes topology, not just "it ran."** 538 assertions,
   including exact vertex-type counts and coordinates.

Now the part that makes this answer land. **HC.jl ships a certifier** — `certify`,
based on interval-arithmetic Krawczyk / Smale α-theory — and I don't call it. I
should, for the critical-point solves: they're square and zero-dimensional, which
is exactly its domain, and it would upgrade stage 1 from "these are the points I
found" to "each of these is provably a true solution, and provably real." That is
a concrete, near-term improvement and I'd take the prompt.

But be honest about what it would *not* fix, because that's the interesting part:
`certify` handles **nonsingular isolated** solutions of square systems. It cannot
certify (i) the *completeness* of the critical-point list — the actual weak link;
(ii) the path-tracked edges and faces, which are positive-dimensional and where
most of the geometry lives; or (iii) the singular points, which are non-isolated
in the relevant system and are precisely the hard part. So certification would
strengthen a real link in the chain without closing it. Certifying a *decomposition*
is a substantially harder open problem than certifying a root.

### Q3. "What happens when your polyhedral solver fails to find all critical points?"

This is the question I'd least like to be asked and the one with the most
documented detail. The answer: **it fails silently, and I can quantify it.**

Four live-confirmed instances:

- **Node curve** `y²−x²`: `compute_critical_points` on `[y²−x², 2y]` returns
  **zero solutions**. Both paths report `excess_solution`; neither converges to the
  double root at the origin.
- **Astroid cusps** (multiplicity 2): 12 raw paths collapse to 4 physical points.
  Measured properly across 130 trials, cusps were lost in 34 (~26%, roughly
  1 in 4) — replaced by `Artificial` fallback vertices. Not coordinate jitter;
  missing features.
- **Cone**: `compute_critical_z_slices` returns `Float64[]` despite a
  mathematically genuine apex. Default call yields a **completely empty
  decomposition — 0 vertices, 0 faces, and no exception.** As the notes put it,
  "no crash this time, just nothing, which is arguably worse to notice."
- **Whitney umbrella as a surface**: same empty critical-z set; the surface fixture
  throws (1 of 24 in the capability survey).

The working hypothesis is that HC.jl's default polyhedral solve struggles when the
augmented system has a solution of **multiplicity ≥ 2** at the sought point —
structurally guaranteed at cusps, reducible crossings, and conical apexes. I'll
flag that this is a *hypothesis*: it has not been verified against HC.jl's
internals.

Three mitigations, honestly ranked:

1. **Generic projection**, when the cause is coordinate alignment. Fixed the torus
   and horn_torus cleanly. Did **not** fix the cone — after 8 perturbed retries it
   fails loudly, which is better than silence but is not a result.
2. **`_robust_slice_at_z`'s perturb-and-retry**, for degenerate slices. Works, and
   the gate threshold is calibrated on one fixture.
3. **The rejected one, worth mentioning because rejecting it was correct.** An
   asymmetric bbox makes the cone *appear* to work: median residual 1.2e-7, but
   **p99 = 12.2, max = 16.0** with a visible seam. The internal note is blunt:
   "more dangerous than the original empty-mesh failure, not a partial fix. Do not
   present this as a usable workaround under any framing." A workaround producing
   plausible wrong answers is worse than a visible failure — that's a principle
   worth stating out loud.

And the real fix, which loops back to Q1: replace stage 1 with a certified exact
solver. Multiplicity-2 roots are precisely what Gröbner methods handle *correctly
and completely* and what polyhedral homotopy handles worst. This is the single
place where this audience's tools would most improve the project.

### Q4. "Does this scale to higher-degree / higher-dimensional systems?"

No, and I'd rather give the shape of the ceiling than a number.

**Demonstrated scale** is trivariate surfaces of low degree (sphere, ellipsoid,
Taubin heart degree 9, torus degree 4), against Bertini_real's reported degree-630
curve in 14 variables. That gap is a difference in development stage, but it should
not be understated.

Four distinct scaling walls, which are *not* the same wall:

1. **Ambient dimension is architecturally capped right now.** Surfaces are
   hypersurfaces in ℝ³, hard-enforced: `compute_critical_z_slices` throws unless
   there are exactly 3 variables and 1 equation (`src/SurfaceDecomposition.jl:49-53`).
   There's no complete-intersection or higher-codimension path at all.
2. **The bounding box scales worse than a sphere — and Bertini_real's own authors
   say so.** HGR uses an AABB, requiring a separate boundary solve per face: 2N
   faces in N variables. The Bertini_real manual (§A.2.2) records that they
   *started* with a bounding box and abandoned it: "for higher dimensions, it meant
   doing more and more curve decompositions, with each of the 2N planes. It got
   messy." A sphere is a single degree-2 equation regardless of N. So this is a
   known-worse choice on scaling grounds, made by their own account.
3. **Deflation blows up combinatorially.** Round k appends C(rows,m)·C(cols,m)
   symbolic determinants. Whitney handle round 5: **~5.3×10⁹ minors, never
   terminates.** The test is hard-capped at 2 rounds for this reason. This is a
   genuinely symbolic bottleneck — same flavour as Gröbner blowup — and it's the
   dominant runtime cost: `test_isosingular_deflation.jl` alone is ~76% of the full
   suite (18m15s of ~24m).
4. **Path count grows with mixed volume**, the mildest wall, and the one where
   homotopy genuinely beats CAD asymptotically.

Also worth stating: Bertini_real's own open challenge is an eight-polynomial,
ten-variable Burmester surface, blocked by the size of the symbolic determinant in
the critical-curve computation. So the *architecture* has a ceiling too, not just
this implementation.

### Q5. "Your test suite has 538 assertions and you say the count varies. Doesn't that mean your tests are nondeterministic?"

Yes, in a specific and understood way, and it's worth separating two things.

The ±1 variance (537 vs 538) is real and traced: the astroid isosingular-deflation
loop fires a variable number of `@test` invocations depending on how many `Singular`
vertices the live solver produces that run (`test/test_isosingular_deflation.jl:126-137`).
So the *assertion count* varies because the *loop trip count* varies. No assertion
changes from pass to fail.

But that's a symptom of something more interesting: cross-process nondeterminism
is real throughout, with an identified cause. HC.jl's polyhedral start-system
construction draws its random lifting from Julia's global `default_rng()`, seeded
from OS entropy per process, with no seeding path exposed through `solve`. So:
torus mesh counts vary by tenths of a percent; naked-edge counts on the Taubin
heart range 32–37 across runs (measured across 28 trials, median/mode 34); and
on `folium_descartes`, **1 of 3 identical runs
orphaned an on-curve `Boundary` vertex referenced by zero edges, with no error**.

The fix is known and one line — `Random.seed!` before the call makes output
bit-identical, confirmed directly — and it's not applied by default. So the
reported figures reflect current default behaviour, not an inherent limit. I'd
rather report the variance than a single run's numbers as if they were invariants.

### Q6. "You say the Taubin heart's cusps are singular points, but your own default output doesn't classify them as `Singular`. Which is it?"

The cusps are genuinely singular, and the default return value does not report them
as such. That's a real usability defect, and the paper states it (§4.2, §6).

Concretely: the plain `decompose_3d_surface` tuple gives 14 vertices — 10
`Critical`, 4 `Artificial`, **zero `Singular`**. To see the cusps you must enable
the surface-level incidence diagnostic and union its reported critical vertices
with the plain set, which yields 20 vertices: 14 `Critical`, 4 `Artificial`, 2
`Singular` at the cusps, at identical mesh statistics.

Why: classification happens per z-slice, against the *2D slice curve*
`f(x,y,z_val)=0`, not against the 3-variable surface. A cusp of the surface need
not be a singular point of a generic 2D slice through it. Asking the
3-dimensional question requires `_surface_critical_vertices`, which is the separate
diagnostic. The design reason is documented (`src/SurfaceDecomposition.jl:78-82`)
and defensible — deflating against the surface would be a different,
higher-dimensional question — but the *default* being the less informative answer
is not defensible, and Bertini_real is better here: singular-curve decomposition is
integrated into its main pipeline, so its users don't have to know to ask twice.

### Q7. "Your paper reports a maximum residual of 1.12×10⁻⁴ on the ellipsoid, but your test gate is 10⁻⁴. Your own example fails your own test."

Correct, and the paper says so rather than rounding it away (§5.1).

The precise situation: the 1e-4 gate is applied to the **coarse** mesh used in
automated testing, where the max is 1.000e-5. The 1.12e-4 figure is at
**production** density (19,500 points), which the gate is not applied to. So no
test is failing — but the gate is also not a claim about arbitrary sampling
densities, and presenting it as one would be wrong.

Three honest observations. First, it is genuinely counterintuitive that *more*
sampling gives a worse maximum: it's consistent with accumulated floating-point
error at higher density in high-curvature regions, and note the *median* stays at
2.79e-8 — only the tail moves. Second, it's 1 or 2 outliers, not a distribution
shift: at coarse density, ~170× above p90. Third, it's four orders of magnitude
below the bounding box's geometric scale and the topology is unchanged (identical
vertex, edge, and face counts).

The fair criticism I'd accept: a pass/fail gate calibrated on the coarse mesh and
not re-validated at production density is a weak gate, and the ellipsoid was never
the fixture the adaptive re-anchoring machinery was validated on — that was the
Taubin heart. That's a documented hole in my own validation coverage, not
something the residual number itself resolves.

### Q8. "Your mesh isn't watertight. For a decomposition into cells, isn't closure the whole point?"

It isn't watertight, and I'd rather give the numbers than the adjective.

Instrumented naked-edge counts on the Taubin heart: **188** with plain
`weld_mesh`; **58** with incidence snapping; **32–37** with the loft stage
(measured across 28 trials, median/mode 34) — and that range is itself not
reproducible run to run. Enabling incidence stitching
closes the fold- and point-type singular features (the two tips) **completely**;
what survives is at boundaries where **more than two faces meet along a shared
edge** — the singular notch and a saddle pair.

Why it's hard: faces are tracked independently and welded afterward by clustering
under `vertex_match_tol`. Where exactly two faces meet, samples from both land
within tolerance and weld. Where three or more meet, the sample points don't
coincide pairwise and no single clustering radius fixes it. Three separate
investigations attributed roughly a third to genuine cross-edge junctions, a third
to seam artifacts, and **a third remains undiagnosed** — confirmed not to be a
resolution artifact, since the `sample_edge` Newton fix left the count unchanged.

Now the sharpest part of the answer: **for a cell decomposition, the combinatorial
incidence structure is the primary object, and the triangle mesh is a rendering of
it.** The vertex/edge/face complex is correct — that's what the topology assertions
check. Naked edges are a defect in the *mesh realization*, not evidence of a wrong
decomposition. That's a real distinction and not a dodge. But I won't overclaim it
either: a rendering with visible holes is a poor advertisement for a correct
complex, and full closure is currently deferred, not solved.

### Q9. "You depend on HomotopyContinuation.jl for the part that's hardest to get right. What happens if it's wrong?"

Then I am wrong, and I have limited ability to detect it. That's a genuine
architectural dependency, not a hedge.

What I *can* do, and do: check residuals against f independently of what the
solver claims (a tracker that lies about success still can't fake a small |f(p)|);
refuse to trust HC.jl's success flag for path termination, since paths legitimately
terminate near branch points without it (paper §3.7); and cross-validate against a
completely independent implementation — Bertini_real, built on Bertini1, sharing no
code — on the astroid, plus against published literature on the Whitney umbrella.

What I *can't* detect: **missing solutions.** If the polyhedral solve silently
returns 3 of 4 cusps, no residual check will notice, because the 3 it returned are
all fine. That's Q3, and it's the failure mode this architecture is worst at. It's
also precisely why an exact certified solver for stage 1 is the most valuable
improvement anyone in this room could suggest.

### Q10. "Your isosingular deflation is 'validated exactly' against Hauenstein–Wampler, but you also say it's diagnostic only. What does the validation actually buy you?"

Fair, and the two are less connected than "validated exactly" makes them sound.

The validation is real and narrow: on the bare equation `x²−y²z`, the corank
sequences match the published ones exactly — `[3,2,0]` at the tip, `[3,1,1]` at the
handle. That's genuine cross-validation against an independent published source,
and it caught real bugs (a prior `deflation_stabilized` wrongly returned `false`
for `[1,1,1]`; a docstring claiming "every case resolved in exactly 1 attempt" was
retracted after measuring 1,2,1,4,1,1,4,1,1,1 across 10 trials).

What it buys: confidence that the *diagnosis* is right. What it does not buy: any
change to the output. The isosingular dimension is stamped into vertex metadata and
nothing downstream reads it. Coordinates don't move; the mesh doesn't change.

So the accurate framing — and the one I should use on stage — is that the deflation
subsystem is a **correctly-implemented and independently-validated component that
is not yet wired into the pipeline's output.** Bertini_real's contribution in Brake
et al. 2014 was precisely to *use* this diagnosis to decompose singular curves as
first-class objects. I've reproduced the theory and not yet the consequence. And
note the `[3,1,1]` assertion is at a **hard 2-round cap**: it's consistent with the
published plateau, but my test does not itself demonstrate the plateau continues,
because round 5 doesn't terminate (§4.4).

---

# Section 4 — Independent audit: where my own read differs from the paper and codebase

Everything here was verified directly in this session against Bertini_real's C++
source, its published manuals, or HGR's own code. None of it is a restatement of
`docs/DESIGN_NOTES.md`. Ordered by how likely it is to matter on stage.

### 4.1 The paper's §2.3 states the plane-curve critical system incorrectly — fix this before the talk

Paper §2.3 (`HomotopyGetsReal_paper_current.md:50`) reads:

> "…this is the same augmentation used in §3.4 (e.g. **{f, ∂ₓf, ∂ᵧf}** for a plane
> curve projected onto x)."

That is wrong, and wrong in a way this audience will catch instantly by counting:
**{f, f_x, f_y} is 3 equations in 2 unknowns** — overdetermined, generically empty.

The correct system is **{f, f_y}**. Derivation: J stacks ∇f = (f_x, f_y) with
dπ = (1,0); det[[f_x, f_y],[1,0]] = −f_y, so the augmented system is {f, f_y}.
Geometrically: critical points of π_x on a plane curve are where the tangent is
vertical, i.e. f_y = 0.

**The code is correct** — `src/Topology.jl:369` is
`F_aug = System([f, differentiate(f, y_var)], F.variables)`, exactly {f, f_y}. So
this is a typo in the paper's exposition, not a bug. It looks as though the
surface case {f, f_x, f_y} (which *is* correct, 3 equations in 3 variables, and is
what `src/Solver.jl:749-752` builds) was pasted into the curve sentence.

Why this matters disproportionately: it is in §2.3, the *mathematical background*
section, in the one parenthetical a symbolic-methods reader will check to see
whether the author understands the geometry. Getting the equation count wrong there
costs more credibility than any of the genuine limitations elsewhere in the paper.
Worth correcting in the manuscript and being ready to state the right version from
memory.

### 4.2 HGR's `_verify_projection_ok` is *narrower* than Bertini_real's — the paper says "different," which is right, but understates how different

I read both implementations directly. The paper (§6, "Incidence and
projection-degeneracy diagnostics") says HGR's check is "a uniformly stricter
enforcement policy" while being careful that the underlying criteria are "related
but distinct." Both halves check out — but the asymmetry is sharper than the
paper's phrasing conveys, and it cuts against HGR.

**Bertini_real** (`codigo_cplusplus/src/decompositions/curve.cpp:1447`): evaluates
the randomized system's Jacobian at a **random point**, stacks the projection rows
beneath it into a square matrix `detme`, and takes its **determinant**. That is a
genuine, if numerical, **transversality test**: it asks whether the projection is
transverse to the variety's tangent spaces at a generic point — exactly the
genericity condition the method needs.

**HGR** (`src/Projection.jl:172`): checks whether ∂f′/∂x or ∂f′/∂y is the
**identically-zero polynomial**, via 3 fixed complex probes. Its own docstring is
candid that this is scoped to "exactly the crash class," deliberately skips
∂f′/∂z, and leaves positive-dimensional critical loci to downstream machinery.

So on the *criterion*, Bertini_real's is strictly stronger: it tests the actual
transversality condition; HGR tests only for a degenerate equation that would crash
the start-system builder. On the *enforcement policy*, HGR is stronger, and I
confirmed the paper's claim exactly: `curve.cpp:68` calls `br_exit(196)` on failure
— a hard abort — while `surface.cpp:409` merely prints a red warning and continues.
That inconsistency in Bertini_real is real.

**Net:** HGR enforces a weaker condition more consistently. The paper is not wrong,
but "uniformly stricter" is the phrase a reader remembers, and the criterion gap is
the thing that actually matters mathematically. Concretely: nothing in HGR checks
that the chosen projection is transverse. §2.4's almost-sure argument is carried
entirely by mathematics and by the empirical fixture record — not by a runtime
check. If pressed on "do you verify your projection is generic?", the answer is
**no, I verify it isn't catastrophically degenerate**, and that is a materially
different statement.

### 4.3 The witness point in Bertini_real's deflation is load-bearing theory, not bureaucracy — and this is the sharpest question available about §3.5

The paper (§6) presents HGR's witness-set-free deflation as a clean win:
`deflate_once` works on the bare defining system, "which we confirmed by exactly
reproducing Hauenstein and Wampler's own published deflation sequences on the
Whitney umbrella without constructing a witness set at all." The reproduction is
real (I verified the assertions). But I read Bertini_real's manual on *why* it uses
a witness point, and the reason is not overhead:

> "sub-determinants of the system are recursively added until the rank of the
> Jacobian stabilizes, where **this rank is computed using a witness point for the
> component being deflated**."

The point of a witness point is that it is a **generic point of the component**.
The rank you want is the *generic* rank on the isosingular set, not the rank at
whichever point you happened to be handed. HGR's `estimate_corank` evaluates at the
specific x₀ under examination. If x₀ is itself special *within* its own isosingular
set — a more degenerate point sitting inside a larger singular stratum — the corank
measured there can exceed the generic corank of the stratum, and `minor_size =
expected_rank − corank + 1` is then computed from a non-generic rank. The deflation
still produces *a* system; whether the resulting sequence is the one the theory
predicts is a different question.

On HGR's fixtures this doesn't bite: the Whitney tip and handle are chosen precisely
as the canonical generic representatives of their strata, which is why the published
sequences reproduce. But that means the validation is **on exactly the inputs where
the distinction is invisible** — the reproduction confirms the minor construction is
right, without exercising the generic-vs-special-point issue at all.

Suggestive supporting evidence, which I'd read as consistent with this rather than
as proof: the Taubin heart shows 5 outlier firings needing 4 rounds
(`[2,1,1,1,0]`), all at the slice-level critical points x=(±1,0) — geometrically
distinguished points. The internal notes call this "a real, non-spurious deviation."
A deeper-than-expected sequence at a geometrically special point is what you'd
expect if corank is being read at a non-generic point of its stratum.

**How to handle this on stage.** Don't volunteer it as a flaw; do have it ready,
because someone who knows the Hauenstein–Wampler paper may ask "at which point do
you compute the rank?" The good answer: "at the singular point itself, not at a
witness point of its isosingular set — which is why I can skip witness-set
construction entirely, and which is a real theoretical difference from
Bertini_real, not just an implementation shortcut. It reproduces the published
sequences on the canonical examples; whether it agrees in general when the query
point is non-generic within its own stratum is something I have not tested, and
those two claims should not be conflated." That is a strictly stronger position than
presenting witness-set-freedom as an unqualified improvement.

### 4.4 The `[3,1,1]` Whitney-handle assertion is a truncation, and the surrounding language slightly oversells it

`CLAUDE.md` describes the handle case as "`[3,1,1,...]` at the handle (0,0,1)
(genuine plateau on the 1-dim singular locus)". The paper (§3.5) says "corank
sequence [3,1,1,…] — a genuine plateau that never reaches 0."

Reading the test directly (`test/test_isosingular_deflation.jl:91-92`), the
assertion is `seq_handle == [3, 1, 1]` computed with **`max_rounds = 2`**. The
sequence is *capped*, not observed to plateau. And the reason for the cap is
recorded right there (`:86-90`): round 5 of this exact case needs ~5.3×10⁹ symbolic
minors and never terminates — it had to be killed.

The plateau is almost certainly genuine — it is what the published literature says,
and the geometry (a point on a 1-dimensional singular locus) demands it. My point
is narrower: **the test does not demonstrate it.** It demonstrates agreement with
the published sequence for two rounds. Those are different evidentiary claims, and
"never reaches 0" is a statement about all rounds that this repo's own evidence
does not reach.

There's also a subtlety in the plateau accounting worth having straight: a genuine
nonzero plateau always costs one deflation round *beyond* the terminal corank,
because `[3,1]` alone contains no repeat — you need the second `1` to see a
plateau. So the round count and the sequence length carry slightly different
information, and `rounds == 2` for the handle is consistent with terminal corank 1
reached after round 1.

Recommended phrasing: "matches the published sequence for as many rounds as is
computationally feasible — round 5 needs ~5×10⁹ symbolic minors and doesn't
terminate" is both more honest and more impressive than "never reaches 0," because
it surfaces the combinatorial wall as a finding rather than hiding it.

### 4.5 A discrepancy inside Bertini_real itself: the manual and the source disagree on the bounding-sphere radius

Independently verified, and useful to have in your pocket because it shows you read
both.

- **Bertini_real's manual**, §A.2.2: the sphere's "radius is arbitrarily chosen to
  be **3 times** the distance from the center to the furthest critical point."
- **Bertini_real's source**,
  `codigo_cplusplus/src/decompositions/decomposition.cpp:417-419`:
  `mpf_set_str(temp_rad->r,"2.0",10); … mul_mp(sphere_radius_,temp_rad,sphere_radius_);`
  with the inline comment **"double the radius to be safe."**

So the implemented factor is **2×**, not 3×. `docs/BERTINIREAL_AUDIT.md:207`
records 2× and is therefore correct against the source; the published manual is
stale or was never updated. (The 3.0 that *does* appear in the code is the fallback
radius for the empty / single-critical-point cases, `decomposition.cpp:333`,
`:345` — possibly the source of the manual's confusion.)

This is Bertini_real's discrepancy, not HGR's, and HGR uses a box anyway. Its value
is defensive: if anyone cites the manual's 3× at you, you know which is right and
why.

### 4.6 `docs/BERTINIREAL_AUDIT.md` is stale on two capabilities the project has since built — flagging, not resolving

`docs/BERTINIREAL_AUDIT.md` is dated 2026-07-15 and states as **NOT IMPLEMENTED**:

- "Random / user projection as sweep axis" (`:352`), rated "Worth adding (high) …
  the biggest 'why doesn't this do X?' vs BR", JSAG reviewer risk **High**.
- "Isosingular deflation / multiplicity split" (`:306`, `:355`) — "Singular verts
  classified via SVD; **no formal deflation**", risk **High** for singular examples.

Both are now implemented: `src/Projection.jl` in full, and `src/Solver.jl:119-633`
(`estimate_corank` at `:119` through `resolve_isosingular_dimension` at `:597`,
plus `deflate_once` `:209` and `verify_isosingular_dimension` `:419`), all exported from
`src/HomotopyGetsReal.jl:74-78`.

So that audit's two highest-rated reviewer risks have been closed, and the document
still reads as though they're open. It's gitignored local working notes, so this
has no external consequence — but it does mean **anyone (including a future agent)
using it as a current capability map will materially understate the project.** Per
`CLAUDE.md`'s ground rules I'm flagging rather than editing it; the decision about
whether to mark it historical the way `ORCHESTRATOR_BRIEFING.md` was is yours.

Related, minor, same class: the paper's §7.1 says "registered in Julia's General
registry as **v0.2.0**," while `CLAUDE.md` and `README.md` say **v0.2.1**. Worth a
pass before submission.

### 4.7 One place I'd frame the framing itself differently

Across the paper and internal docs, the recurring move is to state limitations
plainly and then explain why they're acceptable. That's the right instinct and it's
executed unusually well — §3.9's watertightness paragraph and §4.3's torus
nondeterminism paragraph are genuinely exemplary scientific writing.

The one adjustment I'd make for a *talk*, as opposed to a paper: the paper's
honesty is distributed across nine sections, so no single limitation dominates. In
a 30-minute talk with live questions, that distribution collapses — an audience
member who finds one gap will probe it and it becomes the whole conversation. The
defense is not to hide anything; it's to **pre-empt with the strongest one
yourself**. Specifically: say early, unprompted, that the critical-point solve is
the weakest link, that it fails silently on multiplicity ≥ 2, and that this is
exactly where exact methods would help — then the room's sharpest question becomes
a collaboration proposal instead of a cross-examination. §3 Q1 and Q3 are written
to be delivered that way.

---

# Appendix — numbers worth having cold

| Quantity | Value | Source |
|---|---|---|
| Test suite | 538/538 full; 478/478 fast (~3.5 min); full ~30–34 min | `CLAUDE.md`; paper §5.3 |
| Sphere residual (production, n=19,504) | mean 2.85e-8, p99 7.05e-8 | `paper_artifacts/data/results.json` |
| Ellipsoid residual (production) | mean 2.70e-7, **max 1.12e-4** | paper Table 3 |
| Taubin residual (n=1,638) | mean 1.03e-7, max 2.42e-6 | paper Table 3 |
| Torus residual | mean 2.21e-6, p99 1.008e-5 | `results.json` |
| Taubin critical z | ≈ {−1.0, 1.0, 1.0648, 1.2367} | paper §4.2 |
| Taubin degenerate slab | naive z_mid=0 rejected; z=0.06 accepted after 5 retries | paper §5.2 |
| Taubin naked edges | 188 bare → 58 snapped → 32–37 lofted (28 trials, median/mode 34) | re-measured live this session, 28 trials; paper §3.9 still says 31–35, not yet updated to match |
| HC.jl `solve()` role coverage | `Solver.jl:448` uses parameter-homotopy `solve(...; start_parameters=, target_parameters=)` to verify one flagged vertex's isosingular dimension — doesn't fit either of the paper's stated two roles | flagged live this session; paper §3 (`section3_background.tex`)'s two-role framing (zero-dim vertex solve vs. parameter-homotopy edge/face tracking) has no slot for this hybrid call; paper not updated to match |
| Whitney corank sequences | tip `[3,2,0]`; handle `[3,1,1]` (2-round cap) | `test/test_isosingular_deflation.jl:82,92` |
| Whitney round-5 cost | ~5.3×10⁹ minors, never terminates | same file, `:86-90` |
| Astroid | HGR 4 edges vs BR 6 (2 extra = BR projection artifacts) | paper §6 |
| `sample_edge` regression | max residual 0.4998 → 7.3e-7 | paper §5.1 |
| Capability survey | 24 fixtures: 14 clean, 7 caveats, 1 exception, 2 empty | `dev/scratch/capability_survey/SURVEYOR_DONE.md` |
| Deflation runtime share | ~76% of full suite (18m15s) | `docs/DESIGN_NOTES.md` |
| Key tolerances | `critical_point_tol` 1e-6; `vertex_match_tol` 1e-4; `jacobian_rank_tol` 1e-8; `singular_value_threshold` 1e-6 | `src/Config.jl:50-138` |

**Three sentences to have memorized, in case of a hostile question:**

1. "The critical-point solve is the weakest link in this pipeline: it fails
   silently on multiplicity ≥ 2, I've measured it losing astroid cusps in
   34 of 130 runs (~26%, roughly 1 in 4), and it's exactly the subproblem
   your tools solve better than mine."
2. "I don't certify anything — HomotopyContinuation.jl ships a Krawczyk-based
   certifier and I don't call it. It would certify the critical-point solve, but
   not the path-tracked faces and not the singular points, so it would strengthen
   the chain without closing it."
3. "The deflation subsystem reproduces Hauenstein–Wampler's published sequences
   exactly, and it changes nothing about the output — it's a validated diagnosis
   that isn't wired into the geometry yet."
