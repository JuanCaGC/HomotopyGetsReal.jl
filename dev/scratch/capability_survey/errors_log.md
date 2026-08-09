# Capability Survey — Errors / Discrepancies Log

This file is a running log of any exception, timeout, off-baseline result, or discrepancy
against documented predictions, observed during the 24-fixture capability survey of
HomotopyGetsReal.jl. Entries are appended in real time as fixtures run — this file is a
deliverable log, not a to-do list, and entries are never retroactively cleared or "fixed".

Baseline residual reference: sphere fixture measured median residual 2.14e-8 (confirmed live
prior to this survey run). Residuals far outside this order of magnitude are logged as
off-baseline.

---

## squircle_quartic — genuine smooth critical points misclassified Singular; edges built on off-curve Artificial placeholders instead

**Expected:** `x^4+y^4-1` is a smooth, convex, closed quartic. `decompose_1d_curve`'s x-slicing should find exactly 2 vertical-tangent critical points, at `(1,0)` and `(-1,0)` (by direct symmetry: `f_y=4y^3=0 => y=0`, then `x^4=1 => x=±1`), classified `Critical`, and used as the endpoints of the 2 edges (upper/lower half of the curve) — directly analogous to the `circle` fixture (`x^2+y^2-1`), which produced exactly this: 2 `Critical` vertices, 2 edges, both used as real edge endpoints.

**Actual:** The 2 correct points ARE found, at essentially machine-exact coordinates (`(1.0000000000000009+8e-19im, -0.0+0.0im)` and `(-1.0+3.7e-19im, 0.0+0.0im)`), but classified `v_type=Singular`, not `Critical` — `jacobian_rank=1`, `singular_values=[4.0, 0.0]` on the augmented system `[f, f_y]`. This is because `f_y=4y^3` has a triple (not simple) zero in `y` at these points, making the AUGMENTED 2x2 Jacobian rank-deficient, even though the curve itself is perfectly smooth there (`F_original`'s own gradient is `(4x^3,4y^3)=(±4,0)`, full rank 1 — not actually singular). This is a live, previously-undocumented-for-this-fixture instance of exactly the mechanism `docs/DESIGN_NOTES.md`'s "Stage 4c — `_deflation_applicable`" entry already describes in the abstract (using `f=x-y^3` as its illustrating example): a point can be `Faug`-singular while `F_original` is regular.

Worse, and NOT covered by the existing docs entry: the 2 correctly-located points are then not even used to build the edge graph. Both of the 2 resulting edges have their `left_vertex_id`/`right_vertex_id` pointing at 4 separate `Artificial`, `:endpoint_fallback`-origin vertices at `(-1,1)`, `(1,1)`, `(-1,-1)`, `(1,-1)` — none of which lie on the curve at all (`|f(1,1)| = |1+1-1| = 1`, nowhere near 0). The genuine `Singular`-typed vertices (ids 1,2) are computed but never wired into any edge's endpoint fields.

Net vertex count for a topologically simple smooth quartic: 6 (2 Singular + 4 Artificial), vs. circle's 2 (both Critical, both real edge endpoints) for the analogous topology. Residuals over the actual sampled interior curve points are still in/near baseline (median `3.8e-10`) but `p90`/`p99`/`max` (`9.6e-7`/`9.6e-7`/`9.6e-7`) run ~33x above circle's own max (`2.9e-8`), consistent with degraded Newton conditioning near the misclassified points.

**Evidence:**
```
Total vertices: 6
  id=1 type=Singular coords=ComplexF64[1.0000000000000009 + 8.296841292029852e-19im, -0.0 + 0.0im] metadata=Dict{Symbol, Any}(:tolerance_used => 1.0e-6, :singular_values => [4.000000000000011, 0.0], :jacobian_rank => 1, :cluster_member_ids => [1], :cluster_size => 1)
  id=2 type=Singular coords=ComplexF64[-1.0 + 3.676274987540704e-19im, 0.0 + 0.0im] metadata=Dict{Symbol, Any}(:tolerance_used => 1.0e-6, :singular_values => [4.0, 0.0], :jacobian_rank => 1, :cluster_member_ids => [2], :cluster_size => 1)
  id=3 type=Artificial coords=ComplexF64[-1.0 + 0.0im, 1.0 + 0.0im] metadata=Dict{Symbol, Any}(:tolerance_used => 0.0001, :origin => :endpoint_fallback)
  id=4 type=Artificial coords=ComplexF64[1.0000000000000009 + 0.0im, 1.0 + 0.0im] metadata=Dict{Symbol, Any}(:tolerance_used => 0.0001, :origin => :endpoint_fallback)
  id=5 type=Artificial coords=ComplexF64[-1.0 + 0.0im, -1.0 + 0.0im] metadata=Dict{Symbol, Any}(:tolerance_used => 0.0001, :origin => :endpoint_fallback)
  id=6 type=Artificial coords=ComplexF64[1.0000000000000009 + 0.0im, -1.0 + 0.0im] metadata=Dict{Symbol, Any}(:tolerance_used => 0.0001, :origin => :endpoint_fallback)
Total edges: 2
  id=1 left=3 right=4 n_samples=6 is_singular=false
  id=2 left=5 right=6 n_samples=6 is_singular=false
```
Residuals (over the 12 flattened `sampled_points`, which ARE Newton-projected onto the true curve): mean `3.20e-7`, median `3.81e-10`, p90 `9.60e-7`, p99 `9.60e-7`, max `9.60e-7`. Diagnostic script: `dev/scratch/capability_survey/_diag_squircle.jl`.

**Correction after visual inspection:** `renders/squircle_quartic.png` shows the actual damage is worse than the vertex/edge bookkeeping numbers alone suggest. The render shows only 2 tiny, disconnected arc fragments (near the curve's corner regions close to the bad Artificial anchor points), with the 2 correctly-located `Singular` vertices at `(+-1,0)` plotted as isolated markers with no curve segment reaching them. The bulk of the true squircle boundary -- through `(1,0)`, `(-1,0)`, and the actual top/bottom apexes near `(0,+-1)` -- is absent from the render entirely. Plausible mechanism (not confirmed further): `sample_edge`'s Newton correction, seeded from the far-off-curve Artificial endpoints `(+-1,+-1)`, converges to the NEAREST point on the true curve to each interpolated seed rather than tracing the full intended arc between the named endpoints, collapsing each edge down to a short local arc near the corner instead of spanning the claimed interval.

## folium_descartes — run-to-run nondeterminism: identical fixture/config, 1 of 3 runs silently drops a curve segment

**Expected:** `x^3+y^3-3*x*y` (folium of Descartes) is a single connected real curve through a node at the origin: a bounded loop (`t>=0` in the standard rational parametrization) plus one unbounded branch split into two arcs by the asymptote `x+y+1=0` (`t∈(-1,0)` and `t<-1`), both arcs also passing through the node. `decompose_1d_curve` should produce a fully-connected vertex/edge graph where every located vertex (the node, the loop's turning point, both bbox-exit boundary points) is used as an edge endpoint, with no curve segment dropped.

**Actual:** Running the identical fixture/config three times (once as the officially-captured survey run, twice more as immediate diagnostic re-runs) produced two different topologies due to cross-process HC.jl solver jitter (the same general nondeterminism source documented elsewhere in `docs/DESIGN_NOTES.md`, e.g. the torus's naked-edge-count spread and the 536-vs-537 test count):
- **2 of 3 runs** (the officially-captured survey run + 1 diagnostic rerun): correct, fully-connected 5-vertex/5-edge topology. All 5 vertices (1 Critical loop-turning-point, 1 Singular node, 2 Boundary bbox-exits, 1 Artificial `:endpoint_fallback` bridging vertex) referenced by at least one edge; all vertex residuals ~0 (including the Artificial one, which in this case happens to land exactly on-curve: `(x,y)=(1.5874010519681994, -2.5198420997897464)` satisfies `x^3+y^3-3xy=0` to `residual=0.0`, confirmed both numerically and algebraically -- `x=x_crit=2^{2/3}`, `y=-2*y_crit` where `y_crit=2^{1/3}`, and substituting gives `y_crit^3(y_crit^3-2) = 2*(2-2) = 0` exactly).
- **1 of 3 runs**: broken 4-edge topology. Vertex id=4 (a genuine `Boundary` vertex, `coords=(3.0254906532366572, -4.0)`, essentially exactly on-curve) is never referenced by `left_vertex_id`/`right_vertex_id` of ANY of the 4 edges -- an orphaned vertex. The curve segment that should connect it (the unbounded branch's continuation from `x=1.5874` out to `x=3.0255`, where it exits the bbox at `y=-4`) is silently absent from the edge list. No exception raised; `decompose_1d_curve` returns successfully with this segment simply missing.

The officially-captured survey JSON (`data/folium_descartes.json`) reflects the good topology (hence `outcome` left at `success_with_caveats` rather than something worse), but the fixture is demonstrated live, in this same investigation, to be capable of silently dropping a real curve segment on an unlucky run with no error signal of any kind.

**Evidence:**
Bad run (`dev/scratch/capability_survey/_diag_folium.jl`):
```
Total vertices: 5
  id=1 type=Singular coords=ComplexF64[-0.0 - 0.0im, -0.0 - 0.0im] ...
  id=2 type=Critical coords=ComplexF64[1.5874010519681994 - 3.851859888774472e-34im, 1.259921049894873 + 0.0im] ...
  id=3 type=Boundary coords=ComplexF64[-4.0 + 0.0im, 3.0254906532366572 + 0.0im] ...
  id=4 type=Boundary coords=ComplexF64[3.0254906532366572 + 0.0im, -4.0 + 0.0im] ...
  id=5 type=Artificial coords=ComplexF64[1.5874010519681994 + 0.0im, -2.5198420997897464 - 3.5718355977571093e-102im] metadata=...:origin => :endpoint_fallback
Total edges: 4
  id=1 left=3 right=1 n_samples=6 is_singular=true
  id=2 left=1 right=2 n_samples=6 is_singular=true
  id=3 left=1 right=5 n_samples=6 is_singular=true
  id=4 left=1 right=2 n_samples=6 is_singular=true
```
(vertex id=4 appears in no edge's left/right field above.)

Good run (`dev/scratch/capability_survey/_diag_folium2.jl`, with per-vertex residual + `referenced_by_edge` check added):
```
Total vertices: 5
  id=1 type=Critical coords=[1.5874010519681994, 1.259921049894873] residual=0.0 referenced_by_edge=true
  id=2 type=Singular coords=[0.0, 0.0] residual=0.0 referenced_by_edge=true
  id=3 type=Boundary coords=[-4.0, 3.025490653236657] residual=2.1316282072803006e-14 referenced_by_edge=true
  id=4 type=Boundary coords=[3.0254906532366572, -4.0] residual=0.0 referenced_by_edge=true
  id=5 type=Artificial coords=[1.5874010519681994, -2.5198420997897464] residual=0.0 referenced_by_edge=true
Total edges: 5
  id=1 left=3 right=2 n_samples=6
  id=2 left=2 right=1 n_samples=6
  id=3 left=2 right=1 n_samples=6
  id=4 left=2 right=5 n_samples=6
  id=5 left=5 right=4 n_samples=6
```

## three_concurrent_lines_reducible — vertical line component entirely absent from edge graph; triple point misclassified Artificial not Singular

**Expected:** `x*(x-y)*(x+y)` is 3 lines through the origin (`x=0`, `x=y`, `x=-y`), a triple point at the origin (multiplicity 3). All 3 lines should appear in the output, each split into 2 edges by the origin, with the origin itself as a `Singular` vertex shared by all 3 (6 edges, 1 Singular + 6 Boundary vertices total, by direct analogy with `tangent_parabolas_reducible`'s clean handling of its own reducible tangent singularity 2 fixtures earlier in this same survey).

**Actual:** Only 2 of the 3 lines (`x=y` and `x=-y`) appear in the output as 4 edges, all correctly stitched through the origin. The 3rd line, `x=0` (the one line exactly VERTICAL in the x-slicing sense this pipeline uses), has both its boundary vertices found correctly (`(0,4)`/`(0,-4)`, exactly on-curve, residual `0.0`) but **never referenced by any edge** -- confirmed via explicit `referenced_by_edge` check, both show `false`. Root cause is structural, not a stray bug: `decompose_1d_curve` builds edges by sampling at x-values strictly BETWEEN pairs of `distinct_xs`; a component that exists at exactly one x-value and nowhere else (a vertical line) can never be hit by any such sample, so it is invisible to this pipeline's edge-construction step by construction -- it can only ever appear as isolated, disconnected vertices. This is a real, previously-undocumented (per this survey's own grep of `docs/DESIGN_NOTES.md`) capability boundary: `decompose_1d_curve` cannot represent a genuinely vertical-line real component.

Separately, and consistent with the broader multiplicity>=2 pattern already documented for nodes/cusps/astroid-cusps: the triple point at the origin (multiplicity 3, more degenerate than any of those) is found and correctly used to stitch the 2 captured lines together, but classified `v_type=Artificial` (`:endpoint_fallback` origin), not `Singular` -- 0 `Singular` vertices in the whole result, despite an obvious, unambiguous algebraic singularity being present and located.

**Evidence:**
```
Total vertices: 7
  id=1 type=Boundary coords=[-4.0, 4.0] residual=0.0 referenced_by_edge=true metadata=...cluster_size => 2)   [corner, x=-y line]
  id=2 type=Boundary coords=[-4.0, -4.0] residual=0.0 referenced_by_edge=true metadata=...cluster_size => 2)  [corner, x=y line]
  id=3 type=Boundary coords=[4.0, 4.0] residual=0.0 referenced_by_edge=true metadata=...cluster_size => 2)    [corner, x=y line]
  id=4 type=Boundary coords=[4.0, -4.0] residual=0.0 referenced_by_edge=true metadata=...cluster_size => 2)   [corner, x=-y line]
  id=5 type=Boundary coords=[0.0, -4.0] residual=0.0 referenced_by_edge=false metadata=...:fixed_variable => :y...   [x=0 line -- ORPHANED]
  id=10 type=Boundary coords=[0.0, 4.0] residual=0.0 referenced_by_edge=false metadata=...:fixed_variable => :y...   [x=0 line -- ORPHANED]
  id=11 type=Artificial coords=[0.0, -1.440944494681213e-14] residual=0.0 referenced_by_edge=true metadata=Dict{Symbol, Any}(:tolerance_used => 0.0001, :origin => :endpoint_fallback)   [triple point -- should be Singular]
Total edges: 4
  id=1 left=2 right=11 n_samples=6   [x=y, lower-left arm]
  id=2 left=1 right=11 n_samples=6   [x=-y, upper-left arm]
  id=3 left=11 right=3 n_samples=6   [x=y, upper-right arm]
  id=4 left=11 right=4 n_samples=6   [x=-y, lower-right arm]
```
Diagnostic script: `dev/scratch/capability_survey/_diag_threelines.jl`.

## ellipsoid — max residual ~170x above p90, small number of outlier mesh points

**Expected:** sphere/ellipsoid baseline is ~2e-8 median with residuals in a tight band (sphere itself, this same survey: median `2.14e-8`, max `5.61e-8`, less than a 3x spread top-to-bottom).

**Actual:** `x^2+4*y^2+9*z^2-1` (asymmetric ellipsoid): median `2.573e-8` and p90 `5.975e-8` sit right in the baseline band, but p99/max both jump to `1.0003e-5` -- ~170x above p90, ~180x above sphere's own max. Since `n_points=160`, this is roughly 1-2 outlier mesh points, not a broad degradation. Outcome still clean (no exception), vertex/face/edge counts as expected (2 Critical, 2 faces, 2 edges).

**Evidence:**
```
residuals: mean=5.259e-7, median=2.573e-8, p90=5.975e-8, p99=1.0002641087719644e-5, max=1.0002641210538066e-5, n_points=160
```
Plausibly connected to `docs/DESIGN_NOTES.md`'s own "Adaptive re-anchoring" / "The ellipsoid discovery" entry (Face tracking, Phase 5): a fixed literal patch line is only guaranteed to keep intersecting the level curve as `z` sweeps if the surface's gradient is radially symmetric about the slicing axis -- true for a sphere, false in general for an ellipsoid, with the specific failure mode being loss of transversal intersection nearer a pole. Not confirmed as the same mechanism here (not investigated further, per this survey's log-don't-diagnose scope), but consistent with it.

## whitney_umbrella — CONFIRMATION of documented `OverflowError` prediction (not a new finding)

**Expected:** `docs/DESIGN_NOTES.md` ("Uncaught exception when `compute_critical_z_slices` finds zero critical z-values and the naive bbox midpoint is itself degenerate") documents that a bare, default (`deflate=false`) `decompose_3d_surface` call on `x^2-y^2*z` throws `OverflowError: Cannot compute a start system`, because `compute_critical_z_slices` finds no critical z, the naive bbox-midpoint slab lands exactly on the degenerate `z=0` axis, and `slice_at_z(F,0.0,cfg)` throws with no retry.

**Actual:** Reproduced exactly, live, this run. `OverflowError: Cannot compute a start system.`, raised inside `HomotopyContinuation.polyhedral`'s `PolyhedralStartSolutionsIterator`, reached via `intersect_bounding_object` -> `decompose_1d_curve` -> `slice_at_z` -> `_robust_slice_at_z` -> `decompose_3d_surface`. **This confirms the existing documented entry; it is not a new discrepancy.** Logged per this survey's instructions to record confirmations, not only new findings.

**Evidence:** Full stacktrace captured in `dev/scratch/capability_survey/data/whitney_umbrella.log` / `.json` `error` field; terminal frame: `intersect_bounding_object(F::System, cfg::HomotopyConfig{Float64}; deflate::Bool, F_original::System) @ HomotopyGetsReal ~/HomotopyGetsReal/src/Solver.jl:873`.

## horn_torus — silent completely-empty decomposition, extending the documented "cone" prediction to a NEW fixture

**Expected:** No existing docs entry covers the horn torus specifically. But `docs/DESIGN_NOTES.md`'s backlog item 3/4 ("HC.jl polyhedral-solve reliability...") documents the SAME mechanism for the cone (fixture 22, "untested" per that entry, checked live later in this same survey): a symmetric `bbox_z` puts the naive whole-bbox-midpoint slab's `z_mid` exactly on an undetected degenerate/singular z, and `compute_critical_z_slices` returning empty for that z means `decompose_3d_surface` silently returns a completely empty decomposition, no crash.

**Actual:** `(x^2+y^2+z^2)^2 - 4*(x^2+y^2)` (a horn torus -- degenerate case of a torus where the inner "hole" radius shrinks to 0, producing a genuine self-tangency/pinch point at the origin) reproduces this exact silent-empty pattern live, on a fixture not previously tried against this codebase (per this survey's own grep of `docs/DESIGN_NOTES.md`, `test/`, `dev/scratch/`): `decompose_3d_surface(F, cfg)` (bare, default call, no `deflate`/`projection`) returns `(vertices=[], edges=[], faces=[], mesh)` with all counts 0, no exception. Root cause plausibly the same structural mechanism as the cone/Whitney cases: at `z=0`, the surface equation factors as `rho^2(rho^2-4)=0` (`rho^2=x^2+y^2`), i.e. the origin (the pinch point, a genuine singularity) UNION the circle `rho=2` -- and the default `bbox_z=(-4,4)` is symmetric, so the naive whole-bbox slab lands its `z_mid` exactly at this degenerate `z=0`, same as the cone's apex and Whitney's handle axis.

**Evidence:**
```
vertex_counts: {Critical: 0, Boundary: 0, Singular: 0, Artificial: 0}
edge_count: 0, face_count: 0, mesh_vertex_count: 0, mesh_triangle_count: 0
outcome: no exception raised, wall_time_seconds: 16.66
```

## horn_torus — separate finding: `plot_surface_decomposition` throws on an empty mesh instead of handling it gracefully

**Expected:** Given `decompose_3d_surface` can (as just documented above, and as documented for the cone) legitimately return a completely empty mesh with no exception, a caller following the README's own quick-start pattern (`plot_surface_decomposition(mesh; color_by=:z, cfg=cfg)` immediately after `decompose_3d_surface`) would reasonably expect either a valid (if trivial/empty) figure, or at worst a clear, documented error -- not an uncaught internal `ArgumentError` surfaced from deep inside a colorbar-scaling helper.

**Actual:** Calling `plot_surface_decomposition(mesh; color_by=:z, cfg=cfg)` on the horn torus's empty mesh throws `ArgumentError: reducing over an empty collection is not allowed; consider supplying init to the reducer`, from `Base.extrema` inside `HomotopyGetsReal._near_constant_colorrange` (`src/Visuals.jl:210`), called from `plot_surface_decomposition` (`src/Visuals.jl:282`). This is a distinct, separate robustness gap from the empty-decomposition finding above -- it's specifically that the PLOTTING path has no empty-mesh guard, not that the decomposition itself is empty (which can be a legitimate, if uninformative, result).

**Evidence:**
```
ArgumentError: reducing over an empty collection is not allowed; consider supplying `init` to the reducer
Stacktrace:
  [1] _empty_reduce_error() @ Base ./reduce.jl:311
  ...
  [12] _near_constant_colorrange(colorvals::Vector{Float64}) @ HomotopyGetsReal ~/HomotopyGetsReal/src/Visuals.jl:210
  [13] plot_surface_decomposition(mesh::GeometryBasics.Mesh{...}; color_by::Symbol, ...) @ HomotopyGetsReal ~/HomotopyGetsReal/src/Visuals.jl:282
```
Full text in `dev/scratch/capability_survey/data/horn_torus.log`.

## cone — CONFIRMATION of documented "untested" silent-empty prediction (fixture 22, high priority)

**Expected:** `docs/DESIGN_NOTES.md`'s backlog item 4 states, explicitly flagged as "untested" / "the first live check of it": for `x^2+y^2-z^2` (the cone), `compute_critical_z_slices` returns `Float64[]` (the apex, `x=y=0, z^2=0`, is a genuine isolated critical point mathematically but not detected), and because `bbox_z` is symmetric, the resulting single whole-bbox slab's naive midpoint lands exactly on the undetected apex, causing `decompose_3d_surface` to silently return a **completely empty** decomposition (0 vertices, 0 faces) -- no crash, unlike Whitney.

**Actual:** Reproduced exactly, live, bare/default call (no `deflate`/`projection`/hand-tuned bbox): `decompose_3d_surface(F, cfg)` returns `(vertices=[], edges=[], faces=[], mesh)` with ALL counts 0 (`mesh_vertex_count=0`, `mesh_triangle_count=0`), no exception, wall time 16.77s. **This is the first live confirmation of this previously-untested backlog entry.** Also hit the same separate `plot_surface_decomposition` empty-mesh `ArgumentError` documented above for `horn_torus` (same root cause, same stack).

**Evidence:**
```
vertex_counts: {Critical: 0, Boundary: 0, Singular: 0, Artificial: 0}
edge_count: 0, face_count: 0, mesh_vertex_count: 0, mesh_triangle_count: 0
error: null, outcome: no exception, wall_time_seconds: 16.77
```

## quartic_superellipsoid — surface-level recurrence of the squircle_quartic mechanism, worse: Artificial edge endpoints residual ~1.0 (nowhere near surface)

**Expected:** `x^4+y^4+z^4-1` is smooth and convex (a "superellipsoid"). By direct analogy with `sphere` (2 Critical vertices at the x-extrema, 2 edges/faces, both used as real edge endpoints), `decompose_3d_surface`'s crit-z-slice vertex/edge bookkeeping should find the 2 genuine x-extremal critical points (`(+-1,0,0)`) classified `Critical`, wired as real edge endpoints.

**Actual:** Same root mechanism as the `squircle_quartic` curve finding earlier in this log (`f_x=4x^3` -- and by symmetry `f_y`, `f_z` -- has a higher-order zero at the extremal points, making the augmented critical-point system's Jacobian rank-deficient even though the surface itself is perfectly smooth there): the 2 genuine critical points at `(1,0,0)` and `(-1,0,0)` ARE found, at essentially exact coordinates (residual `2.66e-15` and `4.44e-16`), but classified `v_type=Singular` (`jacobian_rank=1`, `singular_values=[4.0, 0.0]`), and **orphaned from the edge graph** (`referenced_by_edge=false` for both). Worse than the curve case: the 4 `Artificial` `:endpoint_fallback` vertices used as the ACTUAL edge endpoints have residual `~0.9999999999999996` -- i.e. essentially exactly `1.0` off the surface (`f(-1,1,0) = 1+1+0-1 = 1`), not merely "approximately near" like some curve-case fallbacks.

**Correction after visual inspection (this entry was updated after first being written from JSON/diagnostic numbers alone):** `renders/quartic_superellipsoid.png` shows a visible wedge-shaped notch/gouge cut into the surface near one pole -- the mesh geometry itself IS visibly affected, not merely the top-level `vertices`/`edges` bookkeeping as the aggregate residual stats (204 vertices/320 triangles, median `5.26e-8`) suggested in isolation. Contrast with `renders/sphere.png` and `renders/hyperboloid_one_sheet.png` (same pipeline, both fully closed/smooth) confirms this is a real, visible defect specific to this fixture, not a rendering artifact. This matches the empirical-evidence-standard lesson (verify with the actual render, not just aggregate stats) -- the residual summary alone under-reported the severity here.

**Evidence:**
```
Total vertices: 6
  id=1 type=Singular coords=[0.9999999999999993, -0.0, -3.3306690738754696e-16] residual=2.6645352591003757e-15 referenced_by_edge=false metadata=Dict{Symbol, Any}(:tolerance_used => 1.0e-6, :singular_values => [3.999999999999992, 0.0], :jacobian_rank => 1, :cluster_member_ids => [1], :cluster_size => 1)
  id=2 type=Singular coords=[-0.9999999999999999, -0.0, -3.3306690738754696e-16] residual=4.440892098500626e-16 referenced_by_edge=false metadata=...jacobian_rank => 1...
  id=3 type=Artificial coords=[-0.9999999999999999, 1.0, -3.3306690738754696e-16] residual=0.9999999999999996 referenced_by_edge=true metadata=...:origin => :endpoint_fallback
  id=4 type=Artificial coords=[0.9999999999999993, 1.0, -3.3306690738754696e-16] residual=0.9999999999999973 referenced_by_edge=true
  id=5 type=Artificial coords=[-0.9999999999999999, -1.0, -3.3306690738754696e-16] residual=0.9999999999999996 referenced_by_edge=true
  id=6 type=Artificial coords=[0.9999999999999993, -1.0, -3.3306690738754696e-16] residual=0.9999999999999973 referenced_by_edge=true
Total edges: 2
  id=1 left=3 right=4 n_samples=6
  id=2 left=5 right=6 n_samples=6
Total faces: 2
```
Diagnostic script: `dev/scratch/capability_survey/_diag_superellipsoid.jl`.

## 2026-08-06 follow-up — render-count accounting, and the torus/taubin_heart/hyperboloid_one_sheet visual gaps Juan caught by eye

Targeted follow-up, not a re-run of the 24-fixture survey. Triggered by
Juan looking directly at `torus.png`, `taubin_heart.png`, and
`hyperboloid_one_sheet.png` and catching visible defects the original
survey never called out as findings.

### Step 1 — render-count accounting (documentation gap, not a new bug)

**One-line summary:** 20 of 24 renders were produced. The 4 missing:
`whitney_umbrella` (exception before the plotting step is ever reached —
no PNG expected, not a gap); `cone`/`horn_torus`/`empty_surface` (all
three: `decompose_3d_surface` returns a completely empty 0/0/0/0
decomposition, `plot_surface_decomposition` is then attempted and throws
`ArgumentError: reducing over an empty collection is not allowed` from
`_near_constant_colorrange`, `src/Visuals.jl:210`).

**Confirmed against the raw evidence, not just re-asserted:** every one
of these 4 already has a `"PLOT FAILED: ..."` (or, for `whitney_umbrella`,
no plot attempt at all — `outcome=exception` before reaching that code
path) recorded in its own `data/<fixture>.json` `notes` field, and the
`cone`/`horn_torus` instances of the `_near_constant_colorrange` crash are
already logged above (see the `cone` and `horn_torus` entries in this
same file). This was a documentation-completeness gap in
`SURVEYOR_DONE.md`/`summary_report.md` (neither stated the 20/24 count
explicitly), not an unlogged or new discrepancy — nothing here contradicts
anything already recorded.

### Step 2a — torus and taubin_heart: CONFIRMS the known no-incidence naked-edge mechanism

Neither fixture in the original 24-fixture survey was run with
`incidence=true` (confirmed: `run_fixture.jl`'s fixture registry never
passes it for any fixture). `docs/DESIGN_NOTES.md`'s "Watertightness
measurements" entry already documents this exact mechanism for Taubin
heart at (presumably) production/default density: bare `weld_mesh`
leaves 188 naked edges, Phase 9b incidence-stitching closes to 58, full
Phase 9c coordinated loft (same `incidence=true` path) reduces further to
31-35 — with the entry's own text stating fold/point-type boundaries
close COMPLETELY while multi-face edge-type boundaries (singular notch,
saddle pair) reduce substantially WITHOUT fully closing.

Re-ran both fixtures at this survey's own coarse config
(`edge_sample_density=6, midslice_sample_density=8`; torus also with its
established `projection=:random, rng=Xoshiro(42)`), bare and then with
`incidence=true`, and counted naked edges directly via
`HomotopyGetsReal._naked_mesh_edges(mesh)` (the same utility
`docs/DESIGN_NOTES.md`'s own entry uses):

| Fixture | naked edges (bare) | naked edges (incidence=true) | reduction | wall time bare / incidence |
|---|---|---|---|---|
| torus | 80 | 8 | 10.0x | 48.4s / 4.3s |
| taubin_heart | 132 | 29 | 4.6x | 17.7s / 11.6s |

**Visual confirmation, not just the numbers**: `taubin_heart_incidence.png`
shows the top cleft gap from `taubin_heart_bare.png`/the original
survey's `taubin_heart.png` largely closed, matching the documented
"fold/point-type boundaries close COMPLETELY" pattern. `torus_incidence.png`
shows visibly fewer/smaller dark gap regions than `torus_bare.png`/the
original survey's `torus.png`, but **does still show residual visible
gaps — it is not fully closed**, matching the documented "multi-face
edge-type boundaries... reduce substantially WITHOUT fully closing" caveat
exactly. Do not read `torus_incidence.png` as "clean, gap-free" — 8
naked edges remain, and they are visible.

**Explicitly not claimed**: this survey's coarse-density counts (132→29
for taubin_heart) are NOT being equated with the documented production-
density baseline (188→58→31-35) — different `edge_sample_density`/
`midslice_sample_density` produce different absolute triangle/edge
counts, so the absolute numbers aren't comparable across those two
configs. (29 happens to land inside the documented 31-35 range for the
production-density Phase 9c result — this is flagged explicitly as
coincidental-looking, not claimed as a meaningful correspondence, given
the density mismatch.) What IS directly comparable and confirmed: the
qualitative mechanism (bare leaves substantial naked edges;
`incidence=true` reduces them by roughly 5-10x but does not zero them
out) reproduces exactly at this survey's own coarse density too.

Renders saved: `renders/torus_bare.png`, `renders/torus_incidence.png`,
`renders/taubin_heart_bare.png`, `renders/taubin_heart_incidence.png`.
Full run log: `data/_investigate_render_gaps.log`. Script:
`dev/scratch/capability_survey/investigate_render_gaps.jl`.

### Step 2b — hyperboloid_one_sheet: CONFIRMS hypothesis (b), a coarse-sampling rim artifact, NOT the naked-edge/incidence mechanism above

Re-ran at `HomotopyConfig{Float64}()` full defaults (`edge_sample_density=50,
midslice_sample_density=100`, vs. this survey's coarse `6`/`8`) and
compared directly against the original survey's coarse-density render.

**Visual result**: the coarse-density render (`renders/hyperboloid_one_sheet.png`,
from the original survey) shows a large, blocky, angular zigzag at the
bbox-clipped top rim (`z=4`) — exactly what Juan flagged. The
default-density render (`renders/hyperboloid_one_sheet_default_density.png`)
shows that same rim as smooth and regular, with only fine sub-pixel
triangulation texture remaining — no large angular artifacts. **The rim
smooths out at higher density — this points to (b), a coarse-sampling
artifact specific to the bbox-clipped boundary, not (a) the same
naked-edge/incidence mechanism as torus/taubin, and not (c) something
else.**

**A category-error caveat, checked explicitly rather than assumed**: naked
edges at default density measured at 196 — numerically *larger* than
either torus's or taubin_heart's *bare* counts above. This number is
flagged as **not meaningfully comparable** to the torus/taubin figures:
for an unbounded surface clipped by the bounding box, every edge along
that clipped rim is naked *by construction* — there is no surface beyond
the box for it to stitch to, unlike torus/taubin's naked edges, which are
genuine interior stitching gaps between adjacent faces. A high naked-edge
count here does not indicate the same defect class; it's an expected
property of any bbox-clipped unbounded-surface mesh, regardless of
sampling density. This was not re-measured at survey/coarse density since
the raw count isn't the informative signal for this fixture class — the
render comparison is.

Render saved: `renders/hyperboloid_one_sheet_default_density.png`. Same
run log/script as Step 2a.
