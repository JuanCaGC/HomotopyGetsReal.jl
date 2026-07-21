# src/SurfaceDecomposition.jl
#
# Phase 5 (surface-sweeping half): everything expressible in terms of
# EXISTING Phase 2/3 primitives, with no new tracking machinery of its
# own (that lives in FaceTracking.jl, included before this file -- see
# HomotopyGetsReal.jl's include order note):
#   - compute_critical_z_slices: thin wrapper around
#     Solver.compute_critical_points's existing 3-variable branch.
#   - slice_at_z: literal z-substitution + 3D-lifting adapter around
#     Topology.decompose_1d_curve, called AS-IS (confirmed in the Phase 5
#     investigation: zero signature/config changes needed there).
#   - decompose_3d_surface: top-level orchestrator, mirroring
#     decompose_1d_curve's own "the orchestrator runs every pipeline step,
#     including final resampling/welding, so callers never have to
#     remember a manual follow-up call" precedent (decompose_1d_curve
#     itself calls sample_edge internally before returning).
#   - weld_mesh: final cross-face vertex welding + triangle remapping +
#     winding-consistency correction, producing a watertight
#     GeometryBasics.Mesh.
#
# Precision boundary: identical to every prior phase. slice_at_z's
# `subs(f, z_var => Float64(z_val))` mirrors compute_midslice's own
# Float64-substitution-before-solve pattern (HC's polyhedral start
# system has no Complex{BigFloat} method); weld_mesh welds at genuine
# T-precision (NOT truncated to Float32 the way the old prototype's
# weld_faces_to_mesh did at prototipo_viejo_julia/SurfaceTopology.jl:181,
# which silently discarded BigFloat precision before the dedup decision
# itself) and only narrows to Float32/Point3f at the very last step, for
# the GeometryBasics.Mesh container itself.

"""
    compute_critical_z_slices(F::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> Vector{T}

Find critical z-values where the surface cannot be written as z = g(x, y).

`F` must be one equation in three variables ordered `[x, y, z]` (z last).
Delegates to [`compute_critical_points`](@ref) and clusters the resulting
z-coordinates with `cfg.vertex_match_tol`.

# Arguments
- `F::System`: surface system.
- `cfg::HomotopyConfig{T}`: tolerances for solving and clustering.

# Returns
Sorted `Vector{T}` of distinct z-values.
"""
function compute_critical_z_slices(F::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    length(F.variables) == 3 && length(F.expressions) == 1 || throw(ArgumentError(
        "compute_critical_z_slices: expected F with exactly 1 equation in exactly 3 variables " *
        "(a raw surface, variables ordered [x_var, y_var, z_var] -- z LAST); " *
        "got $(length(F.expressions)) equation(s) in $(length(F.variables)) variable(s).",
    ))
    # compute_critical_points differentiates w.r.t. the first two variables only;
    # z must be last (same convention as build_patch_system).
    crit_vertices = compute_critical_points(F, cfg)
    z_values = T[real(v.coordinates[3]) for v in crit_vertices]
    # Same cluster_scalars / vertex_match_tol as decompose_1d_curve uses for x-values.
    return cluster_scalars(z_values, cfg.vertex_match_tol)
end

"""
    slice_at_z(F::System, z_val::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> (vertices_3d::Vector{NativeVertex{T}}, edges_3d::Vector{Edge{T}})

Decompose the plane curve f(x, y, z_val) = 0 and lift it to 3D at fixed z.

Substitutes `z => z_val`, runs [`decompose_1d_curve`](@ref) on the resulting
2-variable system, and appends `z_val` as the third coordinate of every vertex
and edge sample. Vertex and edge ids are local to this call; renumber before
combining slices (as [`decompose_3d_surface`](@ref) does).

# Arguments
- `F::System`: surface system (`length(F.expressions) == 1`, `nvariables(F) == 3`).
- `z_val::T`: slice height.
- `cfg::HomotopyConfig{T}`: curve decomposition settings.

# Returns
A tuple of 3D vertices and edges at the given z.
"""
function slice_at_z(F::System, z_val::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    length(F.variables) == 3 && length(F.expressions) == 1 || throw(ArgumentError(
        "slice_at_z: expected F with exactly 1 equation in exactly 3 variables " *
        "(a raw surface); got $(length(F.expressions)) equation(s) in $(length(F.variables)) variable(s).",
    ))
    x_var, y_var, z_var = F.variables
    f = F.expressions[1]
    # Float64 substitution matches compute_midslice (HC polyhedral start has no Complex{BigFloat} path).
    f_2d = subs(f, z_var => Float64(z_val))
    F_2d = System([f_2d], variables = [x_var, y_var])

    vertices_2d, edges_2d = decompose_1d_curve(F_2d, cfg)

    # decompose_1d_curve restarts ids near 1 each call; callers must offset before concatenating.
    vertices_3d = NativeVertex{T}[
        NativeVertex{T}(
            id = v.id,
            coordinates = vcat(v.coordinates, Complex{T}(z_val)),
            v_type = v.v_type,
            metadata = v.metadata,
        )
        for v in vertices_2d
    ]
    edges_3d = Edge{T}[
        Edge{T}(
            id = e.id,
            left_vertex_id = e.left_vertex_id,
            right_vertex_id = e.right_vertex_id,
            sampled_points = [vcat(p, T(z_val)) for p in e.sampled_points],
            is_singular = e.is_singular,
        )
        for e in edges_2d
    ]
    return vertices_3d, edges_3d
end

"""
    _robust_slice_at_z(F::System, patch::NamedTuple, z_bottom::T, z_top::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> (vertices_3d::Vector{NativeVertex{T}}, edges_3d::Vector{Edge{T}}, z_mid::T)

Chooses which literal `z` to hand to [`slice_at_z`](@ref) for one slab
`(z_bottom, z_top)`, defending against a genuine failure mode discovered
while validating this file against the Taubin heart surface
(`(x^2+(1.2y)^2+z^2-1)^3 - x^2z^3 - 0.1(1.2y)^2z^3`,
`scratch_phase5_taubin_check.jl` section 6): the naive exact midpoint
`(z_bottom+z_top)/2` can coincidentally land on a z-value at which
`f(x,y,z_mid)` is itself a NON-REDUCED (repeated-factor) plane curve --
for that surface, `z_mid=0` makes every z^3 term vanish, leaving
`f(x,y,0) = (x^2+1.44y^2-1)^3`, an exact cube. [`Topology.decompose_1d_curve`](@ref)
has no way to know this ahead of time; what actually happens is that its
true critical points get correctly classified `Singular` (their Jacobian
IS genuinely rank-deficient there), but `connect_the_dots!` then cannot
track paths back onto those `Singular` vertices within `vertex_match_tol`,
so `Topology._resolve_endpoint`'s fallback fabricates new `Artificial`
vertices wherever the paths actually (wrongly) landed -- silently
producing edges that do not follow the true curve at all.

`_resolve_endpoint` tags exactly this fallback provenance with
`metadata[:origin] = :endpoint_fallback` (see its own docstring for why
this is a DIFFERENT, more specific signal than plain `v_type ==
Artificial`, which can also arise from a benign `cluster_vertices`
merge). This function's retry loop keys off that tag, but NOT alone --
see the "why not just Artificial" note below.

1. Try `z_mid = (z_bottom+z_top)/2` (the normal, cheap path -- this is
   the ONLY attempt made for the overwhelming majority of slabs, e.g.
   every sphere/ellipsoid slab and 3 of this Taubin heart's own 5 slabs
   in the regression suite; empty slices, which have no vertices at all,
   trivially pass this check with zero retries).
2. If any returned vertex has `v_type == Artificial &&
   metadata[:origin] == :endpoint_fallback` **AND** at least one
   returned vertex has `v_type == Singular`, retry at a perturbed z_mid:
   for attempt `k = 1, 2, ..., cfg.max_z_mid_retries`, offset
   `= min(ceil(k/2), 0.45/cfg.z_mid_retry_frac) * cfg.z_mid_retry_frac *
   (z_top-z_bottom)` with alternating sign `(-1)^(k+1)` -- i.e. tries
   `+1, -1, +2, -2, ...` step multiples of `cfg.z_mid_retry_frac *
   (z_top-z_bottom)`, small steps first, alternating direction before
   growing. The `0.45` cap is a hardcoded SAFETY BOUND, not a new config
   knob: it just guarantees the perturbed `z_mid` can never reach within
   5% of `z_bottom`/`z_top` themselves (the slab's own critical/bbox
   boundaries, which may be singular in their own right), independent of
   how `cfg.z_mid_retry_frac`/`cfg.max_z_mid_retries` are tuned.
3. Returns as soon as an attempt comes back clean (empirically, on the
   Taubin heart's degenerate slab, attempt 1 already succeeds -- see
   `scratch_phase5_taubin_check.jl`'s retry report). If EVERY attempt
   (the original midpoint plus all `cfg.max_z_mid_retries` perturbations)
   remains suspect, throws an `ErrorException` naming the slab bounds and
   attempt count -- a loud, explicit failure for a genuinely pathological
   slab, never a silent fallback to a slice already known to be wrong.

# Why "`:endpoint_fallback` AND `Singular`", not just `:endpoint_fallback`
The Taubin heart's `[1.0, 1.0648]` slab (a narrow slab immediately below
the singular notch, where the level curve genuinely has 4 x-critical
points -- a small inner sliver plus an outer loop, not a simple oval)
ALSO produces `:endpoint_fallback` vertices, from legitimately more
complex (not degenerate) curve topology: `connect_the_dots!` needs an
extra vertex to stitch together a curve with more than 2 x-extrema, and
none of that slab's critical points are ever classified `Singular` (they
are all ordinary, well-conditioned `Critical` points). Retrying on
`:endpoint_fallback` ALONE made this slab impossible to resolve (every
perturbed `z_mid` nearby has the same 4-x-critical-point topology, so the
retry loop exhausts `cfg.max_z_mid_retries` and throws) even though this
slab's SWEPT output was already fine before this fix existed (measured
max `|f|` along swept points: `2.4e-6`, via `track_face`'s own
`_project_to_slice` Gauss-Newton correction, which converges reliably
here precisely because the local gradient is well-conditioned -- unlike
the true repeated-factor case, where the gradient vanishes identically
along the whole degenerate curve and Newton correction has nothing to
converge on).

A direct residual-magnitude tiebreaker (flag `z_mid` suspect if
`maximum(|_residual_at(patch, p, cfg)| for p in every edge sample point)`
exceeds some threshold, instead of inspecting vertex types) was
empirically tested and rejected: raw `sample_edge` output is LINEARLY
INTERPOLATED between homotopy-tracked points (a known, pre-existing
approximation -- see `Topology.sample_edge`'s own docstring, and exactly
why `track_face` already needs `_project_to_slice`'s correction), so its
residual has a substantial baseline even on completely healthy slices:
a fully clean control slice at `z=0.05` (2 vertices, 0 `Artificial` of
any kind) already measures max raw residual `0.105`; the fully clean
`[1.0648,1.2367]` slab measures `0.144`; the narrow-but-fine
`[1.0,1.0648]` slab measures `0.268`; and the genuinely-degenerate
`z_mid=0` slab measures `1.000`. There is no scale-free gap to split a
threshold into -- the "fine" cases already span nearly 3x among
themselves, and the "bad" case is only ~3.7x above the worst "fine" case,
entirely because of shape/curvature-dependent interpolation error
unrelated to whether `f(x,y,z_mid)` is actually non-reduced. The
`Singular`-typed-vertex co-occurrence check has no such scale problem
(it is a boolean, not a magnitude) and correctly separates all four cases
above (`true` only for the genuinely degenerate slab).

# The vertex-type gate is necessary but NOT sufficient: the gradient gate
Confirmed by re-running the very same Taubin heart regression once the
`Singular` co-occurrence refinement above was in place: the `[-1,1]`
slab's retry now lands on a topologically CLEAN `z_mid=0.02` (2 edges, 0
`Artificial`/`Singular` vertices) -- yet the downstream sweep was STILL
catastrophic (`track_face`'s max `|f|` along swept points: `1.43`, worse
than before this fix). Root cause: near `z=0`, `f(x,y,z) ≈ g(x,y)^3 -
z^3·h(x,y)` for this surface, so along the level curve (where `g≈0` by
definition of it being close to the true repeated-factor curve) EVERY
first partial derivative of `f` -- not just the ones `decompose_1d_curve`
happens to probe -- scales like `O(z^2)`: genuinely tiny, but not
*exactly* zero, so no vertex ever gets classified `Singular` and the
vertex-type gate alone cannot see it. This wrecks `FaceTracking`'s
Newton-based patch tracking (tiny gradient magnitude means the patch
system is nearly singular) even though the 2D curve decomposition itself
looks perfectly clean.

Two normalizations were tried to close this gap, reusing
[`FaceTracking.patch_direction`](@ref) (the exact `(a,b) =
(f_y,-f_x)` pair `track_face` itself seeds each sweep with) rather than
introducing new gradient-computation logic:

- **Same-point ratio** `hypot(patch_direction(patch,x,y,z,cfg)...) /
  |f_z(x,y,z)|` at the candidate's own anchor: REJECTED. Measured on the
  Taubin heart, this ratio does not correlate with sweep quality at all
  -- the two already-healthy reference slabs measured `0.60` and `1.16`,
  while several already-CONFIRMED-BAD candidates near `z=0` measured
  `2.08`-`4.29`, i.e. HIGHER than the healthy baseline. A same-point
  ratio cannot discriminate this failure mode because near `z=0`, `f_z`
  is suppressed by the SAME `O(z^2)` factor as `f_x`/`f_y` (all partials
  vanish together), so the ratio between them stays `O(1)` right through
  the degenerate neighborhood.
- **Cross-z reference ratio** (adopted): compare a candidate `z_mid`'s
  own anchor-gradient magnitude against a reference magnitude measured
  elsewhere ON THE SAME SLAB, at two FIXED locations independent of the
  retry ladder's own step schedule: the slab's quarter-points,
  `z_bottom + 0.25*(z_top-z_bottom)` and `z_top - 0.25*(z_top-z_bottom)`
  (taking the larger `hypot(patch_direction(...)...)` magnitude found
  across every edge's first sample point at either reference `z`, so a
  fluke-small value on one side doesn't spuriously weaken the reference).
  A candidate is gradient-suspect if the MINIMUM
  `hypot(patch_direction(...)...)` across all its own edges' anchors
  falls below `cfg.z_mid_gradient_ratio_tol` times this reference. This
  ratio is genuinely dimensionless (numerator and denominator are the
  same quantity at different `z`), so it self-scales with whatever
  gradient magnitude is "normal" for a given surface rather than needing
  a surface-specific absolute cutoff. Measured gap on the Taubin heart:
  the two already-healthy reference slabs scored `0.82`/`0.98` (i.e.
  nearly as strong as their own slab's quarter-point baseline); the
  `[-1,1]` slab's bad candidates (`z=0.02,-0.02,0.04,-0.04`) scored
  `0.0014`-`0.0057`; its eventually-accepted good candidate (`z=0.06`,
  independently confirmed via `track_face` to give max `|f|=2.4e-6`)
  scored `0.013`. The chosen default, `cfg.z_mid_gradient_ratio_tol =
  0.01`, sits with almost two orders of magnitude of margin on both
  sides of this specific gap -- not a threshold squeezed uncomfortably
  close to either boundary, unlike the rejected residual-magnitude
  approach above.

Reference-scale computation calls [`slice_at_z`](@ref) at the two
quarter-points, so it costs two extra 2D-curve decompositions per
non-empty slab -- memoized once per slab, and paid EAGERLY on the very
first (naive-midpoint) attempt, not only on slabs that go on to retry.
This eagerness is deliberate and load-bearing, not an oversight: the
gradient gate exists precisely to catch candidates that are
topologically CLEAN but numerically degenerate, and the naive midpoint
itself can be such a candidate, because the degenerate z (z=0 for the
Taubin family) is NOT a critical z-value -- it sits invisibly inside a
slab, and nothing pins it to the slab's exact center. Concretely
(measured 2026-07, see dev/scratch/scratch_robust_slice_eagerness_check.jl):
this same Taubin heart with `bbox_z = (-0.96, 1.3)` -- an ordinary
asymmetric bounding-box crop, which drops the critical value z=-1 and
makes the bottom slab [-0.96, 1.0] with naive midpoint z=0.02 -- yields
a topology-clean naive slice (2 ordinary `Critical` vertices, no
`:endpoint_fallback`, no `Singular`) whose downstream sweep measures
max `|f| ≈ 1.6`, versus `2.4e-7` at the gradient-gate-chosen `z ≈ 0.059`.
A retry-armed lazy variant (gradient gate active only from the first
retry onward, restoring the cost model this paragraph originally
claimed) was evaluated and REJECTED for exactly that reason: it
silently accepts that slab. The measured eager cost is a ~3x multiplier
on a healthy slab's slicing time (sphere/ellipsoid: 3.0x), accepted as
the price of the guarantee; cheapening the reference itself (e.g.
Gauss-Newton-projecting the candidate's own edge anchors to the
quarter-point z's via `_project_to_slice` instead of running two full
decompositions there) would change the gate's empirically measured
thresholds above and is left for a future pass. If the reference scale
itself comes back exactly zero (only possible if the quarter-points are
themselves degenerate), the gradient gate is skipped for that slab rather
than dividing by zero -- documented as a known limitation of this
heuristic, not silently miscounted as "healthy".

# Historical note: the min-vs-max false positive, and its Phase 8 resolution
The original gate compared the candidate's MIN anchor gradient against the
reference's MAX, which false-fires on curves with multiple, legitimately
very-differently-conditioned branches (a thin inner sliver's pinch anchor vs.
a strong outer loop). On the fixed-axis `[1.0, 1.0648]` slab this cost one
harmless avoidable retry (candidate min `0.035` vs reference max `3.94`,
ratio `0.0089`, just under the `0.01` threshold; both the rejected and the
accepted `z_mid` sweep fine, max `|f| ~ 2.4e-6`/`3.0e-6`). The Phase 8
five-seed rotated-Taubin regression then showed the FATAL form of the same
artifact: narrow slabs between two close genuine critical values (seed 3:
`[-0.864, -0.858]`, seed 4: `[0.884, 0.900]` -- nowhere near the singular
ellipse, no Singular vertices at all) carry a just-born tiny branch whose
Critical anchor gradient (`4e-4` / `3.5e-3`) is structurally weak across the
ENTIRE slab, so every retry candidate failed identically and the ladder
exhausted -- a hard throw. The structural-heterogeneity skip (refinement 3
in the Phase 8 section below) resolves both: measured `ref_min/ref_max` is
`0.00306` on the `[1.0, 1.0648]` slab (gate skipped, false positive gone,
naive midpoint accepted, sweep max `|f| = 2.4e-6`) and `0.372` / `0.0723` on
the `[-1, 1]` / `[1.0648, 1.2367]` slabs (gate active, all true positives
preserved) -- the `0.01` threshold sits ~3x from either side of this gap.

# Phase 8 amendments: transversal singular curves (generic-projection charts)
Both gates above were originally calibrated on fixed-axis fixtures, where a
positive-dimensional singular curve of the surface (e.g. the Taubin heart's
ellipse `{x^2+1.44y^2=1, z=0}`, along which `∇f` vanishes identically) is
either parallel to the slicing planes -- confined to a single degenerate z,
exactly the case the retry ladder handles -- or absent. Under a generic
rotated chart (`decompose_3d_surface`'s Phase 8 `projection` support), that
curve becomes TRANSVERSAL: it spans a whole range of chart-z values, and every
midslice in that range legitimately contains isolated singular points (nodes).
Measured on the rotated Taubin heart (seed 1, 2026-07): the naive midpoint of
slab `[-0.956, -0.24]` has 2 `Singular` vertices (the ellipse crossings), 2
benign branch-stitching `:endpoint_fallback` vertices, and per-edge
first-sample anchor gradients `{4.01, 4.01, 1.44, 8.6e-15, 1.8e-13, 0.28}`
against reference scale `3.98` -- BOTH gates false-fired on every retry
candidate (the same structure recurs at every nearby z), so 5 of 9 slabs
threw and the decomposition died. Two refinements fix this:

1. **Topology gate**: the co-occurrence only fires when the quarter-point
   reference slices do NOT also contain `Singular` vertices. A structural
   singularity cannot be retried away; the original z=0 catastrophe still
   fires because its references at z=±0.5 are Singular-free.
2. **Gradient gate, exclusion**: edges left-anchored AT a `Singular` vertex
   are excluded from gradient measurements, candidate and reference alike
   (their anchor gradient is ~0 because the anchor IS the node); if every
   candidate edge is excluded, the slab is suspect. A mid-edge anchor was
   evaluated as an alternative and REJECTED with data: `sample_edge` chords
   sit measurably off-curve precisely in the degenerate cases (the z=0.02
   catastrophe candidate reads a healthy `0.353` at its chord midpoint while
   its downstream sweep is provably catastrophic), which would blind the gate
   exactly where it is most needed.
3. **Gradient gate, structural-heterogeneity skip** (added after the 5-seed
   regression exposed the fatal form of the documented min-vs-max false
   positive -- see the Historical note above): if the slab's OWN reference
   slices fail the gradient criterion (`ref_min < z_mid_gradient_ratio_tol *
   ref_max` over eligible anchors), then a legitimately weak-gradient branch
   coexists with strong ones across the whole slab and the comparison cannot
   discriminate candidates -- the gate is skipped instead of failing every
   retry identically. This is the same principle as refinement 1: a signature
   the references share with every candidate carries no information about
   WHICH z_mid to prefer.

Validation of the refined gates (measured 2026-07): fixed-axis Taubin --
`[-1,1] -> 0.06` retried exactly as before (reference heterogeneity `0.372`,
gate active, the z=0 catastrophe ladder unchanged), `[1.0, 1.0648] ->` naive
midpoint (heterogeneity `0.00306`, gate skipped, resolving the documented
false positive; sweep max `|f| = 2.4e-6`), `[1.0648, 1.2367] ->` naive
midpoint (heterogeneity `0.0723`, gate active, no fire); full fixed-axis
decompose max `|f| = 2.4e-6`. Rotated Taubin, seeds 1-5: zero retries, zero
throws; seed 3 (previously fatal) full decompose median world-`|f|` `2.1e-8`
with 21/5677 points above `1e-4`, confined to the singular-curve band. Point
singularities also scatter their chart-z critical-value estimates (~2e-4
observed, beyond `vertex_match_tol`); the resulting sub-resolution sliver
slabs are merged upstream by `cfg.min_slab_width` in `_slab_bounds`, not
handled here.
"""
function _robust_slice_at_z(F::System, patch::NamedTuple, z_bottom::T, z_top::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    width = z_top - z_bottom
    max_multiple = max(1, floor(Int, T(0.45) / cfg.z_mid_retry_frac))

    _cooccurrence(vertices) = any(
        v -> v.v_type == Artificial && get(v.metadata, :origin, nothing) == :endpoint_fallback,
        vertices,
    ) && any(v -> v.v_type == Singular, vertices)

    _anchor_gradient_magnitude(edge::Edge{T}) = begin
        p = edge.sampled_points[1]
        a, b = patch_direction(patch, p[1], p[2], p[3], cfg)
        hypot(a, b)
    end

    # Phase 8 refinement: edges whose LEFT vertex (the anchor, sampled_points[1]) is
    # Singular-classified are excluded from gradient measurements — their anchor gradient
    # is ~0 because the anchor IS the singular point, not because the whole slice curve is
    # degenerate. Applied identically to candidate and reference slices.
    _eligible_edges(vertices, edges) = begin
        sing_ids = Set(v.id for v in vertices if v.v_type == Singular)
        [e for e in edges if !isempty(e.sampled_points) && !(e.left_vertex_id in sing_ids)]
    end

    # Memoized once per slab, computed EAGERLY on the naive-midpoint attempt (see the
    # eagerness section of the docstring). Phase 8 extends it to record (a) whether the
    # reference slices contain Singular vertices (the structural-singularity signal for
    # the refined topology gate) and (b) the min AND max eligible anchor gradients (the
    # structural-heterogeneity signal for the gradient-gate skip below).
    ref_info = Ref{Union{Nothing,Tuple{T,T,Bool}}}(nothing)
    function _reference_info()
        ref_info[] === nothing || return ref_info[]
        ref_max = zero(T)
        ref_min = T(Inf)
        ref_has_singular = false
        for z_ref in (z_bottom + T(0.25) * width, z_top - T(0.25) * width)
            v_ref, edges_ref = slice_at_z(F, z_ref, cfg)
            ref_has_singular |= any(v -> v.v_type == Singular, v_ref)
            for e in _eligible_edges(v_ref, edges_ref)
                g = _anchor_gradient_magnitude(e)
                ref_max = max(ref_max, g)
                ref_min = min(ref_min, g)
            end
        end
        ref_info[] = (ref_max, ref_min, ref_has_singular)
        return ref_info[]
    end

    # Phase 8 refinement: the fallback+Singular co-occurrence only counts against a
    # candidate when the slab's reference slices do NOT also contain Singular vertices.
    # If they do, a singular curve crosses the slab transversally and no z_mid choice
    # avoids it — retrying is definitionally useless (see the docstring's Phase 8 section).
    _topology_suspect(vertices) = _cooccurrence(vertices) && !_reference_info()[3]

    _gradient_suspect(vertices, edges) = begin
        isempty(edges) && return false
        ref_max, ref_min, _ = _reference_info()
        ref_max == zero(T) && return false
        # Phase 8 refinement (structural heterogeneity): if the slab's own reference
        # slices fail the gradient criterion themselves (a legitimately weak-gradient
        # branch coexists with strong ones across the WHOLE slab), the min-candidate-vs-
        # max-reference comparison cannot discriminate candidates here — skip the gate
        # rather than fail every retry identically (see the docstring's Phase 8 section).
        ref_min < cfg.z_mid_gradient_ratio_tol * ref_max && return false
        eligible = _eligible_edges(vertices, edges)
        isempty(eligible) && return true
        cand_scale = minimum(_anchor_gradient_magnitude(e) for e in eligible)
        cand_scale < cfg.z_mid_gradient_ratio_tol * ref_max
    end

    _suspect(vertices, edges) = _topology_suspect(vertices) || _gradient_suspect(vertices, edges)

    z_mid = (z_bottom + z_top) / T(2)
    vertices_3d, edges_3d = slice_at_z(F, z_mid, cfg)
    _suspect(vertices_3d, edges_3d) || return vertices_3d, edges_3d, z_mid

    for k in 1:cfg.max_z_mid_retries
        multiple = T(min(cld(k, 2), max_multiple))
        sign = isodd(k) ? T(1) : T(-1)
        z_try = (z_bottom + z_top) / T(2) + sign * multiple * cfg.z_mid_retry_frac * width
        vertices_3d, edges_3d = slice_at_z(F, z_try, cfg)
        if !_suspect(vertices_3d, edges_3d)
            return vertices_3d, edges_3d, z_try
        end
    end

    throw(ErrorException(
        "_robust_slice_at_z: slab [$z_bottom, $z_top] still fails the vertex-type and/or " *
        "gradient-magnitude gate after the original midpoint plus $(cfg.max_z_mid_retries) " *
        "perturbed retries (cfg.z_mid_retry_frac=$(cfg.z_mid_retry_frac), " *
        "cfg.z_mid_gradient_ratio_tol=$(cfg.z_mid_gradient_ratio_tol)); giving up rather than " *
        "silently proceeding with a slice known to be untrustworthy.",
    ))
end

"""
    _slab_bounds(F::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat} -> Vector{T}

Build the sorted z-slab boundary list for [`decompose_3d_surface`](@ref):
`cfg.bbox_z` endpoints plus in-range critical z-values, with boundaries closer
than `cfg.min_slab_width` merged.

The merge (Phase 8) exists because path endpoints landing on a point
singularity carry ~accuracy^(1/multiplicity) scatter in their z-estimate
(~2e-4 observed on the rotated Taubin heart, beyond `vertex_match_tol`'s
reach), which would otherwise mint sliver slabs centered ON a singular point
that no `z_mid` choice can slice. `cluster_scalars` averages each cluster;
the outermost bounds are clamped back to the exact bbox afterward. All
existing fixed-axis fixtures have critical-z gaps >= 0.065 (65x the default
floor), so their bounds are unchanged. This is a documented resolution limit:
genuinely distinct critical values closer than `cfg.min_slab_width` are not
separated.
"""
function _slab_bounds(F::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    z_bottom_bound, z_top_bound = cfg.bbox_z
    z_crits_raw = compute_critical_z_slices(F, cfg)
    # Critical z outside bbox_z is not a slab boundary for this decomposition.
    z_crits = filter(z -> z_bottom_bound <= z <= z_top_bound, z_crits_raw)
    z_bounds = sort(unique(vcat(T[z_bottom_bound], z_crits, T[z_top_bound])))
    length(z_bounds) > 2 || return z_bounds
    merged = cluster_scalars(z_bounds, cfg.min_slab_width)
    length(merged) < 2 && return T[z_bottom_bound, z_top_bound]
    merged[1] = z_bottom_bound
    merged[end] = z_top_bound
    return merged
end

"""
    weld_mesh(faces::Vector{Face{T}}, patch::NamedTuple, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> GeometryBasics.Mesh

Merge per-face meshes into one watertight `GeometryBasics.Mesh`.

Collects every face's `mesh_vertices` and `mesh_topology`, clusters coincident
vertices across face boundaries with `cfg.vertex_match_tol`, remaps triangle
indices, drops degenerate triangles, and flips winding so normals align with
`∇f`.

# Arguments
- `faces::Vector{Face{T}}`: swept faces with local mesh data.
- `patch::NamedTuple`: surface patch system (for gradient-based winding).
- `cfg::HomotopyConfig{T}`: vertex-matching tolerance.

# Returns
A welded `GeometryBasics.Mesh` with `Point3f` vertices.
"""
function weld_mesh(faces::Vector{Face{T}}, patch::NamedTuple, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    # Adjacent faces share boundary geometry (position-based merge, not id offset).
    all_points = Vector{T}[]
    provenance = Tuple{Int,Int}[]
    for (fi, face) in enumerate(faces)
        for r in 1:size(face.mesh_vertices, 1)
            push!(all_points, face.mesh_vertices[r, :])
            push!(provenance, (fi, r))
        end
    end

    reps, membership = cluster_points_indexed(all_points, cfg.vertex_match_tol)

    lookup = Dict{Tuple{Int,Int},Int}()
    for (k, fr) in enumerate(provenance)
        lookup[fr] = membership[k]
    end

    global_triangles = NTuple{3,Int}[]
    for (fi, face) in enumerate(faces)
        for row in 1:size(face.mesh_topology, 1)
            i1, i2, i3 = face.mesh_topology[row, 1], face.mesh_topology[row, 2], face.mesh_topology[row, 3]
            g1, g2, g3 = lookup[(fi, i1)], lookup[(fi, i2)], lookup[(fi, i3)]
            # Drop pinched triangles (e.g. pole/tip where three indices collapse).
            g1 != g2 && g2 != g3 && g3 != g1 && push!(global_triangles, (g1, g2, g3))
        end
    end

    # Global winding fix: track_face emits no orientation guarantee; align each normal with +∇f.
    fixed_triangles = NTuple{3,Int}[]
    for (g1, g2, g3) in global_triangles
        p1, p2, p3 = reps[g1], reps[g2], reps[g3]
        n = cross(p2 .- p1, p3 .- p1)
        gx, gy, gz = _gradient_at(patch, p1[1], p1[2], p1[3], cfg)
        push!(fixed_triangles, dot(n, T[gx, gy, gz]) < 0 ? (g1, g3, g2) : (g1, g2, g3))
    end

    # Cluster at T-precision; narrow to Float32 only for the GeometryBasics.Mesh container.
    points3 = [GeometryBasics.Point3f(Float32(p[1]), Float32(p[2]), Float32(p[3])) for p in reps]
    tris = [GeometryBasics.TriangleFace{Int}(t[1], t[2], t[3]) for t in fixed_triangles]
    return GeometryBasics.Mesh(points3, tris)
end

"""
    decompose_3d_surface(F::System, cfg::HomotopyConfig{T}; projection = nothing, rng = Random.default_rng()) where {T<:AbstractFloat}
        -> (vertices, edges, faces, mesh)

Decompose a real algebraic surface into vertices, edges, faces, and a welded mesh.

`F` must be a single equation in three variables ordered `[x, y, z]` (z last).
The pipeline finds critical z-slices inside `cfg.bbox_z`, decomposes a mid-z
curve in each slab, sweeps faces between slab bounds, and welds everything into
one `GeometryBasics.Mesh`. This is the 3D analogue of [`decompose_1d_curve`](@ref).

# Arguments
- `F::System`: surface system (`length(F.expressions) == 1`, `nvariables(F) == 3`).
- `cfg::HomotopyConfig{T}`: tolerances, bounding box, and sampling densities.

# Keyword arguments (Phase 8)
- `projection`: `nothing` (default) sweeps along the literal z axis with the
  exact pre-Phase-8 behavior -- no coordinate transform is applied at all.
  `:random` draws a Haar-uniform SO(3) rotation internally (seed via `rng`);
  a user-supplied 3x3 orthonormal, det = +1 matrix `Q` is validated and used
  directly. When a projection is active the pipeline runs on the chart system
  `F'(x') = F(Q x')` over the enclosing axis-aligned box of the rotated world
  bbox (slightly larger than the world box at the same `cfg`);
  [`_verify_projection_ok`](@ref) rejects projections that leave an augmenting
  partial identically zero (the `z - x^2` crash class) with a loud
  `ArgumentError`; and every returned vertex coordinate, edge sample, face
  mesh row, and welded mesh point is mapped back to world coordinates via
  `p = Q p'`. `Face.mid_slice_z` and boundary-vertex
  `:fixed_variable`/`:fixed_value` metadata remain CHART-frame quantities.
  Callers who need `Q` itself should construct it with
  [`random_orthogonal_matrix`](@ref) and pass the matrix explicitly --
  `:random` is convenience for exploratory use.
- `rng`: random source used only by `projection = :random`.

# Known limitation: generic projections over singular curves
When a projection makes a positive-dimensional singular curve of the surface
(`∇f = 0` along a curve, e.g. the Taubin heart's z=0 ellipse) transversal to
the slicing planes, mesh quality degrades in a narrow band around that curve.
Measured (rotated Taubin heart, seed 1, 2026-07): median world-`|f|` residual
`6e-9` over 2617 mesh points, but ~2.5% of points exceed `1e-4` (max `0.26`),
ALL confined to `|z_world| <= 0.14` around the singular plane (the surface
spans `|z| <= 1.24`). Away from the singular locus the decomposition is
unaffected. Decomposing singular curves themselves (BertiniReal's
deflation/singular-curve machinery) is future work; see also
`_robust_slice_at_z`'s Phase 8 amendments section.

# Returns
A 4-tuple `(vertices, edges, faces, mesh)` of types
`Vector{NativeVertex{T}}`, `Vector{Edge{T}}`, `Vector{Face{T}}`, and
`GeometryBasics.Mesh`, in world coordinates.
"""
function decompose_3d_surface(
    F::System,
    cfg::HomotopyConfig{T};
    projection::Union{Nothing,Symbol,AbstractMatrix} = nothing,
    rng::Random.AbstractRNG = Random.default_rng(),
) where {T<:AbstractFloat}
    length(F.variables) == 3 && length(F.expressions) == 1 || throw(ArgumentError(
        "decompose_3d_surface: expected F with exactly 1 equation in exactly 3 variables " *
        "(a raw surface, variables ordered [x_var, y_var, z_var] -- z LAST); " *
        "got $(length(F.expressions)) equation(s) in $(length(F.variables)) variable(s).",
    ))

    # Phase 8: change-of-coordinates wrapper. The chart run recurses with
    # projection = nothing, so the entire validated fixed-axis pipeline below
    # is reused unchanged; transforms live only at this boundary.
    if projection !== nothing
        Q = _resolve_projection(projection, rng)
        F_chart = _rotate_system(F, Q)
        _verify_projection_ok(F_chart, cfg)
        cfg_chart = _chart_config(cfg, Q)
        chart_result = decompose_3d_surface(F_chart, cfg_chart)
        return _map_to_world(chart_result..., Q)
    end

    # bbox_z endpoints + in-range critical z-values, with sub-resolution bounds
    # merged per cfg.min_slab_width (see _slab_bounds).
    z_bounds = _slab_bounds(F, cfg)

    # One patch for the whole surface: reused by _robust_slice_at_z (gradient
    # gate via FaceTracking.patch_direction) and by FaceTracking.track_face.
    patch = build_patch_system(F)

    all_vertices = NativeVertex{T}[]
    all_edges = Edge{T}[]
    all_faces = Face{T}[]
    next_face_id = 1

    for i in 1:(length(z_bounds) - 1)
        z_bottom, z_top = z_bounds[i], z_bounds[i+1]
        # Exact midpoint (same convention as Topology.compute_midslice). The old
        # prototype's 0.4137 skew (prototipo_viejo_julia/SurfaceTopology.jl:132)
        # is deliberately dropped, not ported to HomotopyConfig. _robust_slice_at_z
        # may return a nearby z_mid if the midpoint slice fails its gates; track_face
        # must use that same z_mid, not re-derive (z_bottom+z_top)/2.
        vertices_2d, edges_2d, z_mid = _robust_slice_at_z(F, patch, z_bottom, z_top, cfg)

        # Each slice_at_z / decompose_1d_curve restarts ids near 1. Concatenating
        # slabs without offsets would collide namespaces (not geometric merges —
        # each slab has a distinct z_mid). Offset vertex and edge ids independently
        # (same separation decompose_1d_curve keeps within one call); shift
        # left/right_vertex_id by the vertex offset. Analogue of Topology.jl's
        # `offset = maximum(v.id for v in crit_vertices)` pattern, across slabs.
        v_offset = isempty(all_vertices) ? 0 : maximum(v.id for v in all_vertices)
        e_offset = isempty(all_edges) ? 0 : maximum(e.id for e in all_edges)

        vertices_renumbered = NativeVertex{T}[
            NativeVertex{T}(id = v.id + v_offset, coordinates = v.coordinates, v_type = v.v_type, metadata = v.metadata)
            for v in vertices_2d
        ]
        edges_renumbered = Edge{T}[
            Edge{T}(
                id = e.id + e_offset,
                left_vertex_id = e.left_vertex_id + v_offset,
                right_vertex_id = e.right_vertex_id + v_offset,
                sampled_points = e.sampled_points,
                is_singular = e.is_singular,
            )
            for e in edges_2d
        ]

        append!(all_vertices, vertices_renumbered)
        append!(all_edges, edges_renumbered)

        for edge in edges_renumbered
            # Face.boundary_edges come from already-renumbered edges.
            face = track_face(F, patch, edge, z_mid, z_bottom, z_top, next_face_id, cfg)
            push!(all_faces, face)
            next_face_id += 1
        end
    end

    mesh = weld_mesh(all_faces, patch, cfg)
    return all_vertices, all_edges, all_faces, mesh
end
