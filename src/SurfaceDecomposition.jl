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
    CritSlice{T<:AbstractFloat}

Phase 9a: the 1D decomposition of the surface AT one interior slab boundary
(a critical z-value), as opposed to the mid-slab slices the pipeline sweeps
from. Computed by [`_decompose_crit_slice`](@ref) when
[`decompose_3d_surface`](@ref) runs with `incidence = true`.

# Fields
- `boundary_index`: `j` such that this slice sits at `z_bounds[j]` (interior
  boundaries only; bbox endpoints never get a `CritSlice`).
- `z`: the critical z-value (CHART frame under a `projection`, like
  `Face.mid_slice_z`).
- `vertices`, `edges`: the slice's cells, id-renumbered into the same global
  namespace as the main decomposition's cells but NOT appended to the main
  return (backward compatibility; ids never collide).
- `is_degenerate`: `true` iff the slice has no edges — a fold/point-type
  boundary (measured: sphere/ellipsoid extremes decompose to EMPTY slices,
  Taubin cusp/tips to isolated `Singular` vertices; nodal saddle-type
  boundaries decompose to real curves). Fold-type boundaries are anchored to
  the surface's 3D critical vertices instead (see [`_assign_landings`](@ref)).
"""
struct CritSlice{T<:AbstractFloat}
    boundary_index::Int
    z::T
    vertices::Vector{NativeVertex{T}}
    edges::Vector{Edge{T}}
    is_degenerate::Bool
end

"""
    ColumnLanding{T<:AbstractFloat}

Phase 9a: the incidence assignment of ONE swept column's landing point at one
slab boundary. `kind` is `:edge` (a `CritSlice` edge), `:crit_slice_vertex`,
`:critical_point` (a 3D surface critical vertex — the fold-type anchor),
`:bbox` (the boundary is a bounding-box endpoint, no crit-slice exists), or
`:none` (no candidate cell found); `id` is the cell id in the matching
namespace (0 for `:bbox`/`:none`); `dist` is the Euclidean distance from the
landing point to the assigned cell.

Distances are RECORDED, never thresholded, in Phase 9a: assignment is
nearest-wins, and the printed distance distributions (see
test/test_taubin.jl's 9a section) are the evidence base from which Phase 9b
will derive real acceptance thresholds — deliberately not guessed ahead of
the data (the same discipline the Phase 8 gate history taught).
"""
struct ColumnLanding{T<:AbstractFloat}
    kind::Symbol
    id::Int
    dist::T
end

"""
    SurfaceIncidence{T<:AbstractFloat}

Phase 9a: face/edge incidence data for a surface decomposition — which
crit-slice cells each face's swept columns land on at its slab boundaries.
Returned as the 5th element of [`decompose_3d_surface`](@ref) when
`incidence = true`.

# Fields
- `crit_slices`: one [`CritSlice`](@ref) per interior slab boundary that an
  actual face touches (lazily computed; sorted by `boundary_index`).
- `face_bottom_edges` / `face_top_edges`: face id → ordered unique crit-slice
  edge ids its boundary columns land on.
- `face_bottom_anchor` / `face_top_anchor`: face id → 3D critical-vertex id
  for fold-type boundaries (0 when none).
- `edge_faces`: inverse adjacency, crit-slice edge id → face ids touching it.
- `column_landings_bottom` / `column_landings_top`: face id → per-column
  [`ColumnLanding`](@ref) in curve order.
- `critical_vertices` (Phase 9b): the surface's 3D critical points (fold-type
  landing anchors), stored so `:critical_point`-kind landing ids resolve to
  coordinates — needed for snap-unification, not stored in Phase 9a since
  nothing there looked them up by id after assignment.
- `continuity_ok` (Phase 9b): face id → `true` iff every pair of CONFIDENT
  (see `_landing_confidence`) consecutive columns, on both sides, lands on
  combinatorially adjacent crit-slice cells (see `_cells_adjacent`). Absent
  keys mean the face touched no crit-slice (both sides `:bbox`/empty).
- `continuity_violations_bottom` / `continuity_violations_top` (Phase 9b):
  face id → indices `c` of violating consecutive pairs `(c, c+1)`.

Promoting `continuity_ok` from a flag to a hard error is a deliberate FUTURE
decision contingent on the firing-rate evidence this structure gathers on
real fixtures — not something decided now. Measured on the Taubin heart
(2026-07): fixed-axis 10/14 faces `continuity_ok` (16 flagged pairs),
rotated seed 1 15/22 faces (17 flagged pairs); sphere/ellipsoid (single
fold-type boundary each) show zero. The flagged pairs concentrate at the
SAME multi-face edge-type boundaries `weld_mesh`'s Phase 9b docstring
section documents as an open coverage-gap issue — this is exactly the
evidence-gathering this field exists for, not a surprise finding.
"""
struct SurfaceIncidence{T<:AbstractFloat}
    crit_slices::Vector{CritSlice{T}}
    critical_vertices::Vector{NativeVertex{T}}
    face_bottom_edges::Dict{Int,Vector{Int}}
    face_top_edges::Dict{Int,Vector{Int}}
    face_bottom_anchor::Dict{Int,Int}
    face_top_anchor::Dict{Int,Int}
    edge_faces::Dict{Int,Vector{Int}}
    column_landings_bottom::Dict{Int,Vector{ColumnLanding{T}}}
    column_landings_top::Dict{Int,Vector{ColumnLanding{T}}}
    continuity_ok::Dict{Int,Bool}
    continuity_violations_bottom::Dict{Int,Vector{Int}}
    continuity_violations_top::Dict{Int,Vector{Int}}
end

"""
    _decompose_crit_slice(F, z_c::T, j::Int, cfg, v_offset::Int, e_offset::Int) -> CritSlice{T}

Decompose the surface AT the exact interior boundary value `z_c = z_bounds[j]`
via [`slice_at_z`](@ref), deliberately bypassing `_robust_slice_at_z`: the
value is pinned by definition (it IS the critical value), so perturbed
retries are meaningless here. Cell ids are shifted by the offsets into the
global namespace. Measured on every fixture (fixed-axis and rotated Taubin
included): nodal/saddle-type crit-slices decompose like the nodal-cubic
fixture, fold-type extremes come back empty or as isolated `Singular`
vertices, and nothing throws.
"""
function _decompose_crit_slice(F::System, z_c::T, j::Int, cfg::HomotopyConfig{T}, v_offset::Int, e_offset::Int) where {T<:AbstractFloat}
    vertices_raw, edges_raw = slice_at_z(F, z_c, cfg)
    vertices = NativeVertex{T}[
        NativeVertex{T}(id = v.id + v_offset, coordinates = v.coordinates, v_type = v.v_type, metadata = v.metadata)
        for v in vertices_raw
    ]
    edges = Edge{T}[
        Edge{T}(
            id = e.id + e_offset,
            left_vertex_id = e.left_vertex_id + v_offset,
            right_vertex_id = e.right_vertex_id + v_offset,
            sampled_points = e.sampled_points,
            is_singular = e.is_singular,
        )
        for e in edges_raw
    ]
    return CritSlice{T}(j, z_c, vertices, edges, isempty(edges))
end

"""
    _surface_critical_vertices(F, cfg) -> Vector{NativeVertex{T}}

The surface's 3D critical points (the same `compute_critical_points` call
`compute_critical_z_slices` makes internally), needed as fold-type landing
anchors. Called ONCE per `incidence = true` decomposition — an extra
~seconds-scale solve that keeps the `incidence = false` call chain
byte-identical; deduplicating it with `_slab_bounds`'s internal call is a
deliberately deferred optimization.
"""
_surface_critical_vertices(F::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat} =
    compute_critical_points(F, cfg)

"""
    _assign_landings(landings, cs, crit_vertices, z_c, cfg) -> Vector{ColumnLanding{T}}

Nearest-wins incidence assignment for one face side. `landings` are the
per-column `(x, y, z)` landing points (rows 1 / n_z of `Face.mesh_vertices`);
`cs === nothing` means the boundary is a bbox endpoint (`:bbox` for every
column). Otherwise each landing competes across: the crit-slice's edge
polylines (minimum distance over `sampled_points` — chord-interpolated, hence
the recorded distances carry the known chord error), the crit-slice's
vertices, and the surface's 3D critical vertices within `cfg.min_slab_width`
of `z_c` (the fold-type anchors, which the slice itself represents
unreliably). Ties break by that category order. No candidates at all gives
`:none` with `dist = Inf` — never an error in 9a (pure measurement).
"""
function _assign_landings(
    landings::Vector{Vector{T}},
    cs::Union{CritSlice{T},Nothing},
    crit_vertices::Vector{NativeVertex{T}},
    z_c::T,
    cfg::HomotopyConfig{T},
) where {T<:AbstractFloat}
    cs === nothing && return ColumnLanding{T}[ColumnLanding{T}(:bbox, 0, zero(T)) for _ in landings]

    anchors = [v for v in crit_vertices if abs(real(v.coordinates[3]) - z_c) <= cfg.min_slab_width]

    out = Vector{ColumnLanding{T}}(undef, length(landings))
    for (c, p) in enumerate(landings)
        best_kind = :none
        best_id = 0
        best_dist = T(Inf)
        for e in cs.edges
            isempty(e.sampled_points) && continue
            d = minimum(norm(p .- q) for q in e.sampled_points)
            if d < best_dist
                best_kind, best_id, best_dist = :edge, e.id, d
            end
        end
        for v in cs.vertices
            d = norm(p .- real.(v.coordinates))
            if d < best_dist
                best_kind, best_id, best_dist = :crit_slice_vertex, v.id, d
            end
        end
        for v in anchors
            d = norm(p .- real.(v.coordinates))
            if d < best_dist
                best_kind, best_id, best_dist = :critical_point, v.id, d
            end
        end
        out[c] = ColumnLanding{T}(best_kind, best_id, best_dist)
    end
    return out
end

"""
    _median_spacing(pts::Vector{Vector{T}}) where {T<:AbstractFloat} -> T

Median consecutive-point spacing along an ordered point sequence (a
crit-slice edge's `sampled_points`, or one face side's landing points) --
the local chord/column scale Phase 9b's confidence gate measures distances
against. `< 2` points gives `zero(T)` (no local scale defined).
"""
function _median_spacing(pts::Vector{Vector{T}}) where {T<:AbstractFloat}
    length(pts) < 2 && return zero(T)
    seglens = sort(T[norm(pts[k+1] .- pts[k]) for k in 1:(length(pts)-1)])
    return seglens[cld(length(seglens), 2)]
end

_edge_spacing(edge::Edge{T}) where {T<:AbstractFloat} = _median_spacing(edge.sampled_points)

"""
    _inward_row_points(face::Face{T}, side::Symbol, n_z::Int, n_curve::Int) where {T<:AbstractFloat} -> Vector{Vector{T}}

The row ONE STEP IN from `side`'s boundary row of `face`'s swept grid
(`side = :bottom` -> row 2, `side = :top` -> row `n_z - 1`), read directly
from `face.mesh_vertices` (interior rows are never mutated by snapping).
Used as the local-scale reference for `:crit_slice_vertex`/`:critical_point`
confidence (see [`_landing_confidence`](@ref)): the BOUNDARY row's own
inter-column spacing collapses toward zero at a converging fold point (every
column approaching the SAME limit), which self-referentially defeats a ratio
check against the very distance it is supposed to bound. Measured on the
fixed-axis Taubin cusp (2026-07): boundary-row spacing `0.00355` vs a
landing distance of `0.0037` (fails a `1.5x` ratio against itself) while the
row-inward spacing is `0.0982` (comfortably passes) -- the SAME columns,
one ring earlier, before they have collapsed together.
"""
function _inward_row_points(face::Face{T}, side::Symbol, n_z::Int, n_curve::Int) where {T<:AbstractFloat}
    r = side === :bottom ? 2 : n_z - 1
    return Vector{T}[Vector{T}(face.mesh_vertices[(r-1)*n_curve+c, :]) for c in 1:n_curve]
end

"""
    _landing_confidence(point, landing, edge_spacing, scale_ref, patch, cfg) where {T<:AbstractFloat} -> Bool

Phase 9b: whether `landing` is confident enough to participate in EITHER
snap-unification ([`_snap_boundary_points!`](@ref)) or the continuity check
([`_check_continuity!`](@ref)) -- the SAME gate for both, per the Phase 9b
design. An `:edge` landing is confident within
`cfg.incidence_snap_tol_ratio` times its assigned edge's local chord spacing
(`edge_spacing`) AND within `cfg.critical_point_tol` surface residual --
calibrated against the Phase 8 singular-band evidence, where off-surface
`:edge` landings measure residual `~0.257` (six orders above the default
`critical_point_tol`) while on-surface ones measure `~1e-17`, a clean,
well-separated gap.

A `:crit_slice_vertex`/`:critical_point` landing is confident within the
SAME ratio of `scale_ref` (the [`_inward_row_points`](@ref) spacing -- a
single point has no chord spacing of its own) but is NOT additionally
residual-gated: measured on the fixed-axis Taubin cusp, genuinely converging
columns can carry surface residual `~2e-6` (just over the default
`critical_point_tol = 1e-6`) purely from Newton's convergence rate slowing
near a multiple root (the same phenomenon `_newton_polish`'s own docstring
notes elsewhere in this codebase) -- not an off-surface signal the way it is
for `:edge` landings, since a fold/point cell has no "incomplete branch
representation" failure mode to guard against. Snapping such a point only
IMPROVES its on-surface quality (it moves TO the cell's own, independently
validated, near-machine-precision coordinates), so omitting the residual
gate here carries no correctness risk once the distance ratio itself is
sound. `:bbox`/`:none` landings are never confident.
"""
function _landing_confidence(
    point::AbstractVector{T},
    landing::ColumnLanding{T},
    edge_spacing::Dict{Int,T},
    scale_ref::T,
    patch::NamedTuple,
    cfg::HomotopyConfig{T},
) where {T<:AbstractFloat}
    landing.kind in (:edge, :crit_slice_vertex, :critical_point) || return false
    scale = landing.kind === :edge ? get(edge_spacing, landing.id, zero(T)) : scale_ref
    scale > zero(T) || return false
    landing.dist <= cfg.incidence_snap_tol_ratio * scale || return false
    landing.kind === :edge || return true
    resid = _residual_at(patch, T(point[1]), T(point[2]), T(point[3]), cfg)
    return abs(resid) <= cfg.critical_point_tol
end

"""
    _cells_adjacent(cs::CritSlice{T}, a::ColumnLanding{T}, b::ColumnLanding{T}) where {T<:AbstractFloat} -> Bool

Combinatorial adjacency of two landing assignments on the SAME crit-slice
`cs`: identical cells, two `:edge`s sharing an endpoint vertex id, or an
`:edge` and one of its own endpoint `:crit_slice_vertex`. Every other kind
combination (including any pairing involving `:critical_point`, or two
distinct `:crit_slice_vertex` ids) is NOT adjacent -- deliberately
conservative: at a degenerate (fold-type) boundary distinct vertex/anchor ids
represent genuinely different points, and flagging non-adjacency there is
safe under Phase 9b's warn-not-error policy. Used only by
[`_check_continuity!`](@ref); `:bbox`/`:none` landings never reach this call
(excluded by [`_landing_confidence`](@ref) first).
"""
function _cells_adjacent(cs::CritSlice{T}, a::ColumnLanding{T}, b::ColumnLanding{T}) where {T<:AbstractFloat}
    a.kind === b.kind && a.id == b.id && return true
    _endpoints(eid) = begin
        e = only(filter(x -> x.id == eid, cs.edges))
        (e.left_vertex_id, e.right_vertex_id)
    end
    if a.kind === :edge && b.kind === :edge
        la, ra = _endpoints(a.id)
        lb, rb = _endpoints(b.id)
        return la in (lb, rb) || ra in (lb, rb)
    elseif a.kind === :edge && b.kind === :crit_slice_vertex
        la, ra = _endpoints(a.id)
        return b.id == la || b.id == ra
    elseif a.kind === :crit_slice_vertex && b.kind === :edge
        return _cells_adjacent(cs, b, a)
    else
        return false
    end
end

"""
    _check_continuity!(face_id, points, landings, cs, edge_spacing, scale_ref_points, patch, cfg, ok_dict, viol_dict)

Phase 9b branch-continuity check for ONE face side: for every pair of
consecutive columns that are BOTH confident (see
[`_landing_confidence`](@ref)), their assigned cells must be
[`_cells_adjacent`](@ref); any confirmed non-adjacent pair sets
`ok_dict[face_id] = false` (ANDed with whatever the other side already
recorded) and appends the left column's index to `viol_dict[face_id]`. Never
throws -- Phase 9b's approved policy is warn/flag, not error; promotion to a
hard error is a deliberate future decision. `cs === nothing` (a bbox-type
boundary) or fewer than 2 columns trivially passes.

The confidence gate is what makes this check meaningful rather than noisy:
distances beyond `cfg.incidence_snap_tol_ratio * scale` are exactly the
ambiguous cases (chord error near coarse crit-slice sampling, or the Phase 8
singular-band seams) where a jump/no-jump call would be unfounded, so those
pairs are recorded as inconclusive instead of flagged. This is a documented
resolution limit: two branches closer together than roughly `2 *`
the local chord spacing are not distinguishable from a jump at a given
crit-slice sampling density.
"""
function _check_continuity!(
    face_id::Int,
    points::Vector{Vector{T}},
    landings::Vector{ColumnLanding{T}},
    cs::Union{CritSlice{T},Nothing},
    edge_spacing::Dict{Int,T},
    scale_ref_points::Vector{Vector{T}},
    patch::NamedTuple,
    cfg::HomotopyConfig{T},
    ok_dict::Dict{Int,Bool},
    viol_dict::Dict{Int,Vector{Int}},
) where {T<:AbstractFloat}
    viol_dict[face_id] = get(viol_dict, face_id, Int[])
    if cs === nothing || length(points) < 2
        ok_dict[face_id] = get(ok_dict, face_id, true)
        return nothing
    end
    scale_ref = _median_spacing(scale_ref_points)
    confident = Bool[_landing_confidence(points[c], landings[c], edge_spacing, scale_ref, patch, cfg) for c in eachindex(points)]
    ok = get(ok_dict, face_id, true)
    for c in 1:(length(points) - 1)
        (confident[c] && confident[c+1]) || continue
        if !_cells_adjacent(cs, landings[c], landings[c+1])
            ok = false
            push!(viol_dict[face_id], c)
        end
    end
    ok_dict[face_id] = ok
    return nothing
end

"""
    _map_incidence_to_world(inc::SurfaceIncidence{T}, Q::Matrix{Float64}) -> SurfaceIncidence{T}

World-maps a chart-frame [`SurfaceIncidence`](@ref) under a Phase 8
projection: each `CritSlice`'s vertices/edges and `critical_vertices` rotate
exactly like the main decomposition's cells in `_map_to_world`; ids,
adjacency dicts, landing assignments, and continuity data are frame-free
(`dist` is rotation-invariant, adjacency and continuity are purely
combinatorial) and carried over unchanged. Lives here rather than
Projection.jl because the incidence structs are defined in this file
(include order: Projection.jl precedes it).
"""
function _map_incidence_to_world(inc::SurfaceIncidence{T}, Q::Matrix{Float64}) where {T<:AbstractFloat}
    QT = Matrix{T}(Q)
    crit_slices = CritSlice{T}[
        CritSlice{T}(
            cs.boundary_index,
            cs.z,
            NativeVertex{T}[
                NativeVertex{T}(id = v.id, coordinates = QT * v.coordinates, v_type = v.v_type, metadata = v.metadata)
                for v in cs.vertices
            ],
            Edge{T}[
                Edge{T}(
                    id = e.id,
                    left_vertex_id = e.left_vertex_id,
                    right_vertex_id = e.right_vertex_id,
                    sampled_points = [QT * p for p in e.sampled_points],
                    is_singular = e.is_singular,
                )
                for e in cs.edges
            ],
            cs.is_degenerate,
        )
        for cs in inc.crit_slices
    ]
    critical_vertices = NativeVertex{T}[
        NativeVertex{T}(id = v.id, coordinates = QT * v.coordinates, v_type = v.v_type, metadata = v.metadata)
        for v in inc.critical_vertices
    ]
    return SurfaceIncidence{T}(
        crit_slices,
        critical_vertices,
        inc.face_bottom_edges,
        inc.face_top_edges,
        inc.face_bottom_anchor,
        inc.face_top_anchor,
        inc.edge_faces,
        inc.column_landings_bottom,
        inc.column_landings_top,
        inc.continuity_ok,
        inc.continuity_violations_bottom,
        inc.continuity_violations_top,
    )
end

"""
    _monotone_snap_targets(points::Vector{Vector{T}}, targets::Vector{Vector{T}}) where {T<:AbstractFloat} -> Vector{Int}

Phase 9b: assigns each of `points` (a contiguous, same-crit-slice-edge run of
one face side's swept columns, in curve order) to a target index into
`targets` (that edge's Gauss-Newton re-projected samples, also in curve
order), minimizing total assignment distance SUBJECT TO the assigned indices
being monotone (whichever of non-decreasing / non-increasing gives lower
total cost) in the SAME order as `points`.

Motivated directly by measurement: plain nearest-target assignment (no
monotonicity constraint) was checked on the fixed-axis Taubin `z=1.0`
boundary and visits targets non-monotonically on more than half of same-edge
column runs (4/8 sides monotone without this constraint, 16/56 segments
spanning >= 2 targets) -- which would zigzag the welded boundary. This DP
eliminates that by construction: `O(n*m)` two-pass (one per direction)
dynamic program with a running-prefix-minimum (`n = length(points)`,
`m = length(targets)`), trivial at realistic sizes.

Known scope limit: the DP searches the FULL target list, not a local window,
so a self-intersecting or sharply doubling-back crit-slice edge could in
principle attract a run to a distant, geometrically-unrelated section of the
SAME edge. Not observed on any fixture measured so far; flagged rather than
guarded against, per this project's "don't guess ahead of evidence" practice.
"""
function _monotone_snap_targets(points::Vector{Vector{T}}, targets::Vector{Vector{T}}) where {T<:AbstractFloat}
    n = length(points)
    m = length(targets)
    n == 0 && return Int[]

    function solve(increasing::Bool)
        dp = fill(T(Inf), n, m)
        prev = zeros(Int, n, m)
        for j in 1:m
            dp[1, j] = norm(points[1] .- targets[j])
        end
        for i in 2:n
            jr = increasing ? (1:m) : (m:-1:1)
            runmin = T(Inf)
            runarg = 0
            for j in jr
                if dp[i-1, j] < runmin
                    runmin = dp[i-1, j]
                    runarg = j
                end
                dp[i, j] = runmin + norm(points[i] .- targets[j])
                prev[i, j] = runarg
            end
        end
        jbest = argmin(view(dp, n, :))
        path = Vector{Int}(undef, n)
        path[n] = jbest
        for i in n:-1:2
            path[i-1] = prev[i, path[i]]
        end
        return path, dp[n, jbest]
    end

    path_inc, cost_inc = solve(true)
    path_dec, cost_dec = solve(false)
    return cost_inc <= cost_dec ? path_inc : path_dec
end

"""
    _snap_boundary_points!(all_points, index_of, faces, incidence, patch, cfg) where {T<:AbstractFloat} -> Dict{Tuple{Int,Int},Int}

Phase 9b: mutates `all_points` in place, snapping every CONFIDENT (see
[`_landing_confidence`](@ref)) boundary-column landing onto a shared target:
`:edge`-kind landings in a contiguous same-edge column run snap to that
edge's Gauss-Newton re-projected samples via
[`_monotone_snap_targets`](@ref); `:crit_slice_vertex`/`:critical_point`-kind
landings snap directly to that cell's own coordinates. Two faces whose
columns both land near the SAME target now hold EXACTLY equal coordinates, so
the subsequent `cluster_points_indexed` call merges them -- this is what
closes the tolerance-welding cracks `_naked_mesh_edges` measures.
Non-confident, `:bbox`, and `:none` landings are left untouched.

Returns `target_flat_idx`: for every `(crit-slice edge id, target index)`
pair actually used by some snapped column, ONE representative flat index into
`all_points` -- consumed by [`_split_t_junctions`](@ref) after clustering to
build each edge's used-target polyline.
"""
function _snap_boundary_points!(
    all_points::Vector{Vector{T}},
    index_of::Dict{Tuple{Int,Int},Int},
    faces::Vector{Face{T}},
    incidence::SurfaceIncidence{T},
    patch::NamedTuple,
    cfg::HomotopyConfig{T},
) where {T<:AbstractFloat}
    id_to_fi = Dict(f.id => i for (i, f) in enumerate(faces))

    edge_spacing = Dict{Int,T}(e.id => _edge_spacing(e) for cs in incidence.crit_slices for e in cs.edges)
    edge_targets = Dict{Int,Vector{Vector{T}}}()
    for cs in incidence.crit_slices, e in cs.edges
        tg = Vector{Vector{T}}(undef, length(e.sampled_points))
        for (k, p) in enumerate(e.sampled_points)
            xr, yr = _project_to_slice(patch, p[1], p[2], cs.z, cfg)
            tg[k] = T[xr, yr, cs.z]
        end
        edge_targets[e.id] = tg
    end
    cell_point_vertex = Dict{Int,Vector{T}}()
    for cs in incidence.crit_slices, v in cs.vertices
        cell_point_vertex[v.id] = T[real(v.coordinates[1]), real(v.coordinates[2]), cs.z]
    end
    cell_point_anchor = Dict{Int,Vector{T}}(
        v.id => T[real(v.coordinates[1]), real(v.coordinates[2]), real(v.coordinates[3])]
        for v in incidence.critical_vertices
    )

    target_flat_idx = Dict{Tuple{Int,Int},Int}()

    # Pass 1: snap :crit_slice_vertex/:critical_point immediately (every column
    # confidently assigned to the SAME cell gets the exact same coordinates
    # regardless of which face touches it -- no cross-face coordination needed).
    # :edge-kind confident columns are POOLED by crit-slice edge id, across ALL
    # faces, for a joint assignment in pass 2 (see its own comment for why).
    pooled_edge_cols = Dict{Int,Vector{Tuple{Int,Vector{T}}}}()
    for face in faces
        n_pts = size(face.mesh_vertices, 1)
        n_pts == 0 && continue
        n_z = 2 * cfg.midslice_sample_density + 1
        n_curve = n_pts ÷ n_z
        fi = id_to_fi[face.id]

        for (landings_dict, row_of, side) in (
            (incidence.column_landings_bottom, c -> c, :bottom),
            (incidence.column_landings_top, c -> (n_z - 1) * n_curve + c, :top),
        )
            haskey(landings_dict, face.id) || continue
            landings = landings_dict[face.id]
            length(landings) == n_curve || continue

            flat_idxs = [index_of[(fi, row_of(c))] for c in 1:n_curve]
            points = [all_points[idx] for idx in flat_idxs]
            scale_ref = _median_spacing(_inward_row_points(face, side, n_z, n_curve))

            for c in 1:n_curve
                _landing_confidence(points[c], landings[c], edge_spacing, scale_ref, patch, cfg) || continue
                if landings[c].kind === :crit_slice_vertex
                    all_points[flat_idxs[c]] = cell_point_vertex[landings[c].id]
                elseif landings[c].kind === :critical_point
                    all_points[flat_idxs[c]] = cell_point_anchor[landings[c].id]
                elseif landings[c].kind === :edge
                    entries = get!(pooled_edge_cols, landings[c].id, Tuple{Int,Vector{T}}[])
                    push!(entries, (flat_idxs[c], points[c]))
                end
            end
        end
    end

    # Pass 2: per crit-slice edge, jointly monotone-assign ALL pooled columns
    # from EVERY face at once. A per-face-only pass (each face's own contiguous
    # column run monotone-snapped independently, with no visibility into other
    # faces touching the SAME crit-slice edge) was tried first and found
    # insufficient: measured on the fixed-axis Taubin z=1.0648 boundary, two
    # different faces approaching the SAME crit-slice edge from different
    # directions picked internally-monotone but MUTUALLY INCONSISTENT target
    # subsets (e.g. one face -> targets [1,2,2,3], another -> [6,1]) -- each
    # individually valid, but not zippering together, leaving naked edges
    # despite 100% landing confidence. Pooling first and assigning once fixes
    # this by construction. Columns are rough-presorted by their own nearest
    # (unconstrained) target index before the monotone DP, giving it a
    # sensible starting order across faces with no other ordering relationship.
    for (eid, entries) in pooled_edge_cols
        tg = edge_targets[eid]
        rough_idx(pt) = argmin(T[norm(pt .- t) for t in tg])
        sort!(entries; by = e -> rough_idx(e[2]))
        pts_only = Vector{T}[e[2] for e in entries]
        assign = _monotone_snap_targets(pts_only, tg)
        for (i, (flat_idx, _)) in enumerate(entries)
            k = assign[i]
            all_points[flat_idx] = tg[k]
            haskey(target_flat_idx, (eid, k)) || (target_flat_idx[(eid, k)] = flat_idx)
        end
    end
    return target_flat_idx
end

"""
    _split_t_junctions(triangles::Vector{NTuple{3,Int}}, edge_polylines::Dict{Int,Vector{Int}}) -> Vector{NTuple{3,Int}}

Phase 9b: fixes T-junctions left by [`_snap_boundary_points!`](@ref) --
adjacent faces snapping to different, non-adjacent subsets of the same
crit-slice edge's targets means a triangle's boundary edge can have OTHER
faces' snapped points sitting on it as unconnected mesh vertices, which
`_naked_mesh_edges` would still count as cracks. `edge_polylines` maps each
crit-slice edge id to the ORDERED sequence of welded rep indices its USED
targets collapsed to (built by the caller from
[`_snap_boundary_points!`](@ref)'s `target_flat_idx` after clustering). For
every triangle edge `(p, q)` whose both endpoints lie on the SAME polyline
with intermediate polyline points between them, the triangle is fan-split
against the opposite vertex along that sub-polyline -- preserving the
original triangle's winding (the fan shares its cyclic order with the
original `(p, q, r)`). Iterates to a fixed point (bounded at 4 passes,
generous for the observed single-hop-per-edge case) since a fan split can in
principle expose a new triangle edge that itself needs splitting.
"""
function _split_t_junctions(triangles::Vector{NTuple{3,Int}}, edge_polylines::Dict{Int,Vector{Int}})
    on_polyline = Dict{Int,Tuple{Int,Int}}()
    for (eid, seq) in edge_polylines, (pos, r) in enumerate(seq)
        haskey(on_polyline, r) || (on_polyline[r] = (eid, pos))
    end
    isempty(on_polyline) && return triangles

    function split_one_pass(tris)
        out = NTuple{3,Int}[]
        changed = false
        for (a, b, c) in tris
            replaced = false
            for (p, q, r) in ((a, b, c), (b, c, a), (c, a, b))
                (haskey(on_polyline, p) && haskey(on_polyline, q)) || continue
                eid_p, pos_p = on_polyline[p]
                eid_q, pos_q = on_polyline[q]
                eid_p == eid_q || continue
                lo, hi = pos_p < pos_q ? (pos_p, pos_q) : (pos_q, pos_p)
                hi - lo <= 1 && continue
                seq = edge_polylines[eid_p]
                chain = pos_p < pos_q ? seq[lo:hi] : reverse(seq[lo:hi])
                for k in 1:(length(chain) - 1)
                    push!(out, (chain[k], chain[k+1], r))
                end
                replaced = true
                changed = true
                break
            end
            replaced || push!(out, (a, b, c))
        end
        return out, changed
    end

    tris = triangles
    for _ in 1:4
        tris, changed = split_one_pass(tris)
        changed || break
    end
    return tris
end

"""
    weld_mesh(faces::Vector{Face{T}}, patch::NamedTuple, cfg::HomotopyConfig{T}; incidence = nothing) where {T<:AbstractFloat}
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

# Keyword arguments (Phase 9b)
- `incidence`: `nothing` (default) is the exact pre-Phase-9b behavior --
  tolerance-only welding, byte-identical code path. A [`SurfaceIncidence`](@ref)
  (from `decompose_3d_surface(...; incidence = true)`) additionally snaps
  confident boundary-column landings onto shared crit-slice targets
  ([`_snap_boundary_points!`](@ref)) before clustering (pooled across ALL
  faces touching a given crit-slice edge, then jointly monotone-assigned --
  see that function's own docstring for why a per-face-only pass is
  insufficient), then fixes resulting T-junctions
  ([`_split_t_junctions`](@ref)) after.

  Measured on the fixed-axis Taubin fixture (2026-07): `_naked_mesh_edges`
  drops from 188 to 58 (69%). The two FOLD/POINT-type boundaries (cusp,
  lobe tips -- where all approaching columns converge to a single point)
  close COMPLETELY (0 naked edges, matching the sphere/ellipsoid pole
  controls exactly). The two multi-edge, multi-face EDGE-type boundaries
  (the singular notch, the saddle pair) reduce substantially but do NOT
  fully close: even with fully consistent cross-face assignment, a
  crit-slice edge's own sample points are only PARTIALLY touched by any
  swept column (measured: one edge's 8 samples had only 4 touched, e.g.
  `[1,2,4,7]`), because the crit-slice's sampling and the swept columns'
  arrival positions are two INDEPENDENT discretizations of the same curve --
  snapping onto existing targets cannot fabricate coverage that isn't there.
  Confirmed structural, not a resolution artifact: `edge_sample_density` of
  8/16/24 gave TOTAL naked-edge counts of 58/158/248 -- WORSE with density,
  since more independently-approaching columns spread over more targets
  makes exact multi-face reconciliation harder, not easier. Closing this
  fully requires a genuine conforming-triangulation/boundary-zippering step
  (merge every side's vertex set along the shared crit-slice curve and
  retriangulate against the union) -- deliberately DEFERRED as a future,
  separately-scoped phase rather than attempted here; see
  `test/test_taubin.jl`'s Phase 9b section for the exact measured numbers.

# Returns
A welded `GeometryBasics.Mesh` with `Point3f` vertices.
"""
function weld_mesh(
    faces::Vector{Face{T}},
    patch::NamedTuple,
    cfg::HomotopyConfig{T};
    incidence::Union{Nothing,SurfaceIncidence{T}} = nothing,
) where {T<:AbstractFloat}
    # Adjacent faces share boundary geometry (position-based merge, not id offset).
    all_points = Vector{T}[]
    provenance = Tuple{Int,Int}[]
    index_of = Dict{Tuple{Int,Int},Int}()
    for (fi, face) in enumerate(faces)
        for r in 1:size(face.mesh_vertices, 1)
            push!(all_points, face.mesh_vertices[r, :])
            push!(provenance, (fi, r))
            index_of[(fi, r)] = length(all_points)
        end
    end

    target_flat_idx = Dict{Tuple{Int,Int},Int}()
    if incidence !== nothing
        target_flat_idx = _snap_boundary_points!(all_points, index_of, faces, incidence, patch, cfg)
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

    if incidence !== nothing
        edge_polylines = Dict{Int,Vector{Int}}()
        for cs in incidence.crit_slices, e in cs.edges
            seq = Int[]
            for k in eachindex(e.sampled_points)
                idx = get(target_flat_idx, (e.id, k), 0)
                idx == 0 && continue
                r = membership[idx]
                (isempty(seq) || seq[end] != r) && push!(seq, r)
            end
            isempty(seq) || (edge_polylines[e.id] = seq)
        end
        fixed_triangles = _split_t_junctions(fixed_triangles, edge_polylines)
    end

    # Cluster at T-precision; narrow to Float32 only for the GeometryBasics.Mesh container.
    points3 = [GeometryBasics.Point3f(Float32(p[1]), Float32(p[2]), Float32(p[3])) for p in reps]
    tris = [GeometryBasics.TriangleFace{Int}(t[1], t[2], t[3]) for t in fixed_triangles]
    return GeometryBasics.Mesh(points3, tris)
end

"""
    _naked_mesh_edges(mesh::GeometryBasics.Mesh) -> Vector{Tuple{Int,Int}}

Watertightness audit utility (Phase 9a, permanent): the undirected mesh edges
that belong to exactly ONE triangle. A surface closed within the bounding box
should have none; the plain tolerance-only `weld_mesh` (`incidence = nothing`)
measurably leaves cracks along interior crit-slice boundaries (fixed-axis
Taubin baseline, verified deterministic across sessions: 188 naked edges —
10 at z=-1, 70 at z=1.0, 84 at z=1.0648, 24 at z=1.2367). Phase 9b's
incidence-based stitching (`weld_mesh(...; incidence = ...)`) closes this to
58 (0/26/32/0) on the same fixture: FULLY at the fold/point-type boundaries
(cusp, tips), partially at the multi-face edge-type boundaries (notch,
saddle pair) — see `weld_mesh`'s own Phase 9b docstring section for the
diagnosed root cause and why full closure there needs a separately-scoped
future phase, not more of the same mechanism. Both baselines are asserted
exactly in test/test_taubin.jl so any accidental shift is caught immediately.
"""
function _naked_mesh_edges(mesh::GeometryBasics.Mesh)
    counts = Dict{Tuple{Int,Int},Int}()
    for t in GeometryBasics.faces(mesh)
        a, b, c = Int(t[1]), Int(t[2]), Int(t[3])
        for (i, j) in ((a, b), (b, c), (c, a))
            key = i < j ? (i, j) : (j, i)
            counts[key] = get(counts, key, 0) + 1
        end
    end
    return [k for (k, n) in counts if n == 1]
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
- `incidence` (Phase 9a/9b): `false` (default) keeps the exact 4-tuple return
  and a byte-identical code path (no crit-slice decompositions, no extra
  solves, `weld_mesh` runs its original tolerance-only welding). `true`
  additionally decomposes every interior slab boundary an actual face
  touches (the crit-slices), assigns each face's boundary columns to
  crit-slice cells (nearest-wins, distances recorded), checks branch
  continuity between consecutive columns (flagged, not erred — see
  `SurfaceIncidence.continuity_ok`), and passes the result into `weld_mesh`,
  which additionally SNAPS confident landings onto shared crit-slice targets
  before clustering and fixes any resulting T-junctions afterward (see
  `weld_mesh`'s own Phase 9b section) — measured to reduce the fixed-axis
  Taubin fixture's `_naked_mesh_edges` count from 188 to 58, closing the
  fold/point-type boundaries completely and the multi-face edge-type
  boundaries partially (see `weld_mesh`'s docstring for the diagnosed
  reason full closure needs a separate future phase). Returns a 5-tuple
  whose last element is a [`SurfaceIncidence`](@ref). Composes with
  `projection` (incidence and stitching happen in the chart, then
  world-map). Measured cost on the fixed-axis Taubin fixture: ~+9 s on a
  ~16 s decompose (4 crit-slices + one extra critical-point solve for fold
  anchors; snapping/splitting themselves are comparatively negligible, no
  new solves).

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
    incidence::Bool = false,
) where {T<:AbstractFloat}
    length(F.variables) == 3 && length(F.expressions) == 1 || throw(ArgumentError(
        "decompose_3d_surface: expected F with exactly 1 equation in exactly 3 variables " *
        "(a raw surface, variables ordered [x_var, y_var, z_var] -- z LAST); " *
        "got $(length(F.expressions)) equation(s) in $(length(F.variables)) variable(s).",
    ))

    # Phase 8: change-of-coordinates wrapper. The chart run recurses with
    # projection = nothing, so the entire validated fixed-axis pipeline below
    # is reused unchanged; transforms live only at this boundary. Phase 9a:
    # incidence is computed in the CHART frame inside the recursion (slab
    # structure is chart-frame; Q is not recoverable by the caller for
    # :random), then world-mapped here alongside everything else.
    if projection !== nothing
        Q = _resolve_projection(projection, rng)
        F_chart = _rotate_system(F, Q)
        _verify_projection_ok(F_chart, cfg)
        cfg_chart = _chart_config(cfg, Q)
        if incidence
            vc, ec, fc, mc, inc = decompose_3d_surface(F_chart, cfg_chart; incidence = true)
            wv, we, wf, wm = _map_to_world(vc, ec, fc, mc, Q)
            return wv, we, wf, wm, _map_incidence_to_world(inc, Q)
        else
            chart_result = decompose_3d_surface(F_chart, cfg_chart)
            return _map_to_world(chart_result..., Q)
        end
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

    # Global id high-water marks. With incidence = false these track exactly the
    # `maximum(v.id for v in all_vertices)` values the pre-9a code computed (ids
    # grow monotonically), so slab renumbering is unchanged; with incidence = true
    # they additionally advance past lazily-created crit-slice cell ids, keeping
    # ONE collision-free namespace across both containers regardless of the
    # (lazy) order in which slabs and crit-slices are numbered.
    v_hwm = 0
    e_hwm = 0

    # Phase 9a/9b bookkeeping (all empty and untouched when incidence = false).
    crit_slice_cache = Dict{Int,CritSlice{T}}()
    fold_anchors = incidence ? _surface_critical_vertices(F, cfg) : NativeVertex{T}[]
    face_bottom_edges = Dict{Int,Vector{Int}}()
    face_top_edges = Dict{Int,Vector{Int}}()
    face_bottom_anchor = Dict{Int,Int}()
    face_top_anchor = Dict{Int,Int}()
    edge_faces = Dict{Int,Vector{Int}}()
    column_landings_bottom = Dict{Int,Vector{ColumnLanding{T}}}()
    column_landings_top = Dict{Int,Vector{ColumnLanding{T}}}()
    # Phase 9b: edge id -> local chord spacing (populated alongside crit-slice
    # creation below), and the branch-continuity results.
    edge_spacing = Dict{Int,T}()
    continuity_ok = Dict{Int,Bool}()
    continuity_violations_bottom = Dict{Int,Vector{Int}}()
    continuity_violations_top = Dict{Int,Vector{Int}}()

    # Crit-slices exist only at INTERIOR bounds (j in 2:length(z_bounds)-1);
    # bbox endpoints have no critical value to decompose. Computed lazily, once
    # per boundary, shared by both adjacent slabs.
    function _crit_slice_for!(j::Int)
        haskey(crit_slice_cache, j) && return crit_slice_cache[j]
        cs = _decompose_crit_slice(F, z_bounds[j], j, cfg, v_hwm, e_hwm)
        v_hwm = max(v_hwm, maximum(v.id for v in cs.vertices; init = v_hwm))
        e_hwm = max(e_hwm, maximum(e.id for e in cs.edges; init = e_hwm))
        for e in cs.edges
            edge_spacing[e.id] = _edge_spacing(e)
        end
        crit_slice_cache[j] = cs
        return cs
    end

    # Summarize one side's per-column landings into the face-level dicts.
    function _record_side!(face_id::Int, landings::Vector{ColumnLanding{T}},
                           edges_dict, anchor_dict, landings_dict)
        landings_dict[face_id] = landings
        edge_ids = Int[]
        anchor_id = 0
        for l in landings
            if l.kind === :edge
                l.id in edge_ids || push!(edge_ids, l.id)
            elseif l.kind === :critical_point && anchor_id == 0
                anchor_id = l.id
            end
        end
        edges_dict[face_id] = edge_ids
        anchor_dict[face_id] = anchor_id
        for eid in edge_ids
            faces_for_edge = get!(edge_faces, eid, Int[])
            face_id in faces_for_edge || push!(faces_for_edge, face_id)
        end
        return nothing
    end

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
        v_offset = v_hwm
        e_offset = e_hwm

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
        v_hwm = max(v_hwm, maximum(v.id for v in vertices_renumbered; init = v_hwm))
        e_hwm = max(e_hwm, maximum(e.id for e in edges_renumbered; init = e_hwm))

        append!(all_vertices, vertices_renumbered)
        append!(all_edges, edges_renumbered)

        # Boundary j = i is interior iff i >= 2; boundary j = i+1 is interior
        # iff i+1 <= length(z_bounds) - 1. Only compute crit-slices when this
        # slab actually produced faces.
        bottom_interior = i >= 2
        top_interior = i + 1 <= length(z_bounds) - 1
        cs_bottom = (incidence && bottom_interior && !isempty(edges_renumbered)) ? _crit_slice_for!(i) : nothing
        cs_top = (incidence && top_interior && !isempty(edges_renumbered)) ? _crit_slice_for!(i + 1) : nothing

        for edge in edges_renumbered
            # Face.boundary_edges come from already-renumbered edges.
            face = track_face(F, patch, edge, z_mid, z_bottom, z_top, next_face_id, cfg)
            push!(all_faces, face)

            if incidence
                # Landing rows recovered post-hoc from the swept grid: row 1 is
                # the z_bottom landing, row n_z the z_top landing (track_face's
                # documented layout) — FaceTracking needs no changes for this.
                n_z = 2 * cfg.midslice_sample_density + 1
                n_curve = size(face.mesh_vertices, 1) ÷ n_z
                bottom_landings = [Vector{T}(face.mesh_vertices[c, :]) for c in 1:n_curve]
                top_landings = [Vector{T}(face.mesh_vertices[(n_z - 1) * n_curve + c, :]) for c in 1:n_curve]
                bottom_assigned = _assign_landings(bottom_landings, cs_bottom, fold_anchors, z_bottom, cfg)
                top_assigned = _assign_landings(top_landings, cs_top, fold_anchors, z_top, cfg)
                _record_side!(face.id, bottom_assigned, face_bottom_edges, face_bottom_anchor, column_landings_bottom)
                _record_side!(face.id, top_assigned, face_top_edges, face_top_anchor, column_landings_top)
                _check_continuity!(face.id, bottom_landings, bottom_assigned, cs_bottom, edge_spacing,
                                    _inward_row_points(face, :bottom, n_z, n_curve), patch, cfg,
                                    continuity_ok, continuity_violations_bottom)
                _check_continuity!(face.id, top_landings, top_assigned, cs_top, edge_spacing,
                                    _inward_row_points(face, :top, n_z, n_curve), patch, cfg,
                                    continuity_ok, continuity_violations_top)
            end

            next_face_id += 1
        end
    end

    if !incidence
        mesh = weld_mesh(all_faces, patch, cfg)
        return all_vertices, all_edges, all_faces, mesh
    end

    # Phase 9b: SurfaceIncidence is assembled BEFORE weld_mesh so its snap-
    # unification / T-junction stitching can consume it.
    surface_incidence = SurfaceIncidence{T}(
        sort!(collect(values(crit_slice_cache)); by = cs -> cs.boundary_index),
        fold_anchors,
        face_bottom_edges,
        face_top_edges,
        face_bottom_anchor,
        face_top_anchor,
        edge_faces,
        column_landings_bottom,
        column_landings_top,
        continuity_ok,
        continuity_violations_bottom,
        continuity_violations_top,
    )
    mesh = weld_mesh(all_faces, patch, cfg; incidence = surface_incidence)
    return all_vertices, all_edges, all_faces, mesh, surface_incidence
end
