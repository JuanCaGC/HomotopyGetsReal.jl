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

Finds the z-values at which the surface `F` cannot be locally expressed
as a smooth graph `z = g(x,y)` -- i.e. where `∂f/∂x = ∂f/∂y = 0` along
with `f = 0` -- via `Solver.compute_critical_points`'s existing
3-variable branch (`src/Solver.jl`, `ne == 1 && nv == 3`), which already
builds exactly this augmented system internally. No new solver logic.

**`F.variables` must be ordered `[x_var, y_var, z_var]`, z LAST** --
`compute_critical_points`'s 3-variable branch destructures `x, y, _ =
F.variables` and differentiates w.r.t. the first two, so getting this
order backwards silently computes critical points of the wrong
projection with no error thrown. This is the same requirement
[`build_patch_system`](@ref) has, and both must agree on the same `F`.

Deduplicates the resulting critical points' z-coordinates via
`Clustering.cluster_scalars(z_values, cfg.vertex_match_tol)` -- reusing
the exact same primitive and tolerance
[`Topology.decompose_1d_curve`](@ref) already uses for its own "which
x-values are actually distinct" problem, rather than reintroducing the
old prototype's bare `cluster_points_1d(z_values, 1e-5)` literal
(`prototipo_viejo_julia/Solver.jl:86`).

Returns a sorted `Vector{T}` of distinct z-values.
"""
function compute_critical_z_slices(F::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    length(F.variables) == 3 && length(F.expressions) == 1 || throw(ArgumentError(
        "compute_critical_z_slices: expected F with exactly 1 equation in exactly 3 variables " *
        "(a raw surface, variables ordered [x_var, y_var, z_var] -- z LAST); " *
        "got $(length(F.expressions)) equation(s) in $(length(F.variables)) variable(s).",
    ))
    crit_vertices = compute_critical_points(F, cfg)
    z_values = T[real(v.coordinates[3]) for v in crit_vertices]
    return cluster_scalars(z_values, cfg.vertex_match_tol)
end

"""
    slice_at_z(F::System, z_val::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> (vertices_3d::Vector{NativeVertex{T}}, edges_3d::Vector{Edge{T}})

Substitutes `z_var => Float64(z_val)` into the raw surface `F` to get a
2-variable curve, then calls [`Topology.decompose_1d_curve`](@ref)
**as-is** (confirmed in the Phase 5 investigation: `decompose_1d_curve`'s
own signature already accepts exactly this shape, zero Phase 3 changes
needed), and lifts the resulting 2D vertices/edges into 3-space by
appending `z_val` as a third coordinate to every `NativeVertex`
coordinate and every `Edge` sampled point.

This lifting is possible with zero `Types.jl` changes because
`NativeVertex{T}.coordinates`/`Edge{T}.sampled_points` were already
built length-generic in Phase 1 (confirmed: their keyword constructors
just do `Vector{Complex{T}}(undef, length(coordinates))`/`collect(p)`,
no hardcoded length-2 assumption).

Note: the returned vertex/edge `id`s are LOCAL to this one call (Phase
3's `decompose_1d_curve` restarts its own id counters near 1 every
call) -- callers concatenating results across multiple `slice_at_z`
calls (i.e. [`decompose_3d_surface`](@ref)) MUST renumber past the
running max before combining, exactly as `decompose_3d_surface` does.
"""
function slice_at_z(F::System, z_val::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    length(F.variables) == 3 && length(F.expressions) == 1 || throw(ArgumentError(
        "slice_at_z: expected F with exactly 1 equation in exactly 3 variables " *
        "(a raw surface); got $(length(F.expressions)) equation(s) in $(length(F.variables)) variable(s).",
    ))
    x_var, y_var, z_var = F.variables
    f = F.expressions[1]
    f_2d = subs(f, z_var => Float64(z_val))
    F_2d = System([f_2d], variables = [x_var, y_var])

    vertices_2d, edges_2d = decompose_1d_curve(F_2d, cfg)

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
quarter-points, so it costs two extra 2D-curve decompositions -- but only
lazily, memoized once per slab, and only ever paid on slabs that actually
need at least one retry (the overwhelming majority of slabs, including
every sphere/ellipsoid slab and 3 of this Taubin heart's own 5 slabs,
never retry at all and never pay this cost). If the reference scale
itself comes back exactly zero (only possible if the quarter-points are
themselves degenerate), the gradient gate is skipped for that slab rather
than dividing by zero -- documented as a known limitation of this
heuristic, not silently miscounted as "healthy".

# Known limitation: one confirmed, harmless false positive
Re-validating the FULL Taubin heart regression with this gate in place
surfaced a false positive on the `[1.0, 1.0648]` slab (the narrow,
4-x-critical-point slab discussed above): its naive midpoint's minimum
per-edge anchor gradient (`0.035`, at the thin inner sliver's pinch,
where a small inner loop and the large outer loop meet) compared against
the quarter-point reference's maximum (`3.94`, from the outer loop's own
much stronger part of the SAME curve) gives ratio `0.0089` -- just under
`cfg.z_mid_gradient_ratio_tol = 0.01` -- triggering one avoidable retry.
This is a MIN-candidate-vs-MAX-reference artifact of a curve with
multiple, legitimately very-differently-conditioned branches (inner
sliver vs. outer loop), not a true near-degeneracy: the rejected `z_mid`
already gave a perfectly good downstream sweep (max `|f|=2.4e-6`,
confirmed both before and after this gate existed), and the retry lands
on an equally fine `z_mid` (max `|f|=3.0e-6`), so no correctness issue
results -- just one avoidable `slice_at_z` call on this slab. Flagged
explicitly here rather than tuned away: shrinking the threshold to dodge
this specific ratio (`0.0089`) would erode most of the two-orders-of-
magnitude margin against the genuinely bad candidates (`0.0014`-`0.0057`)
found on the `[-1,1]` slab, and a more surgical fix (matching each
candidate edge to its own corresponding branch in the reference slice,
rather than a single global min/max) is a real design change, not a
threshold tweak, and is left for a future pass if a surface is found
where this false-positive rate actually matters.
"""
function _robust_slice_at_z(F::System, patch::NamedTuple, z_bottom::T, z_top::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    width = z_top - z_bottom
    max_multiple = max(1, floor(Int, T(0.45) / cfg.z_mid_retry_frac))

    _topology_suspect(vertices) = any(
        v -> v.v_type == Artificial && get(v.metadata, :origin, nothing) == :endpoint_fallback,
        vertices,
    ) && any(v -> v.v_type == Singular, vertices)

    _anchor_gradient_magnitude(edge::Edge{T}) = begin
        p = edge.sampled_points[1]
        a, b = patch_direction(patch, p[1], p[2], p[3], cfg)
        hypot(a, b)
    end

    ref_scale = Ref{Union{Nothing,T}}(nothing)
    function _reference_scale()
        ref_scale[] === nothing || return ref_scale[]
        scale = zero(T)
        for z_ref in (z_bottom + T(0.25) * width, z_top - T(0.25) * width)
            _, edges_ref = slice_at_z(F, z_ref, cfg)
            for e in edges_ref
                isempty(e.sampled_points) && continue
                scale = max(scale, _anchor_gradient_magnitude(e))
            end
        end
        ref_scale[] = scale
        return scale
    end

    _gradient_suspect(edges) = begin
        isempty(edges) && return false
        scale = _reference_scale()
        scale == zero(T) && return false
        cand_scale = minimum(_anchor_gradient_magnitude(e) for e in edges if !isempty(e.sampled_points))
        cand_scale < cfg.z_mid_gradient_ratio_tol * scale
    end

    _suspect(vertices, edges) = _topology_suspect(vertices) || _gradient_suspect(edges)

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
    weld_mesh(faces::Vector{Face{T}}, patch::NamedTuple, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> GeometryBasics.Mesh

Final assembly step: merges coincident vertices across every `Face`'s
own, independently-numbered `mesh_vertices`/`mesh_topology` into one
watertight mesh. Each `Face.mesh_topology` holds indices LOCAL to that
face's own `mesh_vertices` -- adjacent slabs' faces genuinely share
their boundary row (one slab's `z_top` sweep and the next slab's
`z_bottom` sweep both land on the same physical curve), so unlike
[`decompose_3d_surface`](@ref)'s vertex/edge id collision (a pure
namespace problem, fixed by renumbering, since no two slabs' curve
vertices are ever the same physical point), this genuinely needs
POSITION-BASED clustering, not just an offset.

# Explicit four-step procedure
1. **Flatten with provenance.** Every face's `mesh_vertices` rows are
   pushed into one flat `all_points::Vector{Vector{T}}`, recording
   `provenance[k] = (face_idx, local_row)` for each flattened point.
2. **Cluster globally.** `reps, membership =
   Clustering.cluster_points_indexed(all_points, cfg.vertex_match_tol)`
   -- `reps` becomes the merged mesh's global vertex list; `membership`
   maps each flattened point to its global index.
3. **Build the `(face_idx, local_row) -> global_idx` lookup** from
   `provenance` + `membership`.
4. **Remap every face's `mesh_topology`** through the lookup, dropping
   any triangle whose three remapped indices are not all distinct
   (legitimate: a swept grid can pinch to a point at a pole/tip, exactly
   as the old prototype already guards against at
   `prototipo_viejo_julia/SurfaceTopology.jl:217-224`).

# Winding-order normalization
Done HERE, globally, after remapping and degenerate-triangle filtering
(deliberately NOT in `FaceTracking.track_face`, which only emits the
fixed `(v1,v2,v3)`/`(v1,v3,v4)` diagonal split with no orientation
guarantee -- different edges are discovered independently by Phase 3's
midslice-first algorithm with no shared global orientation convention,
so a per-face fix could not guarantee cross-face consistency). For each
surviving triangle `(g1,g2,g3)`, computes the raw normal
`n = (p2-p1) x (p3-p1)` and the surface gradient `∇f = (fx,fy,fz)`
evaluated at `p1` (this specific vertex was chosen, not the centroid,
per the confirmed resolution: `p1` is already at hand from the cross
product computation, and a well-resolved mesh makes any one corner
equally representative of the local gradient direction for a pure sign
check). If `n . ∇f < 0`, the triangle is emitted as `(g1,g3,g2)`
(swapped winding); otherwise as `(g1,g2,g3)`. This gives every triangle
in the finished mesh a normal consistently aligned with `+∇f` (an
arbitrary but globally fixed convention).

Welds at genuine `T`-precision (see this file's module docstring for
why, contrasting with the old prototype's `Float32` truncation before
its own dedup decision); only narrows to `Float32`/`Point3f` for the
final `GeometryBasics.Mesh` container.

Not `@inferred`-clean -- but for a different reason than
[`build_face_tracker`](@ref)'s documented `Union`: `GeometryBasics.Mesh`'s
own constructor return type carries an existential `where _A` type
parameter (its internal view-storage type), which `Core.Compiler` cannot
resolve away regardless of how concretely-typed `weld_mesh`'s own inputs
are. Confirmed via `Base.return_types` in `scratch_phase5_check.jl`; an
upstream `GeometryBasics.jl` property, not something fixable here.
"""
function weld_mesh(faces::Vector{Face{T}}, patch::NamedTuple, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
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
            g1 != g2 && g2 != g3 && g3 != g1 && push!(global_triangles, (g1, g2, g3))
        end
    end

    fixed_triangles = NTuple{3,Int}[]
    for (g1, g2, g3) in global_triangles
        p1, p2, p3 = reps[g1], reps[g2], reps[g3]
        n = cross(p2 .- p1, p3 .- p1)
        gx, gy, gz = _gradient_at(patch, p1[1], p1[2], p1[3], cfg)
        push!(fixed_triangles, dot(n, T[gx, gy, gz]) < 0 ? (g1, g3, g2) : (g1, g2, g3))
    end

    points3 = [GeometryBasics.Point3f(Float32(p[1]), Float32(p[2]), Float32(p[3])) for p in reps]
    tris = [GeometryBasics.TriangleFace{Int}(t[1], t[2], t[3]) for t in fixed_triangles]
    return GeometryBasics.Mesh(points3, tris)
end

"""
    decompose_3d_surface(F::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> (vertices::Vector{NativeVertex{T}}, edges::Vector{Edge{T}}, faces::Vector{Face{T}}, mesh::GeometryBasics.Mesh)

Top-level Phase 5 orchestrator, the 3D analogue of
[`Topology.decompose_1d_curve`](@ref): computes critical z-slices, slices
+ decomposes + sweeps each resulting slab, and welds the result into a
single mesh -- running every pipeline step including the final weld, so
callers never need a manual follow-up call, mirroring
`decompose_1d_curve`'s own precedent of calling `sample_edge` internally
before returning.

**`F.variables` must be ordered `[x_var, y_var, z_var]`, z LAST** (see
[`compute_critical_z_slices`](@ref)/[`build_patch_system`](@ref)).

1. `z_crits = compute_critical_z_slices(F, cfg)`, filtered to
   `cfg.bbox_z`'s own range (a critical z-value outside the configured
   bounding box is not a slab boundary this decomposition should ever
   visit).
2. `z_bounds = sort(unique([cfg.bbox_z[1]; z_crits; cfg.bbox_z[2]]))`.
3. `patch = build_patch_system(F)` -- ONCE for the whole surface, passed
   both to [`_robust_slice_at_z`](@ref) (which reuses
   [`FaceTracking.patch_direction`](@ref) for its gradient-magnitude gate
   -- see that function's own docstring) and to
   [`FaceTracking.track_face`](@ref) below.
4. For each adjacent `(z_bottom, z_top)` pair: [`_robust_slice_at_z`](@ref)
   at the exact midpoint `(z_bottom + z_top) / 2` (matching
   [`Topology.compute_midslice`](@ref)'s own convention -- the old
   prototype's `0.4137` skew
   (`prototipo_viejo_julia/SurfaceTopology.jl:132`) is deliberately
   dropped, not ported to a new config field), falling back to a nearby
   perturbed `z_mid` only if the exact midpoint's slice comes back
   untrustworthy by either its vertex-type or gradient-magnitude gate
   (see [`_robust_slice_at_z`](@ref)'s own docstring for when/why), then
   renumber that slab's vertex/edge ids past the running max (see below)
   before calling [`FaceTracking.track_face`](@ref) once per resulting
   edge, using the SAME `z_mid` that `_robust_slice_at_z` actually used
   (not blindly re-deriving `(z_bottom+z_top)/2`), so the z-value sliced
   and the z-value swept-from always agree even when a retry fired.
5. [`weld_mesh`](@ref) over every produced `Face`.

# Vertex/edge id renumbering across slabs (explicit, not implied)
Every [`slice_at_z`](@ref) call restarts `decompose_1d_curve`'s own id
counters near 1, so simply concatenating each slab's `vertices_2d`/
`edges_2d` WOULD produce colliding ids across DIFFERENT slabs -- not a
"these are really the same point" merge candidate the way Phase 3's
`crit_vertices`/`bnd_vertices` were (every slab's `z_mid` is distinct,
so no two slabs' vertices are ever the same physical point), but a pure
NAMESPACE collision that would break any downstream `findfirst(v -> v.id
== ..., vertices)` lookup (the pattern `Topology.jl` itself uses
throughout) against the concatenated list. Fixed by offsetting BOTH id
namespaces independently past their own running max before
concatenating -- vertex ids and edge ids are tracked as two separate
counters, exactly as `decompose_1d_curve` itself already keeps them
separate within one call -- and shifting `left_vertex_id`/
`right_vertex_id` by the SAME vertex offset so they stay consistent with
the renumbered vertices. `Face.boundary_edges` is populated from the
already-renumbered edge, so it dereferences correctly against the
returned `edges`. This is the exact structural analogue of Phase 3's own
`offset = maximum(v.id for v in crit_vertices)` pattern
(`src/Topology.jl:358`), applied across slab calls instead of within one
call, and applied to two namespaces instead of one.
"""
function decompose_3d_surface(F::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    length(F.variables) == 3 && length(F.expressions) == 1 || throw(ArgumentError(
        "decompose_3d_surface: expected F with exactly 1 equation in exactly 3 variables " *
        "(a raw surface, variables ordered [x_var, y_var, z_var] -- z LAST); " *
        "got $(length(F.expressions)) equation(s) in $(length(F.variables)) variable(s).",
    ))

    z_bottom_bound, z_top_bound = cfg.bbox_z
    z_crits_raw = compute_critical_z_slices(F, cfg)
    z_crits = filter(z -> z_bottom_bound <= z <= z_top_bound, z_crits_raw)
    z_bounds = sort(unique(vcat(T[z_bottom_bound], z_crits, T[z_top_bound])))

    patch = build_patch_system(F)

    all_vertices = NativeVertex{T}[]
    all_edges = Edge{T}[]
    all_faces = Face{T}[]
    next_face_id = 1

    for i in 1:(length(z_bounds) - 1)
        z_bottom, z_top = z_bounds[i], z_bounds[i+1]
        vertices_2d, edges_2d, z_mid = _robust_slice_at_z(F, patch, z_bottom, z_top, cfg)

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
            face = track_face(F, patch, edge, z_mid, z_bottom, z_top, next_face_id, cfg)
            push!(all_faces, face)
            next_face_id += 1
        end
    end

    mesh = weld_mesh(all_faces, patch, cfg)
    return all_vertices, all_edges, all_faces, mesh
end
