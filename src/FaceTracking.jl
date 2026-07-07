# src/FaceTracking.jl
#
# Phase 5 (surface-sweeping half): the new tracking primitives needed to
# sweep a single 2D-slice curve across z into a 3D mesh "Face", building
# on top of Phase 4's PathTracking.jl engine but NOT simply reusing
# `track_path`/`track_bidirectional` unchanged -- the Phase 5
# investigation identified two genuine structural gaps:
#
#   1. `f(x,y,z)=0` at fixed z is one equation in two unknowns
#      (underdetermined): tracking a SPECIFIC point of that curve as z
#      sweeps needs a second "patch" equation pinning it, using the
#      local surface gradient `(f_y, -f_x)` at the anchor as the patch
#      direction (confirmed decision: gradient-based patch, not the old
#      prototype's origin-radial line, which degenerates for curves not
#      centered on the origin).
#   2. `track_path`/`track_bidirectional` are ENDPOINT trackers: their
#      adaptive bisection only records a point when it *accepts* a leaf,
#      so a smooth run can produce as few as 1-3 raw points -- fine for
#      curve endpoints, wrong for mesh-building, which needs a densely
#      and uniformly sampled z-grid. `track_dense_path` below is the new
#      sibling that walks a caller-supplied sequence of intermediate
#      z-targets, reusing `_track_path_segment!`'s adaptive
#      retry-on-poor-quality/near-singular logic between EACH consecutive
#      pair (so every hop keeps the same singularity-avoidance safety
#      net) while guaranteeing at least one accepted point per hop.
#
# Confirmed architecture decision (Option B): the per-anchor patch
# equation is built with LITERAL `x0,y0,a,b` coefficients (mirroring how
# `compute_midslice`/`connect_the_dots!` already bake literal `x_mid`
# values per interval), using `build_tracker(...; compile = :none)` --
# NOT a parametric patch with `x0,y0,a,b` threaded as extra
# ParameterHomotopy parameters. This was benchmarked directly: across 15
# structurally distinct systems, `compile = :all` cost ~5,700x what
# `compile = :none` cost (dominated by Julia re-specializing `track`'s
# call graph against each new `CompiledSystem` type), and `:none`'s
# per-step tracking throughput was not worse for the small (1-2
# equation) systems this pipeline builds. See `build_tracker`'s own
# docstring (src/PathTracking.jl) for the full writeup.
#
# ADAPTIVE RE-ANCHORING (discovered empirically while validating this
# file, see `_sweep_direction`'s docstring for the full derivation): a
# FIXED literal patch line is only guaranteed to keep intersecting the
# level curve as z sweeps away from the anchor when the surface's
# gradient is radially symmetric about the axis the curve shrinks toward
# (true for a sphere, false in general -- e.g. an ellipsoid). For a
# general surface, the fixed line can lose transversal intersection with
# the curve entirely before reaching z_bottom/z_top, which NO amount of
# `_track_path_segment!` bisection can fix (the system genuinely has no
# nearby solution beyond that point, not just a hard-to-resolve one).
# `_sweep_direction`/`_sweep_hop!` (used by `sweep_face_bidirectional`)
# guard against this with TWO complementary mechanisms, both needed
# (confirmed empirically -- neither alone was sufficient): (1) a
# dimensionless cosine-similarity check between the anchor's fixed
# gradient direction and the current local gradient
# (`cfg.patch_transversality_cos_tol`, a NEW config field -- reusing
# `jacobian_rank_tol`/`singular_value_threshold` was tried first and
# found to reliably fire too late, see `HomotopyConfig`'s docstring for
# the full empirical justification), which rebuilds a fresh, re-anchored
# `build_face_tracker` at the current landing point whenever the margin
# drifts too far; and (2) a residual-based bisection gate (reusing
# `cfg.critical_point_tol`, no new field), which recursively halves a hop
# whenever its landing fails the ground-truth `f≈0` surface-membership
# check -- needed because right at a genuine z-critical target the level
# curve degenerates to a point and NO choice of patch stays transversal,
# so only shrinking the step size (mirroring `_track_path_segment!`'s own
# bisection for the analogous 1D case) actually converges the tracked
# point to the true limit. Both use `_project_to_slice`/`_residual_at`
# (already needed for other reasons, see their own docstrings) to keep
# re-anchored points genuinely on-curve.
#
# Each of the two directions in `sweep_face_bidirectional` builds its OWN
# tracker from the same starting anchor (NOT shared between directions --
# adaptive re-anchoring in one direction must never affect the other) and
# then rebuilds ADAPTIVELY, and only as needed, via `_sweep_hop!`'s two
# quality gates above -- not unconditionally once per hop, which would
# defeat the whole point of Option B's cheap-per-anchor `compile = :none`
# design.
#
# Precision boundary: identical to every prior phase. `track_dense_path`
# and everything it calls are Float64/ComplexF64-only; `patch_direction`/
# `_gradient_at` are the only genuinely `T`-generic pieces here (pure
# symbolic-expression evaluation, same `bits = precision(T)` convention
# as `Solver.jacobian_rank_info`), since they determine the *direction*
# of a system that will itself be tracked in Float64.

"""
    build_patch_system(F::System)

Phase 5 entry point for the patch machinery: one-time, per-surface,
purely SYMBOLIC preparation -- takes only the raw surface, **no anchor
point** (evaluating the gradient at a specific anchor is a separate,
per-anchor responsibility, see [`patch_direction`](@ref)/
[`_gradient_at`](@ref)).

`F` must be a raw surface: `length(F.variables) == 3 &&
length(F.expressions) == 1`, with **`F.variables` ordered
`[x_var, y_var, z_var]`** (z LAST) -- the same convention
[`compute_critical_z_slices`](@ref)/`Solver.compute_critical_points`'s
3-variable branch already requires (`x, y, _ = F.variables`), since
getting this order backwards silently computes critical points of the
wrong projection with no error thrown.

Computes `∂f/∂x`, `∂f/∂y`, `∂f/∂z` once via `differentiate` (cheap,
symbolic, not per-anchor), and also builds `F_for_tracking` -- the SAME
expression `f`, but wrapped in a *second* `System` with variables
reordered to `[z_var, x_var, y_var]` (z FIRST). This second ordering is
required purely so that [`PathTracking._track_path_segment!`](@ref)'s
hardcoded `point = vcat(ComplexF64(x1), y1)` convention (swept parameter
first, then state -- exactly mirroring the plane-curve case where
`x_var` IS `F.variables[1]`) lines up correctly when the swept parameter
is `z` and the state is `(x,y)`. This is entirely independent from, and
does not change, the `[x_var,y_var,z_var]` order `F` itself must keep
for `compute_critical_z_slices`'s own augmentation convention -- the two
`System`s serve different, non-interacting purposes (`F_for_tracking` is
used *only* by [`is_near_singular`](@ref)'s rank check on the genuine
surface, via [`track_dense_path`](@ref); the per-anchor `H_sys` built by
[`build_face_tracker`](@ref) is what's actually numerically tracked).

Returns a `NamedTuple` with fields `f, fx, fy, fz, x_var, y_var, z_var,
F_for_tracking`.
"""
function build_patch_system(F::System)
    length(F.variables) == 3 && length(F.expressions) == 1 || throw(ArgumentError(
        "build_patch_system: expected F with exactly 1 equation in exactly 3 variables " *
        "(a raw surface, variables ordered [x_var, y_var, z_var] -- z LAST); " *
        "got $(length(F.expressions)) equation(s) in $(length(F.variables)) variable(s).",
    ))
    x_var, y_var, z_var = F.variables
    f = F.expressions[1]
    fx = differentiate(f, x_var)
    fy = differentiate(f, y_var)
    fz = differentiate(f, z_var)
    F_for_tracking = System([f], variables = [z_var, x_var, y_var])
    return (
        f = f, fx = fx, fy = fy, fz = fz,
        x_var = x_var, y_var = y_var, z_var = z_var,
        F_for_tracking = F_for_tracking,
    )
end

"""
    _gradient_at(patch::NamedTuple, x0::T, y0::T, z0::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> (fx_val::T, fy_val::T, fz_val::T)

Shared, `T`-generic gradient evaluator used by both
[`patch_direction`](@ref) (per-anchor patch construction) and
`SurfaceDecomposition.weld_mesh` (per-triangle winding-consistency
check). Evaluates `patch.fx, patch.fy, patch.fz` at `(x0,y0,z0)` via the
same low-level `evaluate(...; bits = precision(T))` convention
`Solver.jacobian_rank_info` already uses, so this is genuine
`T`-precision symbolic evaluation, not a Float64-only convenience
wrapper.

The raw `evaluate(...)` result is wrapped in an explicit outer
`Vector{Complex{T}}(...)` constructor call, mirroring the style
`Solver.jacobian_rank_info` uses for its own analogous
`Matrix{Complex{T}}(evaluate(...))` call. **This does NOT actually
achieve full `@inferred`-cleanliness** -- confirmed empirically while
validating this file: `Base.return_types` shows `_gradient_at` (and
hence [`patch_direction`](@ref), which merely destructures its result)
inferring as `Tuple{Any,Any,Any}`/`Tuple{Any,Any}`, and the SAME check
against `Solver.jacobian_rank_info` itself (never previously verified
with a direct `@inferred` test in `scratch_phase2_check.jl`, only
indirectly via callers that immediately re-seal through a concrete
constructor afterward, e.g. `NativeVertex{T}(...)`) shows the identical
`NamedTuple{...,<:Tuple{Any,Any}}` non-concrete result. This is a
genuine, pre-existing upstream limitation of
`HomotopyContinuation.jl`'s low-level `evaluate(...; bits=...)` (its own
return type is runtime-value-dependent, not just `compile`-keyword-value
-dependent like `build_tracker`'s case -- e.g. it returns `Vector{Float64}`
or `Vector{ComplexF64}` depending on whether the numeric result happens
to be real, which no amount of outer wrapping can retroactively narrow),
not something introduced by or fixable from this file. Downstream
callers are unaffected in practice for the same reason
`compute_critical_points` is unaffected by `jacobian_rank_info`'s
instability: both of `_gradient_at`'s callers immediately consume its
result into an already-concrete context (`patch_direction`'s own
`a_64 = Float64(a)`-style casts in [`build_face_tracker`](@ref), and
`weld_mesh`'s `T[gx, gy, gz]` typed-array-literal construction).
"""
function _gradient_at(patch::NamedTuple, x0::T, y0::T, z0::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    bits = precision(T)
    vars = [patch.x_var, patch.y_var, patch.z_var]
    vals = Complex{T}[Complex{T}(x0), Complex{T}(y0), Complex{T}(z0)]
    grad = Vector{Complex{T}}(evaluate([patch.fx, patch.fy, patch.fz], vars => vals; bits = bits))
    return T(real(grad[1])), T(real(grad[2])), T(real(grad[3]))
end

"""
    _residual_at(patch::NamedTuple, x0::T, y0::T, z0::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat} -> T

Shared, `T`-generic surface-membership residual: `|f(x0,y0,z0)|`,
evaluated via the same `evaluate(...; bits = precision(T))` convention as
[`_gradient_at`](@ref). Factored out of [`_project_to_slice`](@ref) so
that [`_sweep_hop!`](@ref)'s bisection quality gate (does this landing
actually satisfy `f≈0`, the ground-truth surface-membership invariant --
not merely "did the patch stay transversal", which
[`_sweep_hop!`](@ref) checks separately via cosine similarity) can reuse
the identical evaluation path rather than duplicating it.
"""
function _residual_at(patch::NamedTuple, x0::T, y0::T, z0::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    bits = precision(T)
    vars = [patch.x_var, patch.y_var, patch.z_var]
    vals = Complex{T}[Complex{T}(x0), Complex{T}(y0), Complex{T}(z0)]
    return T(real(evaluate(patch.f, vars => vals; bits = bits)))
end

"""
    _project_to_slice(patch::NamedTuple, x0::T, y0::T, z_val::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> (x::T, y::T)

Correctness fix discovered empirically while validating this file against
the unit-sphere test (`scratch_phase5_check.jl`): `Edge.sampled_points`
(as returned by `slice_at_z`, which passes `decompose_1d_curve`'s output
through unchanged) is the output of Phase 3's `sample_edge` --
**pure geometric linear interpolation between raw tracked points, with
zero re-projection onto the curve** (explicitly documented on
`sample_edge` itself: "zero `HomotopyContinuation.jl` involvement"). For
a smoothly curving arc sampled at only a few raw tracking points (e.g. a
quarter-circle with no bisection needed), the resulting chord-interpolated
points can land measurably OFF the actual curve -- fine for Phase 3's own
purposes (a display polyline), but fatal for [`track_face`](@ref), which
needs every `(x0,y0)` anchor to be an actual solution of `f(x,y,z_mid) =
0` before treating it as a `ParameterHomotopy` starting point.

Re-projects `(x0,y0)` onto the slice curve at fixed `z_val` via the
standard single-constraint Gauss-Newton correction (the minimal-norm
Newton step for one equation in two unknowns): repeatedly computes
`step = f(x,y,z_val) / (fx^2 + fy^2)` and updates `(x,y) -= step .*
(fx,fy)` (i.e. moves along `-∇f`, the direction of steepest descent
toward the zero level set) until `|f| <= cfg.critical_point_tol` --
reusing the existing "did the solver converge?" tolerance (see
`HomotopyConfig`'s docstring) rather than introducing a new one, exactly
like `Solver._newton_polish`'s own convergence check. Capped at 50
iterations (matching `_newton_polish`'s own cap).
"""
function _project_to_slice(patch::NamedTuple, x0::T, y0::T, z_val::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    x, y = x0, y0
    for _ in 1:50
        fval = _residual_at(patch, x, y, z_val, cfg)
        abs(fval) <= cfg.critical_point_tol && break
        fx_val, fy_val, _ = _gradient_at(patch, x, y, z_val, cfg)
        denom = fx_val^2 + fy_val^2
        denom == zero(T) && break
        step = fval / denom
        x -= step * fx_val
        y -= step * fy_val
    end
    return x, y
end

"""
    patch_direction(patch::NamedTuple, x0::T, y0::T, z0::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> (a::T, b::T)

Confirmed decision: the per-anchor patch equation `a*(x-x0) + b*(y-y0) =
0` uses the LOCAL SURFACE GRADIENT direction `(a,b) = (f_y(x0,y0,z0),
-f_x(x0,y0,z0))` at the anchor -- i.e. the line through the anchor
orthogonal to `(f_x,f_y)`'s in-plane component -- rather than the old
prototype's origin-radial line `x*y0 - y*x0 = 0`
(`prototipo_viejo_julia/SurfaceTopology.jl:39`), which degenerates for
curves not centered on the origin. This is the standard "coordinate
patch" trick from numerical algebraic geometry (Sommese-Wampler-style
witness-point tracking): a transversal line through the anchor is
guaranteed valid at the anchor itself by construction (the patch
equation is tautologically satisfied there), and the gradient direction
maximizes the local intersection angle with the curve, minimizing the
chance the patch itself introduces spurious near-tangency.
"""
function patch_direction(patch::NamedTuple, x0::T, y0::T, z0::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    fx_val, fy_val, _ = _gradient_at(patch, x0, y0, z0, cfg)
    return fy_val, -fx_val
end

"""
    build_face_tracker(patch::NamedTuple, x0::T, y0::T, z0::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> (H_sys::System, ph::ParameterHomotopy, tracker::Tracker)

Option B (confirmed): builds a FRESH, per-anchor `System` with
`x0,y0,a,b` baked in as LITERAL `Float64` coefficients (`a,b` from
[`patch_direction`](@ref) at `(x0,y0,z0)`), leaving only `z_var` as the
system's `ParameterHomotopy` parameter -- so this requires no changes to
`PathTracking.jl`'s existing scalar `x_start`/`x_target` contract.
Constructed via `build_tracker(...; compile = :none)`
(`src/PathTracking.jl`'s Phase 5 addition): per the benchmark documented
on `build_tracker`, baking literals per-anchor (rather than threading
`x0,y0,a,b` as extra homotopy parameters and building once) is only
feasible with `compile = :none`, since this function is called once per
anchor point (potentially hundreds to thousands of times per surface).
"""
function build_face_tracker(patch::NamedTuple, x0::T, y0::T, z0::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    a, b = patch_direction(patch, x0, y0, z0, cfg)
    x0_64, y0_64, a_64, b_64 = Float64(x0), Float64(y0), Float64(a), Float64(b)
    patch_eq = a_64 * (patch.x_var - x0_64) + b_64 * (patch.y_var - y0_64)
    H_sys = System([patch.f, patch_eq], variables = [patch.x_var, patch.y_var], parameters = [patch.z_var])
    ph, tracker = build_tracker(H_sys, Float64(z0), cfg; compile = :none)
    return H_sys, ph, tracker
end

"""
    track_dense_path(
        F::System, ph::ParameterHomotopy, tracker::Tracker,
        y_start::Vector{ComplexF64}, param_start::Float64, param_targets::Vector{Float64},
        max_steps_per_hop::Int, rank_tol::Float64, sv_thresh::Float64, poor_acc_tol::Float64, min_width::Float64,
    ) -> (y_final::Vector{ComplexF64}, dense_path::Vector{Vector{Float64}})

New sibling to [`PathTracking.track_path`](@ref) (see this file's module
docstring for why `track_path` itself is unsuited to mesh-building):
walks `param_targets` IN ORDER, calling
[`PathTracking._track_path_segment!`](@ref) once per consecutive pair
`(param_start_or_previous_target, param_targets[i])` into a FRESH,
per-hop scratch buffer -- reusing the EXACT SAME adaptive
retry-on-poor-quality/`is_near_singular` bisection `track_path` uses,
just applied between many small hops instead of once across the whole
interval. Regardless of how much intra-hop bisection fires, the LAST
point `_track_path_segment!` ever pushes for a given hop is always
exactly at that hop's own target (bisection only ever inserts EARLIER
points before the final landing, never changes what the final landing
value is computed at) -- so `track_dense_path` keeps only that final
per-hop point, guaranteeing `dense_path` has **exactly**
`length(param_targets)` entries, one per requested target, unlike
`track_path`'s endpoint-only accumulation (which may produce anywhere
from 1 point up to the full recursion depth for a single interval).

`F` is the *raw surface's* Jacobian-check system (i.e.
`patch.F_for_tracking` from [`build_patch_system`](@ref), NOT the
per-anchor augmented `H_sys` from [`build_face_tracker`](@ref)) --
`is_near_singular` checks genuine surface singularity, deliberately
decoupled from whether the artificial patch equation happens to be
poorly conditioned at some intermediate point.
"""
function track_dense_path(
    F::System,
    ph::ParameterHomotopy,
    tracker::Tracker,
    y_start::Vector{ComplexF64},
    param_start::Float64,
    param_targets::Vector{Float64},
    max_steps_per_hop::Int,
    rank_tol::Float64,
    sv_thresh::Float64,
    poor_acc_tol::Float64,
    min_width::Float64,
)
    expected_rank = length(F.expressions)
    y_cur = y_start
    p_cur = param_start
    dense_path = Vector{Vector{Float64}}(undef, length(param_targets))
    for (i, p_target) in enumerate(param_targets)
        hop_path = Vector{Float64}[]
        budget = Ref(max(max_steps_per_hop, 1))
        y_cur = _track_path_segment!(
            hop_path, F, ph, tracker, p_cur, y_cur, p_target, budget, expected_rank, rank_tol, sv_thresh, poor_acc_tol, min_width,
        )
        dense_path[i] = hop_path[end]
        p_cur = p_target
    end
    return y_cur, dense_path
end

"""
    _sweep_hop!(
        hop_path::Vector{Vector{Float64}}, patch::NamedTuple, cfg::HomotopyConfig{T},
        state::NamedTuple, y0::Vector{ComplexF64}, p0::Float64, p1::Float64,
        budget::Base.RefValue{Int}, tols::NamedTuple,
    ) where {T<:AbstractFloat} -> (state′::NamedTuple, y1::Vector{ComplexF64})

The bisection-and-re-anchoring engine behind [`_sweep_direction`](@ref).
`state` is `(ph, tracker, fx0_64, fy0_64)` -- the CURRENT patch's tracker
pair plus its anchor's raw gradient `(fx0,fy0) = _gradient_at(...)` at
the anchor (deliberately NOT reconstructed from
[`patch_direction`](@ref)'s `(a,b)=(fy0,-fx0)` return -- that would
require correctly inverting the swap-and-negate, and a sign error there
was caught empirically during validation: it silently negates the
computed `cos_angle`, causing the transversality check to fire on
literally every hop instead of only when genuinely needed. Calling
[`_gradient_at`](@ref) directly for this one purpose costs one extra
symbolic evaluation per re-anchor event, which is negligible next to
[`build_face_tracker`](@ref)'s own `compile=:none` construction cost).

Two independent, complementary quality gates are checked at every
accepted landing, addressing two DIFFERENT empirical failure modes found
validating this file against an asymmetric ellipsoid
(`x^2+4y^2+9z^2=1`, `scratch_phase5_check.jl` section 7):

1. **Transversality drift** (cosine check, `cfg.patch_transversality_cos_tol`,
   see `HomotopyConfig`'s docstring for the full derivation of why
   `H_sys`'s Jacobian determinant reduces to `-(fx*fx0+fy*fy0)`, the
   negative dot product of the current and anchor gradients restricted to
   `x,y`): tracking anchor `(-0.722,0.346)` toward the ellipsoid's pole at
   `z=-1/3`, the patch line eventually stops intersecting the level curve
   at all once the curve's tangent has rotated far enough relative to the
   anchor's fixed patch line (`|a*x0+b*y0| <= R*sqrt(a^2+b^2/4)` fails, in
   quadric-specific algebra terms). Re-anchoring at a still-good landing
   (`cos_angle` still comfortably positive) BEFORE this crossing avoids it
   entirely in the common case.
2. **Residual quality** (`cfg.critical_point_tol`, reusing the existing
   "did the solver converge?" tolerance -- see `HomotopyConfig`'s three
   easy-to-confuse tolerances section -- no new field): even with
   proactive re-anchoring, a COARSE `z_targets` grid (e.g.
   `midslice_sample_density = 8`) can still take one large single hop
   that crosses the *entire* remaining transversality margin in one step
   (confirmed empirically: with 8 hops toward `z=-1/3`, `cos_angle` is
   `0.72` after hop 6 -- above even a generous re-anchor threshold -- yet
   hop 7 already lands with `residual ≈ 0.076`), and this failure mode
   gets qualitatively WORSE, not better, right at a genuine
   `compute_critical_z_slices` target itself: AT `z=z_bottom` exactly,
   `(fx,fy)` vanishes identically (that is the z-slice's defining
   property), so the level curve has shrunk to a single point and
   `H_sys`'s Jacobian is genuinely singular there for EVERY choice of
   patch, not just a poorly-re-anchored one. No amount of proactive
   re-anchoring alone fixes this -- what actually fixes it is the same
   remedy `_track_path_segment!` already uses for the analogous 1D
   problem: **bisect the hop itself** when the landing's actual `|f|`
   residual exceeds `cfg.critical_point_tol`, inserting a midpoint
   `z`-target and resolving each half recursively (mirroring
   `_track_path_segment!`'s own `xm = (x0+x1)/2` pattern), which
   empirically converges the tracked point arbitrarily close to the true
   limiting point `(0,0)` as the hop width shrinks toward it, exactly as
   `_track_path_segment!` converges toward genuine 1D critical points.

Recursion is bounded by the shared `budget` (seeded from
`cfg.max_path_steps`, one budget per top-level `z_targets` entry, reset
by [`_sweep_direction`](@ref) -- same convention as
[`track_dense_path`](@ref)'s own per-hop budget) and by `tols.min_width64`
(`cfg.vertex_match_tol`, reused rather than a new bare literal, exactly
like `_track_path_segment!`'s own `min_width` parameter): once a hop can
no longer be halved or the budget is exhausted, the landing is accepted
regardless of residual (an explicit, documented fallback -- not a silent
one, since [`_sweep_direction`](@ref)'s docstring records this as a
known boundary case), matching `_track_path_segment!`'s own fallback
philosophy for genuine critical points it cannot fully resolve either.

Only the FINAL point resolved for the top-level `(p0,p1)` interval is
pushed onto `hop_path` -- `_sweep_direction` uses `hop_path[end]` for its
per-target `dense_path` entry -- but bisection may push intermediate
points first (mirroring `_track_path_segment!`/`track_dense_path`'s own
"only the last pushed point matters" discard convention for interior
recursion levels).
"""
function _sweep_hop!(
    hop_path::Vector{Vector{Float64}},
    patch::NamedTuple,
    cfg::HomotopyConfig{T},
    state::NamedTuple,
    y0::Vector{ComplexF64},
    p0::Float64,
    p1::Float64,
    budget::Base.RefValue{Int},
    tols::NamedTuple,
) where {T<:AbstractFloat}
    y1, seg = track_dense_path(
        patch.F_for_tracking, state.ph, state.tracker, y0, p0, Float64[p1],
        cfg.max_path_steps, tols.rank_tol64, tols.sv_thresh64, tols.poor_acc_tol64, tols.min_width64,
    )
    budget[] -= 1
    x_land, y_land = seg[1][2], seg[1][3]
    resid = _residual_at(patch, T(x_land), T(y_land), T(p1), cfg)
    poor = !isfinite(resid) || abs(resid) > cfg.critical_point_tol

    if poor && budget[] > 0 && abs(p1 - p0) > tols.min_width64
        pmid = (p0 + p1) / 2
        state_mid, ymid = _sweep_hop!(hop_path, patch, cfg, state, y0, p0, pmid, budget, tols)
        return _sweep_hop!(hop_path, patch, cfg, state_mid, ymid, pmid, p1, budget, tols)
    end

    push!(hop_path, seg[1])
    fx_val, fy_val, _ = _gradient_at(patch, T(x_land), T(y_land), T(p1), cfg)
    fx_64, fy_64 = Float64(fx_val), Float64(fy_val)
    denom = hypot(fx_64, fy_64) * hypot(state.fx0_64, state.fy0_64)
    cos_angle = denom == 0.0 ? 0.0 : (fx_64 * state.fx0_64 + fy_64 * state.fy0_64) / denom
    if cos_angle < tols.cos_tol64
        anchor_x, anchor_y = _project_to_slice(patch, T(x_land), T(y_land), T(p1), cfg)
        _, ph_new, tracker_new = build_face_tracker(patch, anchor_x, anchor_y, T(p1), cfg)
        fx0_new, fy0_new, _ = _gradient_at(patch, anchor_x, anchor_y, T(p1), cfg)
        new_state = (ph = ph_new, tracker = tracker_new, fx0_64 = Float64(fx0_new), fy0_64 = Float64(fy0_new))
        return new_state, ComplexF64[ComplexF64(anchor_x), ComplexF64(anchor_y)]
    end
    return state, y1
end

"""
    _sweep_direction(patch::NamedTuple, x0::T, y0::T, z_start::T, z_targets::Vector{Float64}, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> Vector{Vector{Float64}}

**Adaptive re-anchoring + bisection**, added after an empirical failure
discovered while validating this file against an asymmetric ellipsoid
(`x^2 + 4y^2 + 9z^2 = 1`, `scratch_phase5_check.jl` section 7):
[`build_face_tracker`](@ref)'s patch line is LITERAL -- fixed at the
anchor's starting `(x0,y0)` and gradient -- and stays valid only as long
as it keeps genuinely intersecting the level curve as `z` sweeps away.
For a sphere this always holds (its gradient is exactly radial, so the
patch line always passes through the axis the circle shrinks toward --
an algebraic coincidence of radial symmetry, not a general guarantee).
For a general surface it does NOT, and worse, the failure gets sharper
rather than gentler right at a genuine `compute_critical_z_slices`
target itself, where the level curve degenerates to a single point. See
[`_sweep_hop!`](@ref)'s docstring for the full two-part derivation (a
transversality-drift cosine check AND a residual-based bisection gate --
genuinely different failure modes, both observed empirically, neither
sufficient alone) and `HomotopyConfig`'s docstring for why the cosine
check specifically needed a new field rather than reusing
`jacobian_rank_tol`/`singular_value_threshold`.

# Algorithm
Walks `z_targets` one at a time (rather than delegating the whole list to
a single [`track_dense_path`](@ref) call, so this function can inspect
and react to each intermediate landing): each top-level hop is resolved
by [`_sweep_hop!`](@ref), which tracks it via [`track_dense_path`](@ref)
and then either accepts the landing (updating the current patch's
tracker state if the transversality margin has drifted) or bisects the
hop and recurses, exactly mirroring `_track_path_segment!`'s own
recursive bisection pattern for the analogous 1D problem. Only
`hop_path[end]` -- the actually-resolved landing at each requested
`z_targets[i]` -- is kept per target, guaranteeing `dense_path` has
**exactly** `length(z_targets)` entries, same contract as
[`track_dense_path`](@ref) itself.

Reconstruction (a fresh [`build_face_tracker`](@ref) call) only happens
when the cosine check actually fires, and bisection only happens when
the residual check fires -- both confirmed empirically RARE in the
well-conditioned interior of a sweep and concentrated near genuine
z-critical targets (poles) on the ellipsoid regression case -- so this
stays a genuinely ADAPTIVE strategy, not unconditional per-hop
rebuilding, keeping total `compile = :none` construction cost close to
the original one-per-anchor design in the common case, per Phase 5's own
`build_tracker` benchmark.
"""
function _sweep_direction(
    patch::NamedTuple,
    x0::T,
    y0::T,
    z_start::T,
    z_targets::Vector{Float64},
    cfg::HomotopyConfig{T},
) where {T<:AbstractFloat}
    tols = (
        rank_tol64 = Float64(cfg.jacobian_rank_tol),
        sv_thresh64 = Float64(cfg.singular_value_threshold),
        poor_acc_tol64 = Float64(cfg.critical_point_tol),
        min_width64 = Float64(cfg.vertex_match_tol),
        cos_tol64 = Float64(cfg.patch_transversality_cos_tol),
    )

    _, ph, tracker = build_face_tracker(patch, x0, y0, z_start, cfg)
    fx0, fy0, _ = _gradient_at(patch, x0, y0, z_start, cfg)
    state = (ph = ph, tracker = tracker, fx0_64 = Float64(fx0), fy0_64 = Float64(fy0))
    y_state = ComplexF64[ComplexF64(x0), ComplexF64(y0)]
    p_cur = Float64(z_start)

    dense_path = Vector{Vector{Float64}}(undef, length(z_targets))
    for (i, p_target) in enumerate(z_targets)
        hop_path = Vector{Float64}[]
        budget = Ref(max(cfg.max_path_steps, 1))
        state, y_state = _sweep_hop!(hop_path, patch, cfg, state, y_state, p_cur, p_target, budget, tols)
        dense_path[i] = hop_path[end]
        p_cur = p_target
    end
    return dense_path
end

"""
    sweep_face_bidirectional(
        F::System, patch::NamedTuple,
        x0::T, y0::T, z_mid::T, z_bottom::T, z_top::T,
        cfg::HomotopyConfig{T},
    ) where {T<:AbstractFloat}
        -> (dense_path_down::Vector{Vector{Float64}}, dense_path_up::Vector{Vector{Float64}})

The z-sweep analogue of [`PathTracking.track_bidirectional`](@ref):
sweeps ONE anchor `(x0,y0)` from `z_mid` toward `z_bottom` and toward
`z_top`, each via [`_sweep_direction`](@ref) with
`cfg.midslice_sample_density` intermediate targets EXCLUDING the shared
start (`z_mid`) and INCLUDING the far end (`z_bottom`/`z_top`
respectively, tracked directly, no artificial pullback margin -- see
`sweep_face_bidirectional`'s dropped-epsilon rationale below).

Each direction independently starts from the SAME anchor `(x0,y0,z_mid)`
and may adaptively re-anchor partway through its own sweep (see
[`_sweep_direction`](@ref)) -- the two directions never share a tracker
instance (each calls [`build_face_tracker`](@ref) freshly), so there is
no stale-state carryover between them.

`_sweep_direction`'s internal re-anchoring calls
[`build_face_tracker`](@ref) (hence `build_tracker(...; compile =
:none)`) an adaptive, usually-small number of times per direction --
confirmed empirically fully `@inferred`-clean end-to-end in
`scratch_phase5_check.jl` despite this, for the same union-splitting
reason documented on `build_tracker` itself.

**No `ε` pullback from `z_bottom`/`z_top`.** Phase 3's
`connect_the_dots!`/`track_bidirectional` already track directly to
genuinely critical/singular target values with no pullback (e.g.
`scratch_phase4_check.jl` tracks directly to a real node at `x=0` and
converges), relying entirely on `is_near_singular`'s bisection retry.
`z_bottom`/`z_top` play exactly the same role here (each is, by
construction, either a genuine z-critical value from
[`compute_critical_z_slices`](@ref) or a `cfg.bbox_z` boundary) -- there
is no principled reason the z-sweep needs a safety margin the
already-proven 1D engine doesn't, and introducing one would reintroduce
exactly the kind of unconfigured, ad-hoc literal
(`prototipo_viejo_julia/SurfaceTopology.jl:46`'s `ε = 1e-4`) this
rebuild exists to eliminate. [`_sweep_direction`](@ref)'s adaptive
re-anchoring is the mechanism that makes tracking all the way to a
genuine limiting point (like a pole) actually work on general surfaces,
rather than merely failing gracefully.
"""
function sweep_face_bidirectional(
    F::System,
    patch::NamedTuple,
    x0::T,
    y0::T,
    z_mid::T,
    z_bottom::T,
    z_top::T,
    cfg::HomotopyConfig{T},
) where {T<:AbstractFloat}
    n_side = cfg.midslice_sample_density
    z_mid64 = Float64(z_mid)

    targets_down = collect(range(z_mid64, Float64(z_bottom); length = n_side + 1))[2:end]
    dense_path_down = _sweep_direction(patch, x0, y0, z_mid, targets_down, cfg)

    targets_up = collect(range(z_mid64, Float64(z_top); length = n_side + 1))[2:end]
    dense_path_up = _sweep_direction(patch, x0, y0, z_mid, targets_up, cfg)

    return dense_path_down, dense_path_up
end

"""
    track_face(
        F::System, patch::NamedTuple, edge::Edge{T},
        z_mid::T, z_bottom::T, z_top::T, face_id::Int, cfg::HomotopyConfig{T},
    ) where {T<:AbstractFloat}
        -> Face{T}

Sweeps an entire `Edge{T}` (already resampled to `cfg.edge_sample_density`
points by Phase 3's `sample_edge`) across the full slab `[z_bottom,
z_top]`, assembling one `Face{T}` whose `mesh_vertices` is a
`(2*cfg.midslice_sample_density + 1) x cfg.edge_sample_density` grid of
points all lying on the surface (`f ≈ 0` by construction, since every
row comes from a tracked point of `patch.f`).

# Row/column indexing convention (fixed, see architecture resolution)
- Columns `c = 1..edge_sample_density` index `edge.sampled_points`
  (curve arclength order, unchanged from Phase 3).
- Rows `r = 1..n_z` index the combined z-sweep, `n_z = 2*n_side + 1`
  where `n_side = cfg.midslice_sample_density`: row 1 is `z_bottom`, row
  `n_side+1` is the EXACT `z_mid` anchor (not a tracked value -- the
  anchor coordinates themselves), row `n_z` is `z_top`, built as
  `reverse(dense_path_down) ++ [mid_row] ++ dense_path_up` -- the direct
  z-axis analogue of `track_bidirectional`'s own
  `reverse(path_left) ++ [midpoint] ++ path_right` combination.
  [`sweep_face_bidirectional`](@ref)'s dense paths are `[z,x,y]`-ordered
  (matching `_track_path_segment!`'s `vcat(swept_param, state)`
  convention); each row is reordered to `[x,y,z]` here, matching this
  codebase's standing coordinate convention, before being stored.
- Flat storage is row-major: point `(r,c)` lives at
  `mesh_vertices[(r-1)*edge_sample_density + c, :]`.

Each `edge.sampled_points[c]` is first re-projected onto the actual
slice curve at `z_mid` via [`_project_to_slice`](@ref) before being used
as a sweep anchor -- see that function's docstring for why this is
necessary (Phase 3's `sample_edge` output is only an approximate,
linearly-interpolated polyline, not literally on the curve).

# Triangulation convention (fixed, see architecture resolution)
For each grid cell with corners `v1=(r,c), v2=(r+1,c), v3=(r+1,c+1),
v4=(r,c+1)`, the diagonal always runs `v1`-`v3`, splitting into
`(v1,v2,v3)` and `(v1,v3,v4)` (matching
`prototipo_viejo_julia/SurfaceTopology.jl:211-224`'s existing
convention). Triangles degenerate under THIS face's own local indices
(e.g. two corners of a swept quad coincide from a pinch) are dropped
here; **winding-order normalization (aligning normals with `∇f`) is
deliberately NOT done here** -- it happens once, globally, in
`SurfaceDecomposition.weld_mesh`, after cross-face index remapping, so
that every triangle in the finished mesh (not just within one face) is
consistently oriented.
"""
function track_face(
    F::System,
    patch::NamedTuple,
    edge::Edge{T},
    z_mid::T,
    z_bottom::T,
    z_top::T,
    face_id::Int,
    cfg::HomotopyConfig{T},
) where {T<:AbstractFloat}
    n_curve = length(edge.sampled_points)
    n_curve >= 1 || throw(ArgumentError("track_face: edge $(edge.id) has no sampled points to sweep"))
    n_side = cfg.midslice_sample_density
    n_z = 2 * n_side + 1

    mesh_vertices = Matrix{T}(undef, n_z * n_curve, 3)

    for (c, pt) in enumerate(edge.sampled_points)
        x0, y0 = _project_to_slice(patch, pt[1], pt[2], z_mid, cfg)
        dense_down, dense_up = sweep_face_bidirectional(F, patch, x0, y0, z_mid, z_bottom, z_top, cfg)
        # track_dense_path guarantees exactly n_side points per direction
        # (one per requested target -- see its own docstring).
        length(dense_down) == n_side && length(dense_up) == n_side || throw(ErrorException(
            "track_face: expected exactly $(n_side) sweep points per direction for edge $(edge.id), " *
            "got down=$(length(dense_down)), up=$(length(dense_up))",
        ))

        rows = Vector{Vector{T}}(undef, n_z)
        for (k, p) in enumerate(reverse(dense_down))
            rows[k] = T[p[2], p[3], p[1]]  # [z,x,y] -> [x,y,z]
        end
        rows[n_side + 1] = T[x0, y0, z_mid]
        for (k, p) in enumerate(dense_up)
            rows[n_side + 1 + k] = T[p[2], p[3], p[1]]
        end

        for (r, row) in enumerate(rows)
            mesh_vertices[(r - 1) * n_curve + c, :] = row
        end
    end

    triangles = NTuple{3,Int}[]
    for r in 1:(n_z - 1), c in 1:(n_curve - 1)
        v1 = (r - 1) * n_curve + c
        v2 = r * n_curve + c
        v3 = r * n_curve + (c + 1)
        v4 = (r - 1) * n_curve + (c + 1)
        v1 != v2 && v2 != v3 && v3 != v1 && push!(triangles, (v1, v2, v3))
        v1 != v3 && v3 != v4 && v4 != v1 && push!(triangles, (v1, v3, v4))
    end

    mesh_topology = Matrix{Int}(undef, length(triangles), 3)
    for (i, t) in enumerate(triangles)
        mesh_topology[i, 1] = t[1]
        mesh_topology[i, 2] = t[2]
        mesh_topology[i, 3] = t[3]
    end

    return Face{T}(
        id = face_id,
        mid_slice_z = z_mid,
        boundary_edges = [edge.id],
        mesh_vertices = mesh_vertices,
        mesh_topology = mesh_topology,
    )
end
