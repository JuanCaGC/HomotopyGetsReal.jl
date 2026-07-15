# src/PathTracking.jl
#
# Phase 4: extraction/hardening of the adaptive bidirectional path
# tracking engine that Phase 3's `connect_the_dots!` (src/Topology.jl)
# originally implemented inline as private helpers
# (`_is_near_singular_f64`, `_track_segment!`, `_track_outward!`).
#
# This phase is a REFACTOR, not new tracking logic: the adaptive
# bisection algorithm (retry near poor-quality/singular landings,
# accept-as-fallback once the step budget or interval width is
# exhausted -- see `_track_path_segment!`'s docstring) is unchanged from
# Phase 3. What changes:
#   - The engine is public (no leading underscore) and lives in its own
#     file, so Phase 5's future face-tracking can reuse it instead of
#     re-deriving a copy inside a new SurfaceDecomposition/FaceTracking
#     file.
#   - `y_start`/`y_final` are generalized from a hardcoded scalar
#     `ComplexF64` to `Vector{ComplexF64}`, so tracking more than one
#     variable at once (as Phase 5's face-tracking will need) does not
#     require a second refactor of this file. Phase 3's curve-tracking
#     call sites in `connect_the_dots!` wrap/unwrap a length-1 vector at
#     the boundary.
#   - `is_near_singular` checks BOTH `rank_tol` and `sv_thresh` (the
#     same two-part rule `Solver._classify_vertex_type` already uses),
#     closing a gap where Phase 3's `_is_near_singular_f64` only checked
#     the latter.
#   - `build_tracker` centralizes `ParameterHomotopy`/`Tracker`
#     construction and is the first place `cfg.path_tracker_precision`
#     (a `HomotopyConfig` field defined since Phase 1, unused until now)
#     actually gets threaded into HomotopyContinuation.jl.
#
# Precision boundary: unchanged from Phase 3. Every function here is
# Float64/ComplexF64-only -- `HomotopyConfig{T}` only ever appears in
# `build_tracker`'s signature (to read `cfg.path_tracker_precision`),
# never as the element type of any tracked quantity. Callers are
# responsible for Newton-polishing landing values to genuine
# `T`-precision (`Solver._newton_polish`) after tracking finishes; this
# file only ever deals with Float64 arithmetic.

"""
    is_near_singular(F::System, point::Vector{ComplexF64}, expected_rank::Int, rank_tol::Float64, sv_thresh::Float64)
        -> Bool

Return whether `point` is numerically singular for adaptive path tracking.

Used during Float64 tracking in [`track_path`](@ref) to trigger step bisection.
`point` must match `F.variables` in length and order. Returns `true` when the
Jacobian rank falls below `expected_rank` or the smallest singular value is below
`sv_thresh`.
"""
function is_near_singular(
    F::System,
    point::Vector{ComplexF64},
    expected_rank::Int,
    rank_tol::Float64,
    sv_thresh::Float64,
)::Bool
    # Float64-only probe (tracking is Float64-only); mirrors _classify_vertex_type's two-part rule.
    # HomotopyContinuation's `jacobian(F, x)` convenience wrapper does not
    # infer to a concrete matrix type on its own (confirmed via @inferred:
    # without these explicit `Matrix{ComplexF64}`/`Vector{Float64}` seals,
    # this function's return type infers as `Any`); sealing the type here
    # keeps everything downstream (and this function's own return type)
    # concretely inferred.
    J = Matrix{ComplexF64}(jacobian(F, point))
    svals = Vector{Float64}(svdvals(J))
    isempty(svals) && return true
    rank_deficient = count(>(rank_tol), svals) < expected_rank
    near_zero_sv = minimum(svals) < sv_thresh
    return rank_deficient || near_zero_sv
end

"""
    build_tracker(H_sys::System, x_start::Float64, cfg::HomotopyConfig{T}; compile::Symbol = :all) where {T<:AbstractFloat}
        -> (ph::ParameterHomotopy, tracker::Tracker)

Build a `ParameterHomotopy` and `Tracker` for [`track_path`](@ref) and [`track_bidirectional`](@ref).

Initialize both parameter slots to `x_start`; update them with `start_parameters!` and
`target_parameters!` on each step. Pass `compile = :none` when constructing many short-lived
systems (e.g. per-anchor face tracking); keep the default `:all` when reusing one system.
`cfg.path_tracker_precision` sets `TrackerOptions(min_step_size)`.
"""
function build_tracker(H_sys::System, x_start::Float64, cfg::HomotopyConfig{T}; compile::Symbol = :all) where {T<:AbstractFloat}
    # Fixed compile literal avoids Union return types from HC's :mixed heuristic.
    # :all amortizes for one H_sys per curve; :none wins for many fresh per-anchor systems.
    # path_tracker_precision maps to min_step_size (closest HC lever, not exact accuracy guarantee).
    ph = ParameterHomotopy(H_sys, [x_start], [x_start]; compile = compile)
    options = TrackerOptions(min_step_size = Float64(cfg.path_tracker_precision))
    tracker = Tracker(ph; options = options)
    return ph, tracker
end

"""
    _track_path_segment!(path, F, ph, tracker, x0, y0, x1, budget, expected_rank, rank_tol, sv_thresh, poor_acc_tol, min_width)
        -> Vector{ComplexF64}

Recursive core of [`track_path`](@ref) (unchanged algorithm from Phase
3's `_track_segment!`, generalized to vector-valued state `y0`/`y1`):
attempts a single `ParameterHomotopy` step from `x0` to `x1`. The
landing value is taken directly from `res.solution` whenever every
component is finite -- **not** gated on
`HomotopyContinuation.is_success(res)`, which is deliberately unreliable
right at a branch point (e.g. approaching a `Critical`-classified
curve vertex where `∂f/∂y = 0`): there, the tracked system has a
genuine multiple root, so the tracker's adaptive internal step size
legitimately shrinks to (near) zero and it terminates with
`:terminated_step_size_too_small` well *after* already landing
extremely close to the true root -- treating that as an unusable
failure would throw away a good answer.

If the landing value is non-finite, or `res.accuracy` (an estimate of
solution quality, backed by `poor_acc_tol` -- typically
`cfg.critical_point_tol`, see its own docstring: "did the solver
converge?") is worse than `poor_acc_tol`, **or** [`is_near_singular`](@ref)
fires at the landing point -- and the step-count `budget` (backed by
`cfg.max_path_steps`) is not exhausted and the interval is still wider
than `min_width` (typically `cfg.vertex_match_tol`) -- bisects and
recurses on each half instead of accepting a potentially-bad jump.
Crucially, a poor-accuracy landing triggers this retry **regardless**
of whether `is_near_singular` also fires (Julia's `||` short-circuits,
so `is_near_singular` is not even evaluated once `poor_quality` is
already `true`): both conditions are independently sufficient triggers,
never a joint requirement.

Once `budget` is exhausted or the interval has been bisected down to
`<= min_width`, the *last attempted* landing value is accepted as a
fallback and appended to `path` -- there is no "discard" outcome; every
call always produces some point, good or bad. The only downstream
defense against a genuinely bad landing point is the caller's own
vertex-matching against `cfg.vertex_match_tol` (see
`Topology._resolve_endpoint`), which will surface a bad landing as an
unexpected new `Artificial` vertex rather than silently corrupting an
existing one.

Appends every *accepted* `[x, real.(y)...]` point to `path` in the
order visited, and returns the landing `y` at `x1`.
"""
function _track_path_segment!(
    path::Vector{Vector{Float64}},
    F::System,
    ph::ParameterHomotopy,
    tracker::Tracker,
    x0::Float64,
    y0::Vector{ComplexF64},
    x1::Float64,
    budget::Base.RefValue{Int},
    expected_rank::Int,
    rank_tol::Float64,
    sv_thresh::Float64,
    poor_acc_tol::Float64,
    min_width::Float64,
)
    start_parameters!(ph, [x0])
    target_parameters!(ph, [x1])
    res = track(tracker, y0, 1.0, 0.0)
    budget[] -= 1

    candidate = res.solution
    finite = all(isfinite, candidate)
    y1 = finite ? candidate : y0
    poor_quality = !finite || res.accuracy > poor_acc_tol
    point = vcat(ComplexF64(x1), y1)
    problematic = poor_quality || is_near_singular(F, point, expected_rank, rank_tol, sv_thresh)

    if problematic && budget[] > 0 && abs(x1 - x0) > min_width
        xm = (x0 + x1) / 2
        ym = _track_path_segment!(
            path, F, ph, tracker, x0, y0, xm, budget, expected_rank, rank_tol, sv_thresh, poor_acc_tol, min_width,
        )
        return _track_path_segment!(
            path, F, ph, tracker, xm, ym, x1, budget, expected_rank, rank_tol, sv_thresh, poor_acc_tol, min_width,
        )
    else
        push!(path, vcat(x1, real.(y1)))
        return y1
    end
end

"""
    track_path(
        F::System, ph::ParameterHomotopy, tracker::Tracker,
        y_start::Vector{ComplexF64}, x_start::Float64, x_target::Float64,
        max_steps::Int, rank_tol::Float64, sv_thresh::Float64, poor_acc_tol::Float64, min_width::Float64,
    ) -> (y_final::Vector{ComplexF64}, path::Vector{Vector{Float64}})

Track one parameter-homotopy step from `x_start` to `x_target` with adaptive bisection.

`F` is the original defining system (for singularity checks only); `ph`/`tracker` carry the
parameterized system built by the caller. Returns the landing `y` at `x_target` and sampled
`[x, real.(y)...]` points along the way.
"""
function track_path(
    F::System,
    ph::ParameterHomotopy,
    tracker::Tracker,
    y_start::Vector{ComplexF64},
    x_start::Float64,
    x_target::Float64,
    max_steps::Int,
    rank_tol::Float64,
    sv_thresh::Float64,
    poor_acc_tol::Float64,
    min_width::Float64,
)
    # expected_rank = length(F.expressions), matching Solver entry points.
    path = Vector{Float64}[]
    budget = Ref(max(max_steps, 1))
    expected_rank = length(F.expressions)
    y_final = _track_path_segment!(
        path, F, ph, tracker, x_start, y_start, x_target, budget, expected_rank, rank_tol, sv_thresh, poor_acc_tol, min_width,
    )
    return y_final, path
end

"""
    track_bidirectional(
        F::System, ph::ParameterHomotopy, tracker::Tracker,
        y_mid::Vector{ComplexF64}, x_mid::Float64, x_left::Float64, x_right::Float64,
        max_steps::Int, rank_tol::Float64, sv_thresh::Float64, poor_acc_tol::Float64, min_width::Float64,
    ) -> (full_path::Vector{Vector{Float64}}, y_land_left::Vector{ComplexF64}, y_land_right::Vector{ComplexF64})

Track from a midslice point toward both interval endpoints and assemble the full path.

Calls [`track_path`](@ref) toward `x_left` and `x_right` independently. Returns
`reverse(path_left) ++ midpoint ++ path_right` plus both landing values for endpoint matching.
"""
function track_bidirectional(
    F::System,
    ph::ParameterHomotopy,
    tracker::Tracker,
    y_mid::Vector{ComplexF64},
    x_mid::Float64,
    x_left::Float64,
    x_right::Float64,
    max_steps::Int,
    rank_tol::Float64,
    sv_thresh::Float64,
    poor_acc_tol::Float64,
    min_width::Float64,
)
    y_land_left, path_left = track_path(
        F, ph, tracker, y_mid, x_mid, x_left, max_steps, rank_tol, sv_thresh, poor_acc_tol, min_width,
    )
    y_land_right, path_right = track_path(
        F, ph, tracker, y_mid, x_mid, x_right, max_steps, rank_tol, sv_thresh, poor_acc_tol, min_width,
    )

    full_path = Vector{Float64}[]
    for p in reverse(path_left)
        push!(full_path, p)
    end
    push!(full_path, vcat(x_mid, real.(y_mid)))
    for p in path_right
        push!(full_path, p)
    end

    return full_path, y_land_left, y_land_right
end
