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

Float64-only singularity probe used *during* Float64 path tracking
(adaptive step-size reduction in [`track_path`](@ref)), as distinct from
[`jacobian_rank_info`](@ref) (used elsewhere for genuine `T`-precision
classification of already-matched endpoints). Since tracking itself is
Float64-only by construction (see module docstring), there is nothing
to gain from a `bits`-generic Jacobian evaluation here, so this uses
HomotopyContinuation's own convenience `jacobian(F, x)` wrapper directly.

`point` must match `F.variables` in length and order (unlike Phase 3's
`_is_near_singular_f64`, which hardcoded a 2-element `[x, y]` point,
this works for any system shape).

Returns `true` (i.e. "treat this as a singularity, refine the step")
iff *either*:
- the numerical rank of `jacobian(F, point)` (columns whose singular
  value exceeds `rank_tol`) falls short of `expected_rank`, or
- the smallest singular value is below `sv_thresh`,

exactly mirroring [`Solver._classify_vertex_type`](@ref)'s two-part
rule (rank-deficiency *or* near-zero smallest singular value), so
mid-tracking singularity detection and post-hoc vertex classification
agree on what "singular" means.
"""
function is_near_singular(
    F::System,
    point::Vector{ComplexF64},
    expected_rank::Int,
    rank_tol::Float64,
    sv_thresh::Float64,
)::Bool
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

Constructs the `ParameterHomotopy`/`Tracker` pair used by
[`track_path`](@ref)/[`track_bidirectional`](@ref), centralizing what
Phase 3's `connect_the_dots!` previously built inline -- so that every
current and future caller (Phase 5's face-tracking included) gets
consistent precision wiring for free.

`ph` is initialized with `start_parameters = target_parameters =
[x_start]` (the caller is expected to move it via `start_parameters!`/
`target_parameters!` on each tracking step, exactly as before), built
with `compile = compile` (default `:all`, matching Phase 3/4's original
choice and every existing call site's behavior unchanged) rather than
leaving `HomotopyContinuation`'s own `COMPILE_DEFAULT[] = :mixed`
heuristic in place: `fixed(F; compile = :mixed)` picks between
`InterpretedSystem`/`MixedSystem`/`CompiledSystem` based on a runtime
system-size heuristic, which makes `ParameterHomotopy(H_sys, ...)`'s own
return type a `Union` of those three (confirmed via `@inferred`) --
forcing *any single fixed literal* (`:all` or `:none`) instead always
produces one single concrete system type, which is what makes
`build_tracker`/`track_path`/`track_bidirectional` all `@inferred`-clean
regardless of which literal is chosen (see the `@inferred` note below).

**Why this is now a caller-chosen keyword (Phase 5 addition):** Phase
3/4's curve-tracking call sites (`connect_the_dots!`) build one `H_sys`
per curve and reuse it for every branch/direction of that same curve, so
`:all`'s one-time `CompiledSystem`-specialization cost is paid once and
amortized -- genuinely negligible, as originally documented. Phase 5's
face-tracking (`FaceTracking.build_face_tracker`) instead builds a
*fresh, literal-coefficient* `H_sys` **per anchor point, per sweep
direction** (potentially hundreds to thousands per surface), so that
one-time cost is paid over and over with no amortization. This was
benchmarked directly (not guessed): across 15 structurally distinct
systems, `compile = :all`'s total construct+first-track cost was
~5,700x that of `compile = :none` (dominated by Julia re-specializing
`track`'s whole call graph against each new `CompiledSystem` type, not
by `HomotopyContinuation`'s own system-compilation step), while
`compile = :none`'s per-step tracking throughput on an already-warm
system was *not* worse than `:all`'s for the small (1-2 equation)
systems this pipeline builds. Hence: Phase 3's existing call site keeps
the default `:all` (unchanged behavior), and `FaceTracking.jl`'s
per-anchor construction explicitly passes `compile = :none`.

Note on `@inferred`: even with `compile` fixed to a literal, `build_tracker`
itself is **not** `@inferred`-clean -- `HomotopyContinuation.jl`'s own
`fixed`/`ParameterHomotopy` dispatch chain branches on the runtime
*value* of the `compile` keyword (`compile == true || compile == :all`,
`elseif compile == false || compile == :none`, ...), not just its type
(`Union{Bool,Symbol}`), so `Core.Compiler` cannot narrow past that
`if`/`elseif` chain regardless of which literal we pass at any given
call site; this is the identical category of upstream instability
`solve()` itself has for the same reason, and it is identical for
`:all` and `:none` alike -- neither is "more inferred" than the other.
This is a HomotopyContinuation.jl-internal limitation, not something
`build_tracker` introduces or can paper over -- see
`scratch_phase4_check.jl` for a documented `Core.Compiler.return_types`
inspection confirming this. Downstream callers are unaffected:
[`track_path`](@ref)/[`track_bidirectional`](@ref) accept an
already-constructed, concretely-typed `tracker::Tracker` instance, so
*those* two functions remain fully `@inferred`-clean at their own call
sites regardless of which `compile` literal built that tracker
(verified below, and re-verified for `compile = :none` in
`scratch_phase5_check.jl`).

`cfg.path_tracker_precision` is threaded into
`TrackerOptions(min_step_size = Float64(cfg.path_tracker_precision))`.
**This is the most defensible available mapping given
`HomotopyContinuation.jl`'s actual `TrackerOptions` fields
(`automatic_differentiation`, `max_steps`, `max_step_size`,
`max_initial_step_size`, `extended_precision`, `min_step_size`,
`min_rel_step_size`, `terminate_cond`, `parameters`) -- NOT an exact
semantic match** for "requested precision": HC exposes no literal
requested-accuracy/requested-precision knob. Lowering `min_step_size`
lowers the floor the tracker is allowed to shrink its adaptive step to
before giving up (`:terminated_step_size_too_small`), which is the
closest available lever to "how precisely this tracker is willing to
resolve a numerically difficult path" -- it does not *guarantee* any
particular accuracy digit count.
"""
function build_tracker(H_sys::System, x_start::Float64, cfg::HomotopyConfig{T}; compile::Symbol = :all) where {T<:AbstractFloat}
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

Public, generalized replacement for Phase 3's `_track_outward!`: a thin
wrapper that seeds a fresh `budget = Ref(max_steps)` (backed by
`cfg.max_path_steps`, the total number of `ParameterHomotopy` steps this
direction is allowed to take) and delegates to
[`_track_path_segment!`](@ref)'s adaptive bisection.

`F` is the *original* defining system (used only for
[`is_near_singular`](@ref)'s Jacobian evaluation -- `ph`/`tracker` track
the reduced/parameterized system built by the caller); `expected_rank`
for that check is `length(F.expressions)`, matching the convention
already established in `Solver.compute_critical_points`/
`intersect_bounding_object`.
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

Public, generalized replacement for the bidirectional walk previously
hand-rolled inline in Phase 3's `connect_the_dots!`: calls
[`track_path`](@ref) from `(x_mid, y_mid)` toward `x_left` and `x_right`
independently, and returns the already-combined
`reverse(path_left) ++ [[x_mid, real.(y_mid)...]] ++ path_right` (the
concatenation itself, not just the two halves, so callers no longer
need to hand-roll it) plus both landing values, which the caller
resolves against known vertices (curve-specific logic that stays in
`Topology._resolve_endpoint`, not here).
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
