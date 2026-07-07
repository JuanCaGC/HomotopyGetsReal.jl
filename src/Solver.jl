# src/Solver.jl
#
# Phase 2: steps 1-2 of the six-step algorithmic framework
# (compute_critical_points, intersect_bounding_object) plus the shared
# Jacobian-rank/singularity-classification utility both steps need.
#
# Precision boundary (see Phase 2 architecture discussion):
#   - Path tracking (`solve`) is genuinely Float64/ComplexF64-only deep
#     inside HomotopyContinuation.jl (its NewtonCache/NewtonResult are
#     hardcoded to ComplexF64) -- there is no way to path-track directly
#     in BigFloat. Every function below therefore path-tracks in
#     Float64 regardless of `T`, exactly like the old prototype.
#   - Symbolic evaluation (`evaluate`, and hence `jacobian`/`F(x)`) IS
#     genuinely T-generic: it bottoms out in SymEngine's `evalf(e, bits)`
#     (arbitrary-precision), but HomotopyContinuation's own convenience
#     wrappers `jacobian(F, x)` / `F(x)` never forward a `bits` keyword
#     and silently default to `bits = 53` (i.e. Float64-equivalent
#     accuracy) no matter what `x`'s element type is. To get genuine
#     T-precision Jacobian evaluation we therefore call the low-level
#     `evaluate(..., F.variables => x; bits = precision(T))` ourselves
#     rather than the `jacobian(F, x)` convenience method.
#   - `LinearAlgebra.svdvals` has no generic fallback for `BigFloat` in
#     Base (verified: throws `MethodError`); `GenericLinearAlgebra.jl`
#     (declared as a package dependency) restores it with matching
#     semantics to the Float64/LAPACK path, so `jacobian_rank_info`
#     below is a single, non-branching, genuinely T-generic
#     implementation.
#   - When `T != Float64`, each accepted Float64 path-tracking solution
#     is refined to genuine T-precision by a small hand-rolled generic
#     Newton corrector (`_newton_polish`), since HomotopyContinuation's
#     own `newton` is likewise hardcoded to ComplexF64 and can't be
#     reused here. The corrector uses `cfg.critical_point_tol` (cast to
#     `T`) as its residual stopping criterion and the same
#     `bits = precision(T)` low-level `evaluate` calls.

"""
    jacobian_rank_info(F::System, point::AbstractVector, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> (rank::Int, singular_values::Vector{T})

Evaluates the Jacobian of `F` at `point` (converted to `Complex{T}`) at
`T`-precision (via low-level `evaluate(...; bits = precision(T))`, see
module docstring above), computes its singular values with
`LinearAlgebra.svdvals` (generic for any `T` once `GenericLinearAlgebra`
is loaded), and reports the numerical rank obtained by counting
singular values strictly greater than `cfg.jacobian_rank_tol`.

This is the shared utility used by both [`compute_critical_points`](@ref)
and [`intersect_bounding_object`](@ref); it does not itself decide
`Singular` vs `Critical`/`Boundary` classification (see
`_classify_vertex_type`, which additionally consults
`cfg.singular_value_threshold` -- a deliberately distinct knob from
`jacobian_rank_tol`, see `HomotopyConfig`'s docstring).
"""
function jacobian_rank_info(F::System, point::AbstractVector, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    x = Complex{T}.(point)
    bits = precision(T)
    Jsym = jacobian(F)
    J = Matrix{Complex{T}}(evaluate(Jsym, F.variables => x; bits = bits))
    svals = Vector{T}(svdvals(J))
    rank = count(>(cfg.jacobian_rank_tol), svals)
    return (rank = rank, singular_values = svals)
end

"""
    _classify_vertex_type(info, cfg::HomotopyConfig{T}, expected_rank::Int, base_type::VertexType) where {T}

Shared classification rule for both Phase 2 entry points: a point is
`Singular` if its numerical Jacobian rank (from [`jacobian_rank_info`](@ref),
driven by `cfg.jacobian_rank_tol`) falls short of `expected_rank` --
the number of defining equations, i.e. full row rank -- **or** if its
smallest singular value is below `cfg.singular_value_threshold` (a
second, independent check using the classification-specific tolerance
documented on `HomotopyConfig`). Otherwise it keeps `base_type`
(`Critical` or `Boundary`, depending on the caller).
"""
function _classify_vertex_type(
    info,
    cfg::HomotopyConfig{T},
    expected_rank::Int,
    base_type::VertexType,
) where {T<:AbstractFloat}
    rank_deficient = info.rank < expected_rank
    near_singular = !isempty(info.singular_values) && minimum(info.singular_values) < cfg.singular_value_threshold
    return (rank_deficient || near_singular) ? Singular : base_type
end

"""
    _newton_polish(F::System, x0::Vector{Complex{T}}, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}

Refines `x0` (a solution of the square system `F`, typically obtained
from Float64 path tracking) to genuine `T`-precision via a hand-rolled
generic Newton iteration, using the low-level `bits = precision(T)`
`evaluate` calls described in the module docstring. A no-op for
`T === Float64` (nothing to gain: `precision(Float64) == 53` is already
what path tracking used). Stops when the residual norm is at most
`cfg.critical_point_tol` or after 50 iterations.
"""
function _newton_polish(F::System, x0::Vector{Complex{T}}, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    T === Float64 && return x0
    bits = precision(T)
    tol = cfg.critical_point_tol
    Jsym = jacobian(F)
    x = copy(x0)
    for _ in 1:50
        Fx = Complex{T}.(evaluate(F.expressions, F.variables => x; bits = bits))
        norm(Fx) <= tol && break
        J = Complex{T}.(evaluate(Jsym, F.variables => x; bits = bits))
        x = x - (J \ Fx)
    end
    return x
end

"""
    compute_critical_points(F::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> Vector{NativeVertex{T}}

Step 1 of the six-step framework.

`F` must be either:
- already square (`length(F.expressions) == nvariables(F)`), i.e. an
  already-augmented, 0-dimensional critical-point system -- this is the
  curve case, matching the old prototype's convention where the caller
  pre-builds e.g. `System([f, differentiate(f, y)], [x, y])` before
  calling in; or
- a single equation in exactly 3 variables (a raw surface), in which
  case the augmented z-projection critical-point system
  `{f, ∂f/∂x, ∂f/∂y}` is built internally (mirrors
  `prototipo_viejo_julia/Solver.jl`'s `find_critical_z_slices`).

Any other shape throws `ArgumentError`.

Pipeline: path-track the (possibly augmented) system in Float64: filter
to solutions whose imaginary part is within `cfg.critical_point_tol`
(NOT `cfg.vertex_match_tol`); optionally Newton-polish to `T`-precision
when `T != Float64`; classify each surviving point via
[`jacobian_rank_info`](@ref) / `_classify_vertex_type` as `Critical` or
`Singular`; deduplicate via `Clustering.cluster_vertices` using
`cfg.vertex_match_tol` (NOT `cfg.critical_point_tol` -- see
`HomotopyConfig`'s docstring for why these two must not be confused).
"""
function compute_critical_points(F::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    nv = length(F.variables)
    ne = length(F.expressions)

    Faug = if ne == nv
        F
    elseif ne == 1 && nv == 3
        x, y, _ = F.variables
        f = F.expressions[1]
        System(vcat(F.expressions, [differentiate(f, x), differentiate(f, y)]), F.variables)
    else
        throw(ArgumentError(
            "compute_critical_points: expected F with length(F.expressions) == nvariables(F) " *
            "(pre-augmented, 0-dimensional system) or a single equation in 3 variables " *
            "(raw surface, auto-augmented internally); got $(ne) equation(s) in $(nv) variable(s).",
        ))
    end

    result = solve(Faug; show_progress = false)
    raw_sols = solutions(result; only_nonsingular = false)

    crit_tol64 = Float64(cfg.critical_point_tol)
    expected_rank = length(Faug.expressions)

    candidates = NativeVertex{T}[]
    next_id = 1
    for s in raw_sols
        maximum(abs, imag.(s)) <= crit_tol64 || continue

        x = _newton_polish(Faug, Complex{T}.(s), cfg)
        info = jacobian_rank_info(Faug, x, cfg)
        v_type = _classify_vertex_type(info, cfg, expected_rank, Critical)
        metadata = Dict{Symbol,Any}(:jacobian_rank => info.rank, :singular_values => info.singular_values)

        push!(candidates, NativeVertex(cfg, next_id, x, v_type; metadata = metadata))
        next_id += 1
    end

    return cluster_vertices(candidates, cfg.vertex_match_tol)
end

"""
    intersect_bounding_object(F::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> Vector{NativeVertex{T}}

Step 2 of the six-step framework.

`F` must define a 1-dimensional variety in `nvariables(F)`-space, i.e.
`length(F.expressions) == nvariables(F) - 1` with `nvariables(F) in (2, 3)`
-- a plane curve (2 variables, matching the old prototype's
`find_boundary_points` exactly) or a space curve (3 variables). A raw
single-equation *surface* (1 equation, 3 variables) is out of scope
here: a surface meets a bounding box face in a curve, not isolated
points, so standard 0-dimensional homotopy continuation cannot return
it as a `Vector{NativeVertex}`; that reduction belongs to per-slice
processing in a later phase.

For each variable and each of its two `cfg.bbox_*` bounds, the variable
is fixed at that bound, the resulting square system in the remaining
variables is path-tracked in Float64, and solutions are kept if they
are real (within `cfg.boundary_point_tol`, an imaginary-part cutoff)
and their remaining coordinates fall within their own bounding
intervals up to `cfg.boundary_point_tol` slack. Kept points are
optionally Newton-polished (see `_newton_polish`) in the reduced
variables, reassembled with the exact fixed coordinate, classified via
the shared Jacobian-rank utility (evaluated against the *original* `F`,
using all `nvariables(F)` variables, so the fixed coordinate's partial
derivative is not lost), tagged `Boundary` (or `Singular`), and
deduplicated via `Clustering.cluster_vertices` using
`cfg.vertex_match_tol`.
"""
function intersect_bounding_object(F::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    nv = length(F.variables)
    ne = length(F.expressions)
    nv in (2, 3) && ne == nv - 1 || throw(ArgumentError(
        "intersect_bounding_object: expected F with length(F.expressions) == nvariables(F) - 1 " *
        "and nvariables(F) in (2, 3) (a plane or space curve); got $(ne) equation(s) in $(nv) variable(s).",
    ))

    bboxes = nv == 2 ? (cfg.bbox_x, cfg.bbox_y) : (cfg.bbox_x, cfg.bbox_y, cfg.bbox_z)
    btol64 = Float64(cfg.boundary_point_tol)
    expected_rank = ne

    candidates = NativeVertex{T}[]
    next_id = 1
    for fixed_idx in 1:nv
        fixed_var = F.variables[fixed_idx]
        remaining_idxs = filter(!=(fixed_idx), 1:nv)
        remaining_vars = F.variables[remaining_idxs]

        for bound_val in bboxes[fixed_idx]
            # Two variants of the same reduced system: `Fsub` keeps the exact
            # T-precision `bound_val` (needed by `_newton_polish` below, which
            # genuinely benefits from it), while `Fsub_solve` always substitutes
            # a Float64 copy of `bound_val` before calling `solve` -- HC's
            # polyhedral start-system construction (`ToricHomotopy`) has no
            # method for `Complex{BigFloat}` coefficients, so path-tracking a
            # system built from literal BigFloat coefficients throws a
            # MethodError. This keeps `solve` itself strictly Float64-only,
            # consistent with the precision boundary documented at the top of
            # this file (T-precision is recovered afterward via polishing).
            exprs_sub = [subs(e, fixed_var => bound_val) for e in F.expressions]
            Fsub = System(exprs_sub, remaining_vars)
            Fsub_solve = if T === Float64
                Fsub
            else
                System([subs(e, fixed_var => Float64(bound_val)) for e in F.expressions], remaining_vars)
            end

            result = solve(Fsub_solve; show_progress = false)
            raw_sols = solutions(result; only_nonsingular = false)

            for s in raw_sols
                maximum(abs, imag.(s)) <= btol64 || continue

                inside = true
                for (k, ridx) in enumerate(remaining_idxs)
                    lo, hi = bboxes[ridx]
                    val = real(s[k])
                    if !(Float64(lo) - btol64 <= val <= Float64(hi) + btol64)
                        inside = false
                        break
                    end
                end
                inside || continue

                x_sub = _newton_polish(Fsub, Complex{T}.(s), cfg)

                full = Vector{Complex{T}}(undef, nv)
                full[fixed_idx] = Complex{T}(bound_val)
                for (k, ridx) in enumerate(remaining_idxs)
                    full[ridx] = x_sub[k]
                end

                info = jacobian_rank_info(F, full, cfg)
                v_type = _classify_vertex_type(info, cfg, expected_rank, Boundary)
                metadata = Dict{Symbol,Any}(
                    :jacobian_rank => info.rank,
                    :singular_values => info.singular_values,
                    :fixed_variable => Symbol(fixed_var),
                    :fixed_value => T(bound_val),
                )

                push!(candidates, NativeVertex(cfg, next_id, full, v_type; metadata = metadata))
                next_id += 1
            end
        end
    end

    return cluster_vertices(candidates, cfg.vertex_match_tol)
end
