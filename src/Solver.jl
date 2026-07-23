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

Evaluate the Jacobian of `F` at `point` and return its numerical rank and singular values.

Pass a `HomotopyConfig` to supply `jacobian_rank_tol`, the cutoff used when counting
singular values as nonzero. Used by [`compute_critical_points`](@ref) and
[`intersect_bounding_object`](@ref); vertex-type classification also consults
`cfg.singular_value_threshold` separately.
"""
function jacobian_rank_info(F::System, point::AbstractVector, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    # T-precision Jacobian via low-level evaluate(...; bits = precision(T)); svdvals is
    # T-generic once GenericLinearAlgebra is loaded (see module header).
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
    estimate_corank(F::System, point::AbstractVector, cfg::HomotopyConfig{T};
                     expected_rank::Int) where {T<:AbstractFloat}
        -> Int

Isosingular-deflation Stage 1 primitive (2026-07 investigation): the Julia/HC.jl
analogue of BertiniReal's `nullSpace` corank estimate from `isosingular_deflation`
(`src/symbolics/isosingular.cpp:5-134`, upstream `bertini_real` 1.9.0). BertiniReal
gets this number from Bertini1's own `witnessGeneration` under `TrackType: 6` (an
internal numerical-rank/SVD test, not visible in `bertini_real`'s own source); this
reimplements the same quantity directly on top of [`jacobian_rank_info`](@ref),
which already computes a genuinely `T`-generic SVD-based rank.

`corank = expected_rank - rank(J)`. `expected_rank` is a MANDATORY keyword, no
default (2026-07 hardening, following a witness-slice-anchoring investigation):
this function now serves two genuinely different, co-equal conventions with no
single sensible default between them --
[`intersect_bounding_object`](@ref)/`_classify_vertex_type`'s "full ROW rank ==
smooth point" convention (`expected_rank = length(F.expressions)`, meaningful on
a bare, possibly-underdetermined curve/surface equation) versus
[`deflate_once`](@ref)'s "full COLUMN rank == isolated point in ambient space"
convention (`expected_rank = length(F.variables)`, the one confirmed correct
against the Hauenstein-Wampler `D_det` construction). These coincide only when
`F` is square (as `Faug` always is) and diverge for every bare, single-equation
curve/surface -- exactly what deflation is called on. A prior version of this
function defaulted to the row-rank convention; that default was silently wrong
for column-rank callers relying on it, discovered via a verification run whose
own print statement used the default instead of matching `deflate_once`'s
internal choice. Every call site must now state its convention explicitly.

Note the scope boundary: this is a purely first-order (single-Jacobian-evaluation)
quantity. It cannot by itself distinguish singularities that share the same
first-order corank (e.g. a node vs. a cusp both report corank 1 under the
row-rank convention) -- doing so needs either a further deflation iteration
(appending the minor equations and re-evaluating) or higher-order data. This is
expected behavior, not a bug -- see the Stage 1 verification in `test_solver.jl`.
"""
function estimate_corank(
    F::System,
    point::AbstractVector,
    cfg::HomotopyConfig{T};
    expected_rank::Int,
) where {T<:AbstractFloat}
    info = jacobian_rank_info(F, point, cfg)
    return expected_rank - info.rank
end

"""
    deflation_stabilized(corank_sequence::AbstractVector{<:Integer}) -> Bool

Isosingular-deflation Stage 1 primitive: mirrors the two-part check
`isosingular_deflation`'s outer loop performs on BertiniReal's corank sequence
(`src/symbolics/isosingular.cpp:53-97`):

1. The sequence must be nonincreasing -- BertiniReal treats an increase as a
   hard error (`isosingular.cpp:80-85`, "the deflation sequence must be a
   nonincreasing sequence"), not merely "not yet stabilized". Violated here by
   throwing `ArgumentError`, for the same reason: an increase signals a
   numerical failure in the corank estimate itself, not a legitimate
   deflation state.
2. Deflation has succeeded once the corank reaches `0` -- a regular
   (full-rank) point of the current (possibly already-deflated) system, at
   which Newton's method is once again quadratically convergent. This is the
   standard Leykin-Verschelde-Zhao termination criterion; BertiniReal's own
   `success` flag is computed inside Bertini1's closed internals
   (`witnessGeneration`/`isosingular_summary`, not visible in `bertini_real`'s
   own source), so this equivalence is asserted from the published theory, not
   verified line-for-line against Bertini1 itself.

Returns `false` (not yet stabilized) for an empty or nonzero-terminated
sequence.
"""
function deflation_stabilized(corank_sequence::AbstractVector{<:Integer})
    for i in 2:length(corank_sequence)
        corank_sequence[i] <= corank_sequence[i-1] || throw(ArgumentError(
            "corank sequence must be nonincreasing (BertiniReal isosingular.cpp:80-85); " *
            "got $(corank_sequence[i-1]) -> $(corank_sequence[i]) at step $i",
        ))
    end
    return !isempty(corank_sequence) && last(corank_sequence) == 0
end

"""
    deflate_once(F::System, x0::AbstractVector, cfg::HomotopyConfig{T};
                  expected_rank::Int = length(F.variables)) where {T<:AbstractFloat}
        -> (F_new::System, corank_new::Int)

Isosingular deflation Stage 2 primitive (2026-07): a single deflation iteration,
the Julia/HC.jl-native analogue of `createMatlabDeflation`/`deflate_no_subst.m`
(`src/symbolics/isosingular.cpp:315-478`, `matlab_codes/deflate_no_subst.m:110-173`)
-- symbolic minor construction via `HomotopyContinuation.jacobian`/`LinearAlgebra.det`
instead of shelling out to MATLAB/Python.

Appends every `minorSize`-sized minor of `jacobian(F)` that is not the
identically-zero polynomial (`iszero(expand(det(...)))` -- a structural
zero-polynomial test, NOT a numerical test at `x0`: by definition of rank,
*every* `minorSize`-sized minor vanishes numerically at a corank-`c`
singular point when `minorSize = rank(J)+1`, so filtering on value-at-`x0`
would exclude everything; BertiniReal's own filter, `simplify(det(...)) ~= 0`
in Matlab, is the same structural test, just via a different CAS). `minorSize`
is derived from the current corank at `x0` via [`estimate_corank`](@ref)
(`minorSize = expected_rank - corank + 1`, matching `isosingular.cpp:195`'s
formula with no homogenizing-patch correction, since HC.jl works affinely).

Throws `ArgumentError` if `x0` is already full rank at `F` (corank 0 --
nothing to deflate) or if `minorSize` exceeds `F`'s available rows or
columns (`length(F.expressions)` / `length(F.variables)`) -- the latter is
not a bug to work around, it is `deflate_once` telling you `F` is not yet a
valid witness-point system.

**Design implication for future pipeline wiring (Stages 3-4), noted now so
it isn't rediscovered from scratch later**: `deflate_once` assumes `F`
already includes enough generic slicing hyperplanes to be a 0-dimensional
witness-point system at `x0` -- BertiniReal's own convention (see the `P`
patch equation in its Griffis-Duffy fixture, `test/curve/griffisduffy/input`)
and the reason this function's docstring calls `F` a witness system, not
just "the curve/surface equation". `compute_critical_points` and
`decompose_1d_curve`/`decompose_3d_surface` only ever hand this function the
BARE curve/surface equation(s) -- when a `Singular`-classified vertex is
routed into deflation, a new step, not yet implemented, must first construct
a generic slicing hyperplane through that point and augment `F` with it
*before* calling `deflate_once`. Confirmed necessary by hand: a bare single
plane-curve equation has only 1 Jacobian row, so it can never supply the
`minorSize=2` this function's own default `expected_rank = length(F.variables)`
demands away from the singular point itself (at a SMOOTH point, `rank(J)=1`,
so `corank=1` and `minorSize=2`; `test_solver.jl`'s guard-check confirms the
`ArgumentError` fires there) -- note this guard does NOT fire at the singular
point itself with this default, since `rank(J)=0` there collapses
`minorSize` back down to a trivially-satisfiable `1`; that degenerate case
happens to still "deflate" (via two independent 1x1 minors, i.e. the bare
partials) but is not genuine isosingular deflation and should not be relied
upon -- another reason a real slicing hyperplane is needed before wiring
this into the pipeline, not just for the row-count guard.
"""
function deflate_once(
    F::System,
    x0::AbstractVector,
    cfg::HomotopyConfig{T};
    expected_rank::Int = length(F.variables),
) where {T<:AbstractFloat}
    corank = estimate_corank(F, x0, cfg; expected_rank = expected_rank)
    corank > 0 || throw(ArgumentError(
        "deflate_once: system is already full rank (corank 0) at this point -- nothing to deflate.",
    ))

    ne = length(F.expressions)
    nv = length(F.variables)
    minor_size = expected_rank - corank + 1

    minor_size <= ne && minor_size <= nv || throw(ArgumentError(
        "deflate_once: minorSize=$minor_size exceeds available rows ($ne equations) or " *
        "columns ($nv variables) -- F is not yet a 0-dimensional witness-point system; " *
        "see the docstring's note on generic slicing hyperplanes.",
    ))

    J = jacobian(F)
    new_eqs = Expression[]
    for rows in combinations(1:ne, minor_size), cols in combinations(1:nv, minor_size)
        minor = expand(det(J[rows, cols]))
        iszero(minor) && continue
        push!(new_eqs, minor)
    end

    F_new = System(vcat(F.expressions, new_eqs), variables = F.variables)
    corank_new = estimate_corank(F_new, x0, cfg; expected_rank = expected_rank)
    return F_new, corank_new
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

Find and classify critical points of a polynomial system as `NativeVertex` records.

Call with either a square 0-dimensional system (e.g. a pre-augmented curve critical-point
system) or a single equation in three variables (a raw surface; the z-projection critical
system is built internally). Solutions are path-tracked, classified as `Critical` or
`Singular`, and deduplicated with `cfg.vertex_match_tol`.
"""
function compute_critical_points(F::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    # Accept square F (caller-pre-augmented curve case) or 1 eq / 3 vars (auto-augment surface).
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

    # Path-track in Float64; filter by critical_point_tol (not vertex_match_tol); polish when T != Float64.
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

Find curve–bounding-box intersection points and return them as `NativeVertex` records.

`F` must define a plane or space curve (`length(F.expressions) == nvariables(F) - 1`
with two or three variables). Each variable is fixed at its `cfg.bbox_*` bounds in turn;
real solutions inside the box are classified as `Boundary` or `Singular` and
deduplicated with `cfg.vertex_match_tol`.
"""
function intersect_bounding_object(F::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    # Raw surfaces (1 eq, 3 vars) are out of scope: a face intersection is a curve, not isolated points.
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
