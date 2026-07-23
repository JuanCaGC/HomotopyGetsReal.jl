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
    _corank_plateau_hint(corank_sequence::AbstractVector{<:Integer}) -> Bool

Isosingular-deflation Stage 3 primitive (2026-07, renamed from
`deflation_stabilized`): a cheap, NECESSARY-BUT-NOT-SUFFICIENT pre-filter over
a corank sequence -- not an authoritative stabilization test. That role now
belongs to [`verify_isosingular_dimension`](@ref), which actually verifies
smoothness at a candidate dimension via slice-moving continuation, matching
Hauenstein-Wampler's Definition 5.18 ("the isosingular local dimension of `x`
is the LIMIT of the deflation sequence") and the paper's own explicit warning
(Section 6): *"a necessary condition for stabilization is that two
consecutive terms in the deflation sequence must be equal, but this is not
sufficient"* -- their own `f_{k,l}=[x^k,y^l]` counterexample family proves a
sequence can plateau for multiple rounds and still decrease further.

This function checks only the necessary half: Lemma 3.4's nonincreasing
invariant (`N >= d_k >= d_{k+1} >= 0`; thrown as `ArgumentError` on violation,
matching `isosingular.cpp:80-85`'s identical check), plus a plateau signal
(ends in `0`, or its last two entries are equal). A `true` result means "worth
spending a real verification on this round," never "confirmed stabilized."
Intended use: skip [`verify_isosingular_dimension`](@ref) entirely when this
returns `false` -- the corank strictly decreased from the prior round, so by
Lemma 3.4 alone we already know this round isn't terminal, with no need to pay
for a homotopy track to learn what the integer sequence already told us.

Corrected 2026-07: a prior version of this function (`deflation_stabilized`)
returned `false` for `[1,1,1]`, treating "ends in 0" as the only success
signal. That was wrong on two counts -- retracted after finding it
contradicted both Definition 5.18 directly and this project's own
Whitney-umbrella-handle data (a real `deflate_once` trace holding at corank 1
for 4 consecutive rounds, `[3,1,1,1,1]`). `[1,1,1]` now correctly returns
`true`.
"""
function _corank_plateau_hint(corank_sequence::AbstractVector{<:Integer})
    for i in 2:length(corank_sequence)
        corank_sequence[i] <= corank_sequence[i-1] || throw(ArgumentError(
            "corank sequence must be nonincreasing (Hauenstein-Wampler Lemma 3.4; " *
            "isosingular.cpp:80-85); got $(corank_sequence[i-1]) -> $(corank_sequence[i]) at step $i",
        ))
    end
    isempty(corank_sequence) && return false
    last(corank_sequence) == 0 && return true
    return length(corank_sequence) >= 2 && last(corank_sequence) == corank_sequence[end-1]
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

**Correction, 2026-07 (retracts this docstring's own prior claim)**: an
earlier version of this note assumed `F` must already include generic
slicing hyperplanes before `deflate_once` is usable -- a "construct a
witness slice first" pipeline stage was proposed on that basis and then
directly investigated. It was confirmed unnecessary: `deflate_once` called
directly on the BARE, unsliced curve/surface equation (Whitney umbrella
`x^2-y^2*z`, both the handle `(0,0,1)` and tip `(0,0,0)`) reproduces the
Hauenstein-Wampler `D_det` operator's own published deflation sequences
exactly (`{3,1,1,...}` and `{3,2,0,...}`), using this function's own default
`expected_rank = length(F.variables)` with no slice, no witness-point
construction, and no other change. The `minorSize=2`-exceeds-1-row failure
this note used to cite as evidence a slice was required is real (that guard
still fires, correctly, at a SMOOTH point of a bare curve, `rank(J)=1`), but
it does not mean a slice is needed in general -- it means `deflate_once` was
being asked a question (is this an isolated point?) that a smooth point of a
positive-dimensional variety can never answer yes to, independent of
slicing. `compute_critical_points`, `decompose_1d_curve`, and
`decompose_3d_surface` may hand this function the bare curve/surface
equation directly; no intermediate slicing step is needed before Stage 3/4
wiring.
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
    IsosingularVerdict

Three-way outcome of [`verify_isosingular_dimension`](@ref) -- deliberately
not a `Bool`. Collapsing "confirmed non-terminal" and "no clean answer either
way" into a single `false` would hide exactly the distinction that mattered
in this feature's own investigation: a demonstrated 4-of-6 false-positive
rate under naive `R(f)=A*f` randomization, on a case (the cusp curve's round
1) with a mathematically certain non-terminal answer.

- `Verified`: `d` is confirmed as the isosingular local dimension at this
  round -- either trivially (`d == 0`, an already-isolated point) or via a
  slice-moving homotopy track that both completed successfully AND landed on
  a point satisfying the ORIGINAL system to within tolerance.
- `NotTerminal`: every retry attempt cleanly failed to track -- no attempt
  succeeded, and none reported an inconsistent success either. Real,
  repeated evidence this round is not yet stabilized; keep deflating.
- `Inconclusive`: the retry budget was exhausted without a clean verdict.
  Any single attempt reporting tracker-success with a residual outside
  tolerance against the original system forces this outcome immediately,
  regardless of how many clean rejections surround it -- one inconsistent
  result is itself the signal, not noise to retry away.
"""
@enum IsosingularVerdict Verified NotTerminal Inconclusive

"""
    VerifyResult{T<:AbstractFloat}

Return type of [`verify_isosingular_dimension`](@ref).

- `verdict::IsosingularVerdict`
- `endpoint::Union{Nothing,Vector{Complex{T}}}`: the verified point, set only
  when `verdict == Verified`; `nothing` otherwise.
- `attempts::Int`: retry attempts actually used (`0` for the `d == 0`
  trivial case).
- `inconsistent_count::Int`: number of attempts that reported tracker-success
  but failed the original-system residual check -- a count, not just a flag,
  so a caller can see how bad the inconsistency was, not merely that it
  happened.
"""
struct VerifyResult{T<:AbstractFloat}
    verdict::IsosingularVerdict
    endpoint::Union{Nothing,Vector{Complex{T}}}
    attempts::Int
    inconsistent_count::Int
end

"""
    verify_isosingular_dimension(F_current::System, F_original::System, x0::AbstractVector,
                                   d::Int, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> VerifyResult{T}

Isosingular-deflation Stage 3 primitive (2026-07): verifies whether `d` (the
current corank of `F_current` at `x0`) is the TERMINAL isosingular local
dimension, via a Julia/HC.jl-native analogue of Bertini1's `isosingularDimTest`
(`isosingular.c:328-450`) and Hauenstein-Wampler's Algorithm 6.3 -- a
slice-moving parameter-continuation track, not a numeric-sequence pattern
match (see [`_corank_plateau_hint`](@ref) for why that alone is insufficient,
per Definition 5.18 and the paper's own `f_{k,l}` counterexample family).

`F_current` and `F_original` are deliberately separate, both required,
never inferred from one another: `F_current` (the accumulated, possibly
already-deflated system) is what the moving hyperplanes get appended to for
tracking; `F_original` (the bare, undeflated defining equation(s)) is what
the final acceptance residual is checked against. Conflating these into one
argument was exactly the class of bug this signature is designed to make
impossible to write by accident.

Construction (validated against all three ground-truth cases available in
this project -- node/cusp, Whitney umbrella handle and tip, every round
tested, terminal and non-terminal): build `d` generic real hyperplanes
`L_i(x) = sum_j a_ij*(x_j - x0_j)` through `x0` (random real Gaussian
`a_ij`, matching `random_orthogonal_matrix`'s precedent for genericity, not
complex -- `x0` is always real by construction here, coming from a
real-classified vertex), append them to the FULL, UNMODIFIED `F_current` --
no equation discarding, no `R(f)=A*f` randomization. Both alternatives were
tried and rejected: discarding equations before tracking/certifying produces
false positives (demonstrated directly on the Whitney umbrella handle);
`R(f)=A*f` randomization produces a demonstrated 4-of-6 false-positive rate
on the cusp's known-non-terminal round 1 (residuals 0.08-1.79 against the
original system, on tracks that reported `return_code == :success`). Tracks
via a parameter homotopy from `start_parameters = 0` (the hyperplanes' value
at `x0`, by construction) to a fixed nonzero `target_parameters`.

Acceptance requires BOTH conditions, checked explicitly and never silently
combined into one boolean:
1. The tracked path's `return_code == :success`.
2. The tracked endpoint satisfies `F_original` to within
   `cfg.critical_point_tol` -- checked against the ORIGINAL system, never
   `F_current` or the augmented tracking system.

Condition 2 is not optional decoration -- see the false-positive rate cited
above. An implementation that checked only condition 1 would have wrongly
verified the cusp's round 1 as terminal on 4 of 6 random attempts.

Retries: up to `cfg.isosingular_verify_retries` attempts, each with a FRESH
random hyperplane draw (per Hauenstein-Wampler's own stated reliability
recommendation for Algorithm 6.3 -- "perform this test multiple times using
different `L` and `lambda`" -- NOT Bertini1's literal same-draw retry,
`isosingular.c:389-397`, already found unsound in fixed Float64 precision
with no AMP escalation). Every real case checked in this project's own
investigation resolved in exactly 1 attempt; see `cfg.isosingular_verify_retries`'s
own docstring for why the default is nonetheless kept larger, as an
untested-case safety margin rather than an empirically-required minimum.

`d == 0` short-circuits to `Verified` with no tracking at all -- validated
identically on both the Whitney umbrella tip's and the cusp's own round-2
states (`estimate_corank == 0`, confirmed by two independent methods in each
case, `F_original` residual exactly `0.0` in both).
"""
function verify_isosingular_dimension(
    F_current::System,
    F_original::System,
    x0::AbstractVector,
    d::Int,
    cfg::HomotopyConfig{T},
) where {T<:AbstractFloat}
    if d == 0
        return VerifyResult{T}(Verified, Complex{T}.(x0), 0, 0)
    end

    nv = length(F_current.variables)
    orig_tol = Float64(cfg.critical_point_tol)
    x0_c = ComplexF64.(x0)
    inconsistent_count = 0

    for attempt in 1:cfg.isosingular_verify_retries
        A = randn(d, nv)
        s = Variable.(:__verify_dimension_s, 1:d)
        hyps = [
            sum(A[i, j] * (F_current.variables[j] - x0_c[j]) for j in 1:nv) - s[i]
            for i in 1:d
        ]
        Faug = System(vcat(F_current.expressions, hyps), variables = F_current.variables, parameters = s)

        start_par = zeros(ComplexF64, d)
        target_par = fill(ComplexF64(1.0), d)

        result = solve(
            Faug, [x0_c];
            start_parameters = start_par, target_parameters = target_par, show_progress = false,
        )
        pr = only(path_results(result))

        if pr.return_code == :success
            resid = maximum(abs, HomotopyContinuation.evaluate(
                F_original.expressions, F_original.variables => pr.solution,
            ))
            if resid <= orig_tol
                return VerifyResult{T}(Verified, Complex{T}.(pr.solution), attempt, inconsistent_count)
            else
                inconsistent_count += 1
                return VerifyResult{T}(Inconclusive, nothing, attempt, inconsistent_count)
            end
        end
        # clean rejection (return_code != :success) -- retry with a fresh draw
    end

    return VerifyResult{T}(NotTerminal, nothing, cfg.isosingular_verify_retries, inconsistent_count)
end

"""
    ResolveVerdict

Outcome of [`resolve_isosingular_dimension`](@ref) -- the multi-round
analogue of [`IsosingularVerdict`](@ref), which only ever describes a
single round.

Deliberately uses `Ambiguous`, not `Inconclusive`, for its middle value --
`Solver.jl` is a single flat namespace (everything here is `include`d into
one module, no submodule separation), and a first version of this enum
reused `Inconclusive` verbatim from [`IsosingularVerdict`](@ref). That
silently shadowed the earlier binding at module scope: every bare
`Inconclusive` reference in the file, including inside the
already-shipped [`verify_isosingular_dimension`](@ref), started
resolving to *this* enum's value instead, and the `vr.verdict ==
Inconclusive` comparison inside this file's own orchestration loop
below would have silently always been `false` (comparing across two
different enum types, which Julia allows and just returns `false` for --
no error, no warning). Caught before ever being tested, not after --
renamed instead of qualified, so the possibility can't recur.

- `Resolved`: the isosingular local dimension was found -- either
  trivially (corank reached `0`) or via a [`verify_isosingular_dimension`](@ref)
  call that returned `Verified`.
- `Ambiguous`: some round's verification returned `Inconclusive` --
  propagated up immediately, not retried at the orchestration level (a
  single inconsistent round is a reason to stop and look, not to keep
  deflating past it).
- `Exhausted`: `cfg.max_deflations` rounds were spent without ever
  reaching `Resolved` or `Ambiguous` -- distinct from `Ambiguous` on
  purpose: this means the corank sequence never stabilized long enough
  to even trigger a verification attempt, not that a verification
  attempt gave a mixed signal.
"""
@enum ResolveVerdict Resolved Ambiguous Exhausted

"""
    ResolveResult{T<:AbstractFloat}

Return type of [`resolve_isosingular_dimension`](@ref).

- `verdict::ResolveVerdict`
- `isosingular_dimension::Union{Nothing,Int}`: set only when `verdict == Resolved`.
- `endpoint::Union{Nothing,Vector{Complex{T}}}`: set only when `verdict == Resolved`.
- `corank_sequence::Vector{Int}`: the full `d_0, d_1, ..., d_rounds` trace --
  kept for diagnostics regardless of verdict, not just on success.
- `rounds::Int`: number of `deflate_once` calls actually performed.
- `F_final::System`: the system in force when the loop stopped (`F_original`
  itself if `rounds == 0`). Not consumed anywhere yet -- Stage 4a is
  diagnostic-only -- but kept on the result rather than discarded, since
  a future precision-improvement stage would need exactly this.
"""
struct ResolveResult{T<:AbstractFloat}
    verdict::ResolveVerdict
    isosingular_dimension::Union{Nothing,Int}
    endpoint::Union{Nothing,Vector{Complex{T}}}
    corank_sequence::Vector{Int}
    rounds::Int
    F_final::System
end

"""
    resolve_isosingular_dimension(F_original::System, x0::AbstractVector, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> ResolveResult{T}

Isosingular-deflation Stage 4a primitive (2026-07): the orchestrator Stage
3 deliberately didn't build -- ties [`estimate_corank`](@ref),
[`_corank_plateau_hint`](@ref), [`deflate_once`](@ref), and
[`verify_isosingular_dimension`](@ref) together into a single call that
resolves a point's isosingular local dimension (Hauenstein-Wampler
Definition 5.18) from a corank sequence's first entry through to a
verified terminal round.

Loop, once per round: if the current corank is `0`, or
`_corank_plateau_hint` says the sequence might have stabilized, call
`verify_isosingular_dimension` on the CURRENT (possibly already-deflated)
system against the ORIGINAL `F_original` -- `Verified` returns `Resolved`
immediately; `Inconclusive` returns `Ambiguous` immediately, not
retried at this level. Otherwise (`_corank_plateau_hint` says "definitely
not yet", or verification returned `NotTerminal`), call `deflate_once`
and continue, up to `cfg.max_deflations` rounds before returning
`Exhausted`.

**The one-extra-round cost of a real (nonzero) plateau, confirmed
directly, not assumed**: `_corank_plateau_hint` needs to see a repeated
value to recognize a plateau, so a genuine nonzero-limit case always
costs one deflation round beyond where the corank first reaches its
true limit. Concretely, on the Whitney umbrella handle: `corank_sequence
== [3, 1, 1]`, `rounds == 2` -- corank already reached its terminal value
`1` after round 1 (`F1`, 4 equations), but the loop cannot know that yet
(`[3,1]` has no repeat) and must deflate once more (`F2`, 8 equations)
before `_corank_plateau_hint([3,1,1])` finally sees the repeat and
triggers verification. This is not a cost that appears for a `d==0`
resolution (node/cusp/tip below) -- corank hitting `0` always triggers
verification on the SAME round it's first observed, no repeat needed.

**Round counts, measured directly against every ground-truth case
available in this project, not estimated**:

| case | corank_sequence | rounds | verdict |
|---|---|---|---|
| node (bare curve) | `[2, 0]` | 1 | Resolved, dim 0 |
| cusp (bare curve) | `[2, 1, 0]` | 2 | Resolved, dim 0 |
| Whitney umbrella tip | `[3, 2, 0]` | 2 | Resolved, dim 0 |
| Whitney umbrella handle | `[3, 1, 1]` | 2 | Resolved, dim 1 |

Maximum observed: 2 rounds. `cfg.max_deflations = 10` is therefore, on
the identical footing as `cfg.isosingular_verify_retries`, an explicit
safety margin for cases this project hasn't exercised yet (deeper
singularities, real surfaces under Stage 4c) -- not a value any
ground-truth case has ever needed more than a fifth of.
"""
function resolve_isosingular_dimension(
    F_original::System,
    x0::AbstractVector,
    cfg::HomotopyConfig{T},
) where {T<:AbstractFloat}
    nv = length(F_original.variables)
    F_current = F_original
    corank_seq = Int[estimate_corank(F_current, x0, cfg; expected_rank = nv)]
    rounds = 0

    while true
        d = corank_seq[end]
        if d == 0 || _corank_plateau_hint(corank_seq)
            vr = verify_isosingular_dimension(F_current, F_original, x0, d, cfg)
            if vr.verdict == Verified
                return ResolveResult{T}(Resolved, d, vr.endpoint, corank_seq, rounds, F_current)
            elseif vr.verdict == Inconclusive
                return ResolveResult{T}(Ambiguous, nothing, nothing, corank_seq, rounds, F_current)
            end
            # NotTerminal (only reachable when d > 0, since d == 0 always verifies
            # trivially) -- fall through and deflate further.
        end
        rounds >= cfg.max_deflations && return ResolveResult{T}(Exhausted, nothing, nothing, corank_seq, rounds, F_current)
        F_current, d_new = deflate_once(F_current, x0, cfg; expected_rank = nv)
        push!(corank_seq, d_new)
        rounds += 1
    end
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
