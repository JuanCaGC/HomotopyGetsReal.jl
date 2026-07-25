@testset "Isosingular deflation (Stage 5, 2026-07)" begin

# Permanent coverage for decompose_3d_surface(...; deflate=true) and a direct
# unit test of _deflation_applicable -- until now covered only indirectly
# (via compute_critical_points/decompose_1d_curve's own deflate paths) and
# not at all (direct _deflation_applicable calls), respectively. Design
# proposed and approved 2026-07; built on the properties already established
# as normal in dev/scratch/scratch_stage4c_*.jl and
# dev/scratch/scratch_whitney_umbrella_*.jl, not redesigned from scratch.
#
# Fast/slow split follows test_visuals.jl's own precedent (per-testset
# gating within an always-included file), not test_taubin.jl's whole-file
# gate. Only the two cheapest, purely-algebraic testsets (_deflation_applicable
# direct unit test, Whitney umbrella) stay in the fast suite; historical
# curves and the Taubin heart both measured genuinely slow live (see each
# testset's own comment) and are gated.

@testset "_deflation_applicable direct unit test" begin
    # None of the scratch investigation scripts call the real
    # _deflation_applicable directly -- only indirectly via
    # decompose_3d_surface. This closes that gap on cheap fixtures.
    @var xu yu zu
    F_umbrella = System([xu^2 - yu^2 * zu], variables = [xu, yu, zu])
    cfg = HomotopyConfig{Float64}()
    er = 3

    # d==0 short-circuit: needs a genuinely full-rank Jacobian at
    # expected_rank=3, which F_umbrella (1 equation, 3 variables) can never
    # reach on its own -- a single-row Jacobian has rank <= 1, so its corank
    # against expected_rank=3 is always >= 2 (confirmed live: an initial
    # attempt using F_umbrella at (1,1,1) here measured corank=2, not 0, and
    # correctly failed this test until corrected). Use a trivial square,
    # full-rank system instead to exercise the d==0 branch honestly.
    F_square = System([xu - 1, yu - 1, zu - 1], variables = [xu, yu, zu])
    x_smooth = ComplexF64[1, 1, 1]
    d_smooth = estimate_corank(F_square, x_smooth, cfg; expected_rank = er)
    @test d_smooth == 0
    @test HomotopyGetsReal._deflation_applicable(F_square, x_smooth, cfg; expected_rank = er) == true

    # Umbrella tip (0,0,0): bare 1x3 Jacobian is identically zero there, so
    # corank == er == 3. Cross-check the gate's true/false outcome against
    # its own documented formula (minor_size = expected_rank - d + 1, gated
    # on minor_size <= n_eqs && minor_size <= n_vars), not just re-derive it
    # -- the point is a regression check on the real function's behavior at
    # a real fixture, not a tautology.
    x_tip = ComplexF64[0, 0, 0]
    d_tip = estimate_corank(F_umbrella, x_tip, cfg; expected_rank = er)
    @test d_tip == 3
    gate_tip = HomotopyGetsReal._deflation_applicable(F_umbrella, x_tip, cfg; expected_rank = er)
    minor_size_tip = er - d_tip + 1
    @test gate_tip == (minor_size_tip <= length(F_umbrella.expressions) && minor_size_tip <= length(F_umbrella.variables))
    @test gate_tip == true
end

@testset "Whitney umbrella (N=3): corank sequence matches published reference" begin
    # Reuses dev/scratch/scratch_whitney_umbrella_validation.jl /
    # scratch_whitney_umbrella_tip.jl's exact fixture and calls. NOTE: that
    # first scratch script calls a function named `deflation_stabilized`,
    # which no longer exists in src/ (renamed to `_corank_plateau_hint`,
    # src/Solver.jl -- confirmed by reading the current source before
    # writing this test); this test uses the current name.
    @var xw yw zw
    F0 = System([xw^2 - yw^2 * zw], variables = [xw, yw, zw])
    cfg = HomotopyConfig{Float64}()
    er = 3

    function _corank_sequence(x0; max_rounds::Int)
        seq = Int[estimate_corank(F0, x0, cfg; expected_rank = er)]
        Fcur = F0
        round = 0
        while seq[end] > 0 && round < max_rounds
            round += 1
            Fcur, d = deflate_once(Fcur, x0, cfg; expected_rank = er)
            push!(seq, d)
        end
        return seq
    end

    # Tip (0,0,0): an isolated point (published [3,2,0,...]), terminates
    # cleanly and cheaply -- safe to run to actual completion.
    seq_tip = _corank_sequence(ComplexF64[0, 0, 0]; max_rounds = 4)
    @test seq_tip == [3, 2, 0]
    @test HomotopyGetsReal._corank_plateau_hint(seq_tip)

    # Handle (0,0,1): sits on the 1-dim singular locus (published
    # [3,1,1,...]). SAFETY: hard-capped at max_rounds=2, and must STAY
    # capped -- scratch_whitney_umbrella_validation.jl's own header records
    # that round 5 of this exact case needs ~5.3e9 symbolic minors and never
    # terminates (had to be killed). Do not turn this into an open-ended
    # `while d > 0` loop the way the tip case above safely can.
    seq_handle = _corank_sequence(ComplexF64[0, 0, 1]; max_rounds = 2)
    @test seq_handle == [3, 1, 1]
end

if get(ENV, "HOMOTOPYGETSREAL_RUN_SLOW_TESTS", "0") == "1"

using Random

@testset "Historical curves (N=2): deflate=true red-flag baseline" begin
    # Reuses dev/scratch/scratch_stage4c_historical_curves.jl's three
    # fixtures, cfgs, and red-flag logic verbatim.
    #
    # Moved into the slow gate after an initial fast-suite run measured this
    # testset alone at 6m24.8s (curve 3, degree 8, needs up to 5 real
    # deflation rounds with genuine hyperplane-tracking verify attempts per
    # round) -- not "cheap N=2" as originally assumed in the Stage 5 design
    # proposal. Confirmed live, not reasoned from the equation's degree alone.
    @var xv yv
    f1 = (xv^2 + yv^2 - 1)^3 + 27 * xv^2 * yv^2
    f2 = (xv^3 - xv * yv^2 + yv + 1)^2 * (xv^2 + yv^2 - 1) + yv^2 - 5
    f3 = xv^8 - 28 * xv^6 * yv^2 + 70 * xv^4 * yv^4 - 28 * xv^2 * yv^6 + yv^8 + 15 * xv^4 * yv^2 - 15 * xv^2 * yv^4

    F1 = System([f1], variables = [xv, yv])
    F2 = System([f2], variables = [xv, yv])
    F3 = System([f3], variables = [xv, yv])

    cfg_default = HomotopyConfig{Float64}()
    cfg_wide = HomotopyConfig{Float64}(bbox_x = (-8.0, 8.0), bbox_y = (-8.0, 8.0))

    for (label, F, cfg) in (
        ("Curve 1 (default bbox)", F1, cfg_default),
        ("Curve 2 (default bbox)", F2, cfg_default),
        ("Curve 3 (default bbox)", F3, cfg_default),
        ("Curve 3 (bbox=(-8,8))", F3, cfg_wide),
    )
        vertices, edges = decompose_1d_curve(F, cfg; deflate = true)
        sing = filter(v -> v.v_type == Singular, vertices)
        for v in sing
            if haskey(v.metadata, :isosingular_verdict)
                verdict = v.metadata[:isosingular_verdict]
                @test verdict == Resolved
                rounds = length(v.metadata[:isosingular_deflation_sequence]) - 1
                if !(rounds in 1:4)
                    @info "$label id=$(v.id): rounds=$rounds outside established [1,4]" v.metadata[:isosingular_dimension]
                end
            end
        end
    end
end

# Taubin heart fixture, verbatim from test_taubin.jl (same equation, same
# tightened bbox and sample densities).
@var xh yh zh
F_heart = System(
    [(xh^2 + (1.2 * yh)^2 + zh^2 - 1)^3 - xh^2 * zh^3 - 0.1 * (1.2 * yh)^2 * zh^3],
    variables = [xh, yh, zh],
)
cfg_taubin = HomotopyConfig{Float64}(
    bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3),
    edge_sample_density = 8, midslice_sample_density = 8,
)

(v0, _, _, _, _) = decompose_3d_surface(F_heart, cfg_taubin; incidence = true)
(v1, _, _, _, _) = decompose_3d_surface(F_heart, cfg_taubin; incidence = true, deflate = true)
fold_anchors_deflated = HomotopyGetsReal._surface_critical_vertices(F_heart, cfg_taubin; deflate = true)

@testset "Taubin heart: deflate=true preserves slice-level vertex set" begin
    sing0 = filter(v -> v.v_type == Singular, v0)
    sing1 = filter(v -> v.v_type == Singular, v1)
    @test length(sing0) == length(sing1)
    @test sort(collect(v.id for v in sing0)) == sort(collect(v.id for v in sing1))
end

@testset "Taubin heart: isosingular verdicts, zero Ambiguous/Exhausted" begin
    # This is the ONE testset that needs the full per-firing trace (not just
    # the subset that survives to final decomposition output), matching
    # dev/scratch/scratch_stage4c_taubin_resolve_trace.jl's own rationale:
    # candidates absorbed during clustering/edge-assembly are invisible from
    # the final vertex list alone. Getting attempts (not stored in vertex
    # metadata -- only verdict/dimension/corank_sequence are, see
    # src/Solver.jl's compute_critical_points) requires instrumenting
    # resolve_isosingular_dimension itself.
    #
    # 2026-07: uses the REAL, unmodified resolve_isosingular_dimension via
    # its own private _resolve_isosingular_trace ScopedValue hook (see that
    # function's docstring) -- no monkeypatching, no method-table
    # redefinition. Base.ScopedValues.with is task-local and restores the
    # previous value automatically when the block exits (even on error), so
    # this cannot leak into any other test regardless of include order.
    resolve_log = NamedTuple[]
    Base.ScopedValues.with(HomotopyGetsReal._resolve_isosingular_trace => resolve_log) do
        decompose_3d_surface(F_heart, cfg_taubin; incidence = true, deflate = true)
        HomotopyGetsReal._surface_critical_vertices(F_heart, cfg_taubin; deflate = true)
    end

    n_ambiguous = count(r -> r.verdict == Ambiguous, resolve_log)
    n_exhausted = count(r -> r.verdict == Exhausted, resolve_log)
    n_high_attempts = count(r -> !ismissing(r.attempts) && r.attempts >= 15, resolve_log)

    @test !isempty(resolve_log)
    @test n_ambiguous == 0
    @test n_exhausted == 0
    @test n_high_attempts == 0

    # Soft/informational only -- logged, not failed. Range is [1,4], not
    # [1,2]: docs/DESIGN_NOTES.md's Stage 4c entry documents 5 of 22 real
    # Taubin firings needing exactly 4 rounds at the slice-level silhouette
    # points x=(±1,0) (corank_seq=[2,1,1,1,0]) as confirmed-normal, not an
    # anomaly -- flagging those every run would be noise, not signal.
    # attempts==0 occurs legitimately for every d==0 resolution (by
    # construction: verify_isosingular_dimension needs zero hyperplanes when
    # d==0, per docs/DESIGN_NOTES.md's Stage 4c entry) and is excluded from
    # the attempts check entirely, not flagged.
    for r in resolve_log
        if !(r.rounds in 1:4)
            @info "Taubin isosingular soft flag: rounds=$(r.rounds) outside established [1,4]" r.dim r.corank_seq
        end
        if !ismissing(r.attempts) && r.attempts != 0 && !(r.attempts in 1:4)
            @info "Taubin isosingular soft flag: attempts=$(r.attempts) outside established [1,4]" r.dim r.rounds
        end
    end
end

end # HOMOTOPYGETSREAL_RUN_SLOW_TESTS

end
