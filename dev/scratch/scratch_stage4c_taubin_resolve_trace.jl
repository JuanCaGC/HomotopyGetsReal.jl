# dev/scratch/scratch_stage4c_taubin_resolve_trace.jl
#
# Stage 4c validation gap-closer (2026-07): the earlier gate-call-count
# instrumentation (scratch_stage4c_taubin_validation.jl) showed 16 of 24
# _deflation_applicable calls returned true during the Taubin heart run, but
# only the 4 surface-level firings had a traceable per-vertex verdict -- the
# other 12 (crit-slice level) get absorbed during clustering/edge-assembly
# before their result can be inspected from the final decomposition output.
# The hard red-flag criteria (Ambiguous, Exhausted, attempts>=15) apply to
# the resolution process itself, not just to vertices that survive to final
# output, so this script instruments resolve_isosingular_dimension directly
# (session-only monkeypatch of the exact committed logic, src/Solver.jl,
# with a log line added at each of the 3 return points), to get a verdict/
# dim/rounds/attempts row for EVERY firing, not just the 4 that persist.
# Investigation script only, not committed as a permanent test.

using HomotopyContinuation
using LinearAlgebra
using GeometryBasics
using GLMakie
using HomotopyGetsReal

@var x y z
f = (x^2 + (1.2 * y)^2 + z^2 - 1)^3 - x^2 * z^3 - 0.1 * (1.2 * y)^2 * z^3
F_heart = System([f], variables = [x, y, z])

cfg = HomotopyConfig{Float64}(
    bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3),
    edge_sample_density = 8, midslice_sample_density = 8,
)

# --- Instrumentation: exact mirror of src/Solver.jl's resolve_isosingular_dimension,
# with a log entry pushed at each of the 3 original return points. Logic is
# unchanged from the committed version -- only the push!(_resolve_log, ...) calls
# and the `level`/`coords0` bookkeeping are new.
const _resolve_log = NamedTuple[]

function HomotopyGetsReal.resolve_isosingular_dimension(
    F_original::System,
    x0::AbstractVector,
    cfg::HomotopyGetsReal.HomotopyConfig{T},
) where {T<:AbstractFloat}
    nv = length(F_original.variables)
    level = nv == 3 ? :surface : (nv == 2 ? :slice : Symbol("nv=$nv"))
    coords0 = round.(real.(x0); digits = 6)

    F_current = F_original
    corank_seq = Int[HomotopyGetsReal.estimate_corank(F_current, x0, cfg; expected_rank = nv)]
    rounds = 0
    last_verify_attempts = missing
    while true
        d = corank_seq[end]
        if d == 0 || HomotopyGetsReal._corank_plateau_hint(corank_seq)
            vr = HomotopyGetsReal.verify_isosingular_dimension(F_current, F_original, x0, d, cfg)
            last_verify_attempts = vr.attempts
            if vr.verdict == HomotopyGetsReal.Verified
                push!(_resolve_log, (level = level, coords = coords0, verdict = HomotopyGetsReal.Resolved,
                    dim = d, rounds = rounds, corank_seq = copy(corank_seq), attempts = vr.attempts))
                return HomotopyGetsReal.ResolveResult{T}(HomotopyGetsReal.Resolved, d, vr.endpoint, corank_seq, rounds, F_current)
            elseif vr.verdict == HomotopyGetsReal.Inconclusive
                push!(_resolve_log, (level = level, coords = coords0, verdict = HomotopyGetsReal.Ambiguous,
                    dim = missing, rounds = rounds, corank_seq = copy(corank_seq), attempts = vr.attempts))
                return HomotopyGetsReal.ResolveResult{T}(HomotopyGetsReal.Ambiguous, nothing, nothing, corank_seq, rounds, F_current)
            end
            # else vr.verdict == NotTerminal: fall through, keep deflating (matches committed logic exactly)
        end
        if rounds >= cfg.max_deflations
            push!(_resolve_log, (level = level, coords = coords0, verdict = HomotopyGetsReal.Exhausted,
                dim = missing, rounds = rounds, corank_seq = copy(corank_seq), attempts = last_verify_attempts))
            return HomotopyGetsReal.ResolveResult{T}(HomotopyGetsReal.Exhausted, nothing, nothing, corank_seq, rounds, F_current)
        end
        F_current, d_new = HomotopyGetsReal.deflate_once(F_current, x0, cfg; expected_rank = nv)
        push!(corank_seq, d_new)
        rounds += 1
    end
end

println("="^70)
println("Step 1: baseline (no deflate)")
println("="^70)
t0 = @elapsed (v0, e0, f0, m0, inc0) = decompose_3d_surface(F_heart, cfg; incidence = true)
println("elapsed: $(round(t0; digits=2))s   total vertices: $(length(v0))")

println()
println("="^70)
println("Step 2: deflate=true, incidence=true -- full diagnostic run (resolve_isosingular_dimension traced)")
println("="^70)
t1 = @elapsed (v1, e1, f1, m1, inc1) = decompose_3d_surface(F_heart, cfg; incidence = true, deflate = true)
println("elapsed: $(round(t1; digits=2))s")
sing1 = filter(v -> v.v_type == Singular, v1)
println("total vertices: $(length(v1))   Singular: $(length(sing1))")

println()
println("="^70)
println("Step 3: surface-level diagnostic (_surface_critical_vertices, deflate=true)")
println("="^70)
fold_anchors_deflated = HomotopyGetsReal._surface_critical_vertices(F_heart, cfg; deflate = true)
surface_singular = filter(v -> v.v_type == Singular, fold_anchors_deflated)
println("surface-level total vertices: $(length(fold_anchors_deflated))   Singular: $(length(surface_singular))")

println()
println("="^70)
println("FULL resolve_isosingular_dimension call log ($(length(_resolve_log)) total calls)")
println("="^70)
println(rpad("level", 8), rpad("coords", 34), rpad("verdict", 12), rpad("dim", 6), rpad("rounds", 8), rpad("attempts", 10), "corank_seq")
for (i, r) in enumerate(_resolve_log)
    println(rpad(string(r.level), 8), rpad(string(r.coords), 34), rpad(string(r.verdict), 12),
        rpad(string(r.dim), 6), rpad(string(r.rounds), 8), rpad(string(r.attempts), 10), string(r.corank_seq))
end

println()
println("="^70)
println("Red-flag audit over ALL $(length(_resolve_log)) real resolve_isosingular_dimension firings")
println("(hard: Ambiguous, Exhausted, or attempts>=15  |  soft: rounds outside [1,2] or attempts outside [1,4])")
println("="^70)
hard_flags = String[]
soft_flags = String[]
for r in _resolve_log
    tag = "level=$(r.level) coords=$(r.coords)"
    if r.verdict == Ambiguous
        push!(hard_flags, "$tag verdict=Ambiguous (attempts=$(r.attempts))")
    elseif r.verdict == Exhausted
        push!(hard_flags, "$tag verdict=Exhausted (rounds=$(r.rounds), last_attempts=$(r.attempts))")
    elseif !ismissing(r.attempts) && r.attempts >= 15
        push!(hard_flags, "$tag verdict=$(r.verdict) attempts=$(r.attempts) >= 15")
    else
        if !(r.rounds in 1:2)
            push!(soft_flags, "$tag rounds=$(r.rounds) outside [1,2]")
        end
        if !ismissing(r.attempts) && !(r.attempts in 1:4)
            push!(soft_flags, "$tag attempts=$(r.attempts) outside [1,4]")
        end
    end
end

println("HARD flags: $(length(hard_flags))")
for f in hard_flags
    println("  FLAG: ", f)
end
println("SOFT flags: $(length(soft_flags))")
for f in soft_flags
    println("  flag: ", f)
end
if isempty(hard_flags) && isempty(soft_flags)
    println("  none -- all $(length(_resolve_log)) firings Resolved, in-range rounds and attempts")
end
