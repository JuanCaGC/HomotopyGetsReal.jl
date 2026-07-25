# dev/scratch/scratch_whitney_umbrella_validation.jl
#
# Item 3 of the pre-Stage-3 validation request (2026-07): exercise the
# UNMODIFIED deflate_once/estimate_corank (Stage 2, committed a76f38b,
# hardened ec175eb) against the N=3 Whitney umbrella f(x,y,z)=x^2-y^2*z,
# the worked example quoted from Hauenstein-Wampler by the audit session.
# Not a permanent test -- scratch validation only, per explicit instruction.
#
# Published reference sequences (Hauenstein-Wampler, via audit quote):
#   handle point (0,0,z), z != 0  -> {3, 1, 1, ...}
#   tip (0,0,0)                   -> {3, 2, 0, ...}
#
# expected_rank = length(F.variables) = 3 explicitly throughout (per the
# estimate_corank hardening: no default exists anymore, and this is the
# column-rank/ambient-corank convention deflate_once itself uses).

using HomotopyContinuation, LinearAlgebra, HomotopyGetsReal

cfg = HomotopyConfig{Float64}()
@var x y z
f = x^2 - y^2 * z
F0 = System([f], variables = [x, y, z])
er = 3

function run_sequence(label, x0; max_rounds = 6)
    println("="^70)
    println(label, "  x0 = ", x0)
    println("="^70)
    seq = Int[]
    Fcur = F0
    d = estimate_corank(Fcur, x0, cfg; expected_rank = er)
    push!(seq, d)
    println("  round 0: ", length(Fcur.expressions), " eq -> corank = ", d)
    round = 0
    while d > 0 && round < max_rounds
        round += 1
        Fcur, d = deflate_once(Fcur, x0, cfg; expected_rank = er)
        push!(seq, d)
        println("  round $round: ", length(Fcur.expressions), " eq -> corank = ", d)
        println("    system: ", Fcur)
    end
    println("  corank sequence: ", seq, "  stabilized = ", deflation_stabilized(seq))
    return seq
end

seq_handle = run_sequence("HANDLE (0,0,1)", ComplexF64[0, 0, 1])
println()
seq_tip = run_sequence("TIP (0,0,0)", ComplexF64[0, 0, 0])

println()
println("="^70)
println("Comparison against published sequences")
println("="^70)
println("  handle: code = ", seq_handle, "   published = [3,1,1,...]")
println("  tip:    code = ", seq_tip, "   published = [3,2,0,...]")
