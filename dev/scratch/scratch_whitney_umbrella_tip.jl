# dev/scratch/scratch_whitney_umbrella_tip.jl
#
# Item 5 of the pre-Stage-3 validation (2026-07): the TIP (0,0,0) of the
# Whitney umbrella f=x^2-y^2*z, run SEPARATELY from the handle case (whose
# run had to be killed -- round 5 there needed ~5.3B symbolic minors and
# never terminates). Published reference: {3, 2, 0, ...} -- an isolated
# point (unlike the handle, which sits on the 1-dim singular locus), so
# this SHOULD reach corank 0 and stop after round 2, not runaway.
# Capped at max_rounds=4 as a safety bound in case that expectation is wrong.

using HomotopyContinuation, LinearAlgebra, HomotopyGetsReal

cfg = HomotopyConfig{Float64}()
@var x y z
f = x^2 - y^2 * z
F0 = System([f], variables = [x, y, z])
er = 3

println("="^70)
println("TIP (0,0,0)")
println("="^70)
x0 = ComplexF64[0, 0, 0]
seq = Int[]
Fcur = F0
d = estimate_corank(Fcur, x0, cfg; expected_rank = er)
push!(seq, d)
println("  round 0: ", length(Fcur.expressions), " eq -> corank = ", d)
round = 0
max_rounds = 4
while d > 0 && round < max_rounds
    global round, Fcur, d
    round += 1
    Fcur, d = deflate_once(Fcur, x0, cfg; expected_rank = er)
    push!(seq, d)
    println("  round $round: ", length(Fcur.expressions), " eq -> corank = ", d)
end
println("  tip corank sequence: ", seq, "   published = [3,2,0,...]")
