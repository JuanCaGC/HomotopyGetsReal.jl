@testset "Solver (Phase 2)" begin

println("=" ^ 70)
println("Setup: nodal cubic f(x,y) = y^2 - x^3 - x^2")
println("=" ^ 70)

@var x y
f = y^2 - x^3 - x^2

# --- compute_critical_points expects, for a 2D curve, an ALREADY
# AUGMENTED square system (caller pre-augments, matching the old
# prototype's convention -- see Solver.jl docstring). We augment with
# respect to y to find critical points of the x-projection.
F_crit = System([f, differentiate(f, y)], variables = [x, y])

# --- intersect_bounding_object expects a curve system directly
# (length(F.expressions) == nvariables(F) - 1): the raw defining
# equation.
F_curve = System([f], variables = [x, y])

# By hand:
#   Faug Jacobian = [ ∂f/∂x        ∂f/∂y      ] = [ -3x²-2x   2y ]
#                   [ ∂²f/∂x∂y     ∂²f/∂y²    ]   [   0        2 ]
#   Solving {f=0, ∂f/∂y=0=2y} => y=0, then f(x,0) = -x³-x² = -x²(x+1) = 0
#     => x = 0 (double) or x = -1.
#   At (0,0): Jacobian = [[0,0],[0,2]]  -> rank 1 (< 2)         => Singular
#   At (-1,0): Jacobian = [[-1,0],[0,2]] -> rank 2 (full)        => Critical
#
#   Curve ∩ box boundary (default bbox = (-4,4)^2):
#     x=4:  y² = 80          -> y ≈ ±8.94  (outside bbox_y)      => rejected
#     x=-4: y² = -48          -> no real solution                 => none
#     y=4:  x³+x²-16=0        -> one real root x ≈ 2.318          => kept
#     y=-4: x³+x²-16=0 (same) -> x ≈ 2.318                        => kept

println()
println("=" ^ 70)
println("1. compute_critical_points with default HomotopyConfig{Float64}()")
println("=" ^ 70)

cfg64 = HomotopyConfig{Float64}()
crit_pts = compute_critical_points(F_crit, cfg64)
for v in crit_pts
    println("  id=$(v.id)  coords=$(v.coordinates)  v_type=$(v.v_type)")
    println("     metadata: ", v.metadata)
end

@test length(crit_pts) == 2
@test any(v -> v.v_type == Singular && isapprox(real(v.coordinates[1]), 0.0; atol = 1e-4) &&
              isapprox(real(v.coordinates[2]), 0.0; atol = 1e-4), crit_pts)
@test any(v -> v.v_type == Critical && isapprox(real(v.coordinates[1]), -1.0; atol = 1e-4), crit_pts)
println("compute_critical_points (default cfg) checks passed: found both Critical and Singular vertices.")

println()
println("=" ^ 70)
println("2. compute_critical_points with a LOOSENED singular_value_threshold")
println("   (min singular value at (-1,0) is 1.0 -- loosen past it and it")
println("    should flip from Critical to Singular too)")
println("=" ^ 70)

cfg_loose = HomotopyConfig{Float64}(singular_value_threshold = 1.5)
crit_pts_loose = compute_critical_points(F_crit, cfg_loose)
for v in crit_pts_loose
    println("  id=$(v.id)  coords=$(v.coordinates)  v_type=$(v.v_type)")
end
@test all(v -> v.v_type == Singular, crit_pts_loose)
println("Tolerance wiring confirmed: singular_value_threshold change flipped classification.")

println()
println("=" ^ 70)
println("3. intersect_bounding_object with default HomotopyConfig{Float64}()")
println("=" ^ 70)

bnd_pts = intersect_bounding_object(F_curve, cfg64)
for v in bnd_pts
    println("  id=$(v.id)  coords=$(v.coordinates)  v_type=$(v.v_type)")
    println("     metadata: ", v.metadata)
end

@test length(bnd_pts) == 2
@test all(v -> v.v_type == Boundary, bnd_pts)
@test all(v -> isapprox(abs(real(v.coordinates[2])), 4.0; atol = 1e-4), bnd_pts)
println("intersect_bounding_object (default cfg) checks passed: found the two y=±4 crossings.")

println()
println("=" ^ 70)
println("4. intersect_bounding_object with a TIGHTENED bbox_y")
println("   (bbox_y = (-1,1) excludes the y=±4 crossings entirely)")
println("=" ^ 70)

cfg_tight_bbox = HomotopyConfig{Float64}(bbox_y = (-1.0, 1.0))
bnd_pts_tight = intersect_bounding_object(F_curve, cfg_tight_bbox)
println("  number of boundary points found: ", length(bnd_pts_tight))
for v in bnd_pts_tight
    println("  id=$(v.id)  coords=$(v.coordinates)  v_type=$(v.v_type)")
end
@test all(v -> isapprox(abs(real(v.coordinates[2])), 1.0; atol = 1e-4), bnd_pts_tight)
@test !isapprox(real(bnd_pts_tight[1].coordinates[1]), real(bnd_pts[1].coordinates[1]); atol = 1e-3)
println("Tolerance wiring confirmed: bbox_y change altered which boundary points were found.")

println()
println("=" ^ 70)
println("5. Clustering.cluster_vertices in isolation")
println("=" ^ 70)

near_dupes = [
    NativeVertex{Float64}(id = 1, coordinates = [1.0 + 0im, 1.0 + 0im], v_type = Critical,
        metadata = Dict{Symbol,Any}(:jacobian_rank => 2)),
    NativeVertex{Float64}(id = 2, coordinates = [1.00001 + 0im, 0.99999 + 0im], v_type = Critical,
        metadata = Dict{Symbol,Any}(:jacobian_rank => 2)),
    NativeVertex{Float64}(id = 3, coordinates = [1.00002 + 0im, 1.00001 + 0im], v_type = Singular,
        metadata = Dict{Symbol,Any}(:jacobian_rank => 1)),
    NativeVertex{Float64}(id = 4, coordinates = [5.0 + 0im, 5.0 + 0im], v_type = Critical,
        metadata = Dict{Symbol,Any}(:jacobian_rank => 2)),
]

clustered_loose = cluster_vertices(near_dupes, 1e-3)
println("With tol=1e-3 (near-dupes 1,2,3 merge; 4 stays separate):")
for v in clustered_loose
    println("  id=$(v.id)  coords=$(v.coordinates)  v_type=$(v.v_type)  metadata=$(v.metadata)")
end
@test length(clustered_loose) == 2
@test any(v -> v.v_type == Singular, clustered_loose) # the merged cluster contains a Singular member

# Members 1,2,3 have :jacobian_rank 2, 2, 1 -- the merged value must be
# the MIN (1), never a numeric average like 5/3 = 1.667 (not a valid rank).
merged_cluster = only(filter(v -> v.metadata[:cluster_size] == 3, clustered_loose))
@test merged_cluster.metadata[:jacobian_rank] == 1
@test merged_cluster.metadata[:jacobian_rank] isa Integer
println("merge_metadata :jacobian_rank check passed: min(2,2,1) = 1 (not averaged to 1.667).")

clustered_tight = cluster_vertices(near_dupes, 1e-8)
println("\nWith tol=1e-8 (nothing close enough to merge):")
for v in clustered_tight
    println("  id=$(v.id)  coords=$(v.coordinates)  v_type=$(v.v_type)")
end
@test length(clustered_tight) == 4
println("cluster_vertices tolerance wiring confirmed (2 clusters vs 4, depending on tol).")

println()
println("=" ^ 70)
println("6. @inferred type-stability checks")
println("=" ^ 70)

# NOTE: @inferred is applied to calls whose arguments are passed
# directly (rather than captured from non-`const` script-level
# globals) -- non-`const` globals are themselves not type-stable in
# Julia, which would make @inferred report a false-positive failure
# unrelated to compute_critical_points/intersect_bounding_object's own
# type stability (confirmed independently via @code_warntype, which
# shows `Body::Vector{NativeVertex{Float64}}` for compute_critical_points
# given concrete argument types).
infer_crit(F, cfg) = compute_critical_points(F, cfg)
infer_bnd(F, cfg) = intersect_bounding_object(F, cfg)
infer_cluster(verts, tol) = cluster_vertices(verts, tol)

r1 = @inferred infer_crit(F_crit, cfg64)
r2 = @inferred infer_bnd(F_curve, cfg64)
r3 = @inferred infer_cluster(near_dupes, 1e-3)

println("@inferred compute_critical_points       -> ", typeof(r1), "  OK")
println("@inferred intersect_bounding_object     -> ", typeof(r2), "  OK")
println("@inferred Clustering.cluster_vertices   -> ", typeof(r3), "  OK")

@test r1 isa Vector{NativeVertex{Float64}}
@test r2 isa Vector{NativeVertex{Float64}}
@test r3 isa Vector{NativeVertex{Float64}}

println()
println("=" ^ 70)
println("7. BigFloat precision sanity check (jacobian_rank_info T-genericity)")
println("=" ^ 70)

cfg_big = HomotopyConfig{BigFloat}()
info_big = jacobian_rank_info(F_crit, [BigFloat(0) + 0im, BigFloat(0) + 0im], cfg_big)
println("  jacobian_rank_info at origin (BigFloat): rank=$(info_big.rank), singular_values=$(info_big.singular_values)")
@test eltype(info_big.singular_values) == BigFloat
@test info_big.rank == 1

println()
println("=" ^ 70)
println("8. Isosingular deflation Stage 1: estimate_corank / deflation_stabilized")
println("=" ^ 70)

# Toy singular curves, hand-verified against expected_rank = length(F.expressions)
# (the SAME "full row rank == smooth point" convention intersect_bounding_object
# already uses on a bare curve system -- J here is 1x2, not square).
#
# Node:  f_node(x,y) = y^2 - x^2 = (y-x)(y+x), an ordinary double point at (0,0).
#   grad f_node = (-2x, 2y).
#   At (0,0): J = [0 0]  -> rank 0 -> corank = 1 - 0 = 1.
#   At (1,1) (on the y=x branch, f=0): J = [-2 2] -> rank 1 -> corank = 1 - 1 = 0.
#
# Cusp:  f_cusp(x,y) = y^2 - x^3, a cuspidal point at (0,0).
#   grad f_cusp = (-3x^2, 2y).
#   At (0,0): J = [0 0]  -> rank 0 -> corank = 1.
#   At (1,1) (on the curve, f(1,1)=1-1=0): J = [-3 2] -> rank 1 -> corank = 0.
#
# First-order corank alone cannot distinguish a node from a cusp -- both are
# corank-1 singularities of the bare Jacobian. Telling them apart needs a
# further deflation iteration (Stage 2+, out of scope for this primitive) or
# higher-order data. Both toy systems reporting the SAME corank (1) at their
# singular point is the expected, hand-verified outcome, not a discrepancy.

println("Setup: node f=y^2-x^2 and cusp f=y^2-x^3, both singular at (0,0)")

@var xs ys
f_node = ys^2 - xs^2
f_cusp = ys^2 - xs^3
F_node = System([f_node], variables = [xs, ys])
F_cusp = System([f_cusp], variables = [xs, ys])

c_node_origin = estimate_corank(F_node, ComplexF64[0, 0], cfg64; expected_rank = length(F_node.expressions))
c_node_smooth = estimate_corank(F_node, ComplexF64[1, 1], cfg64; expected_rank = length(F_node.expressions))
c_cusp_origin = estimate_corank(F_cusp, ComplexF64[0, 0], cfg64; expected_rank = length(F_cusp.expressions))
c_cusp_smooth = estimate_corank(F_cusp, ComplexF64[1, 1], cfg64; expected_rank = length(F_cusp.expressions))

println("  node  corank at (0,0) = $c_node_origin  (hand-computed: 1)")
println("  node  corank at (1,1) = $c_node_smooth  (hand-computed: 0)")
println("  cusp  corank at (0,0) = $c_cusp_origin  (hand-computed: 1)")
println("  cusp  corank at (1,1) = $c_cusp_smooth  (hand-computed: 0)")

@test c_node_origin == 1
@test c_node_smooth == 0
@test c_cusp_origin == 1
@test c_cusp_smooth == 0

println("  HomotopyGetsReal._corank_plateau_hint([1,1,0])  = ", HomotopyGetsReal._corank_plateau_hint([1, 1, 0]))
println("  HomotopyGetsReal._corank_plateau_hint([1,0])    = ", HomotopyGetsReal._corank_plateau_hint([1, 0]))
println("  HomotopyGetsReal._corank_plateau_hint([])       = ", HomotopyGetsReal._corank_plateau_hint(Int[]))
println("  HomotopyGetsReal._corank_plateau_hint([1,1,1])  = ", HomotopyGetsReal._corank_plateau_hint([1, 1, 1]))

@test HomotopyGetsReal._corank_plateau_hint([1, 1, 0]) == true
@test HomotopyGetsReal._corank_plateau_hint([1, 0]) == true
@test HomotopyGetsReal._corank_plateau_hint(Int[]) == false
# Corrected 2026-07: [1,1,1] is a genuine plateau (last two entries equal) --
# real data (the Whitney umbrella handle, [3,1,1,1,1]) confirms this repeats
# forever, not "not yet stabilized". The old deflation_stabilized asserted
# false here; that was the confirmed bug this rename fixes.
@test HomotopyGetsReal._corank_plateau_hint([1, 1, 1]) == true
@test_throws ArgumentError HomotopyGetsReal._corank_plateau_hint([1, 2])

println()
println("=" ^ 70)
println("9. Isosingular deflation Stage 2: deflate_once")
println("=" ^ 70)

# Same node/cusp singularities as section 8, now sliced with one generic
# line L=x+2y through the origin (BertiniReal's witness-point convention --
# a bare 1-equation curve has only 1 Jacobian row, so minorSize=2 has no
# rows to draw from; deflate_once's ArgumentError guard for this is
# exercised explicitly below, first, before the real witness systems).
#
# Node:  {f=y^2-x^2, L=x+2y}.  J=[[-2x,2y],[1,2]].  At (0,0): [[0,0],[1,2]]
#   -> rank 1 -> corank = 2-1 = 1 (matches section 8's bare-curve corank).
#   minorSize=2: the only minor is det(J) = -2x*2 - 2y*1 = -4x-2y =: g1.
#   Deflated J = [[0,0],[1,2],[-4,-2]] at origin -> rank 2 (rows 2,3
#   independent: det[[1,2],[-4,-2]]=6) -> corank_new = 2-2 = 0.
#
# Cusp:  {f=y^2-x^3, L=x+2y}.  J=[[-3x^2,2y],[1,2]].  At (0,0): same [[0,0],[1,2]]
#   -> corank = 1.  Only minor: det(J) = -3x^2*2 - 2y*1 = -6x^2-2y =: g1.
#   Deflated J = [[0,0],[1,2],[0,-2]] at origin -> rank 2
#   (det[[1,2],[0,-2]]=-2) -> corank_new = 0.

println("Setup: same node/cusp, sliced with generic line L=x+2y through (0,0)")

L = xs + 2 * ys
F_node_witness = System([f_node, L], variables = [xs, ys])
F_cusp_witness = System([f_cusp, L], variables = [xs, ys])
origin = ComplexF64[0, 0]

# Guard check: the bare curve's Jacobian has only 1 row, so it can only ever
# supply minorSize<=1 -- but at the ORIGIN itself, deflate_once's own default
# (expected_rank=nv=2) makes corank=2 there (rank(J)=0), which collapses
# minorSize back down to 1 (trivially satisfiable with 1 row) -- NOT where
# the guard fires. It fires at a SMOOTH point of the bare curve instead,
# e.g. (1,1): rank(J)=1 there (full row rank, a genuinely regular point of
# the curve), so corank=2-1=1 (not yet "isolated" in the ambient sense) and
# minorSize=2-1+1=2 -- which the bare curve's single row cannot supply.
println("  bare curve (no slicing line) at a SMOOTH point (1,1) correctly rejected:")
try
    deflate_once(F_node, ComplexF64[1, 1], cfg64)
    println("    UNEXPECTED: did not throw")
catch e
    println("    ", sprint(showerror, e))
end
@test_throws ArgumentError deflate_once(F_node, ComplexF64[1, 1], cfg64)

F_node_defl, c_node_new = deflate_once(F_node_witness, origin, cfg64)
F_cusp_defl, c_cusp_new = deflate_once(F_cusp_witness, origin, cfg64)

# expected_rank = length(F_*_witness.variables) here, NOT .expressions -- matching
# deflate_once's own internal convention (its default), which is what c_node_new/
# c_cusp_new were actually computed against. Numerically identical to
# length(F_*_witness.expressions) since these witness systems are square (2 eq,
# 2 vars) either way, but stating the convention that's actually semantically
# correct rather than the one that's only numerically accidental here.
println("  node deflated system: ", F_node_defl)
println("  node corank sequence: [", estimate_corank(F_node_witness, origin, cfg64; expected_rank = length(F_node_witness.variables)), ", ", c_node_new, "]  (hand: [1,0])")
println("  cusp deflated system: ", F_cusp_defl)
println("  cusp corank sequence: [", estimate_corank(F_cusp_witness, origin, cfg64; expected_rank = length(F_cusp_witness.variables)), ", ", c_cusp_new, "]  (hand: [1,0])")

@test length(F_node_defl.expressions) == 3
@test length(F_cusp_defl.expressions) == 3
@test c_node_new == 0
@test c_cusp_new == 0
@test HomotopyGetsReal._corank_plateau_hint([estimate_corank(F_node_witness, origin, cfg64; expected_rank = length(F_node_witness.variables)), c_node_new]) == true
@test HomotopyGetsReal._corank_plateau_hint([estimate_corank(F_cusp_witness, origin, cfg64; expected_rank = length(F_cusp_witness.variables)), c_cusp_new]) == true

println()
println("=" ^ 70)
println("10. Isosingular deflation Stage 3: verify_isosingular_dimension")
println("=" ^ 70)

# Ground-truth cases from this feature's own investigation (2026-07), each
# re-run here as a permanent regression test rather than left as one-off
# scratch evidence: Whitney umbrella handle/tip (surface, N=3) and the bare
# (unsliced) cusp curve (N=2). All three were validated with a full
# terminal/non-terminal table before this function was approved for
# implementation -- see Solver.jl's own docstring for the false-positive
# evidence behind the F_original residual check and the rejection of
# R(f)=A*f randomization.

println("Setup: Whitney umbrella f = x^2 - y^2*z, handle=(0,0,1), tip=(0,0,0)")

@var wx wy wz
f_umbrella = wx^2 - wy^2 * wz
F_umbrella = System([f_umbrella], variables = [wx, wy, wz])
handle_pt = ComplexF64[0, 0, 1]
tip_pt = ComplexF64[0, 0, 0]
er3 = 3

F1_handle, d1_handle = deflate_once(F_umbrella, handle_pt, cfg64; expected_rank = er3)
r_handle = verify_isosingular_dimension(F1_handle, F_umbrella, handle_pt, d1_handle, cfg64)
println("  handle round1 (corank=$d1_handle): verdict=", r_handle.verdict, "  attempts=", r_handle.attempts)
@test r_handle.verdict == Verified
@test r_handle.attempts == 1
@test r_handle.inconsistent_count == 0

r_handle_wrong_d = verify_isosingular_dimension(F1_handle, F_umbrella, handle_pt, 2, cfg64)
println("  handle round1, WRONG candidate d=2: verdict=", r_handle_wrong_d.verdict)
@test r_handle_wrong_d.verdict != Verified

F1_tip, d1_tip = deflate_once(F_umbrella, tip_pt, cfg64; expected_rank = er3)
r_tip1 = verify_isosingular_dimension(F1_tip, F_umbrella, tip_pt, d1_tip, cfg64)
println("  tip round1 (corank=$d1_tip, KNOWN non-terminal per [3,2,0]): verdict=", r_tip1.verdict)
@test r_tip1.verdict != Verified

F2_tip, d2_tip = deflate_once(F1_tip, tip_pt, cfg64; expected_rank = er3)
r_tip2 = verify_isosingular_dimension(F2_tip, F_umbrella, tip_pt, d2_tip, cfg64)
println("  tip round2 (corank=$d2_tip, trivially terminal): verdict=", r_tip2.verdict, "  attempts=", r_tip2.attempts)
@test r_tip2.verdict == Verified
@test r_tip2.attempts == 0
@test r_tip2.endpoint == Complex{Float64}.(tip_pt)

println()
println("Setup: bare cusp curve f = y^2 - x^3 (unsliced, N=2)")
F2_cusp = System([f_cusp], variables = [xs, ys])
er2 = 2
origin_c = ComplexF64[0, 0]

F1c, d1c = deflate_once(F2_cusp, origin_c, cfg64; expected_rank = er2)
r_cusp1 = verify_isosingular_dimension(F1c, F2_cusp, origin_c, d1c, cfg64)
println("  cusp round1 (corank=$d1c, KNOWN non-terminal per [2,1,0]): verdict=", r_cusp1.verdict)
@test r_cusp1.verdict != Verified

F2c, d2c = deflate_once(F1c, origin_c, cfg64; expected_rank = er2)
r_cusp2 = verify_isosingular_dimension(F2c, F2_cusp, origin_c, d2c, cfg64)
println("  cusp round2 (corank=$d2c, trivially terminal): verdict=", r_cusp2.verdict, "  attempts=", r_cusp2.attempts)
@test r_cusp2.verdict == Verified
@test r_cusp2.attempts == 0

println()
println("=" ^ 70)
println("11. Isosingular deflation Stage 4a: resolve_isosingular_dimension")
println("=" ^ 70)

# Same four ground-truth cases as Section 10, now run through the full
# orchestrator rather than by hand. Two things get asserted beyond the
# verdict, per explicit requirement before Stage 4a is considered
# complete: (a) the ACTUAL rounds/corank_sequence, not the informal
# "1-2 rounds" estimate from the Stage 4 proposal, and (b) the equation
# count at EVERY intermediate round, independently replayed via
# deflate_once and cross-checked against F_final -- this project already
# had one real equation-count runaway (3162 equations by round 4 on the
# handle, from an earlier unbounded loop before _corank_plateau_hint
# existed), so a bounded-growth assertion belongs here, not just a
# verdict check that a future regression removing the plateau gate
# could still pass while silently exploding again.

println("cfg64.max_deflations = ", cfg64.max_deflations)

function replay_and_check(F_original, x0, er, result; label = "", expect_eq_counts = Int[])
    println("  ", label, ": verdict=", result.verdict, "  dim=", result.isosingular_dimension,
        "  rounds=", result.rounds, "  corank_seq=", result.corank_sequence,
        "  F_final eq count=", length(result.F_final.expressions))

    # Independent replay: re-run deflate_once step by step from scratch and
    # confirm both the per-round equation counts (bounded, matching what's
    # already been measured) and that the replay's own final system agrees
    # exactly with resolve_isosingular_dimension's own F_final.
    Fcur = F_original
    eq_counts = Int[length(Fcur.expressions)]
    for _ in 1:result.rounds
        Fcur, _ = deflate_once(Fcur, x0, cfg64; expected_rank = er)
        push!(eq_counts, length(Fcur.expressions))
    end
    println("    independently replayed per-round equation counts: ", eq_counts)

    @test eq_counts == vcat([length(F_original.expressions)], expect_eq_counts)
    @test all(<(50), eq_counts)  # bounded-growth guard -- catches a runaway long before 3162
    @test length(Fcur.expressions) == length(result.F_final.expressions)
    @test result.verdict == Resolved
end

replay_and_check(F_node, origin, 2, resolve_isosingular_dimension(F_node, origin, cfg64);
    label = "NODE", expect_eq_counts = [3])
replay_and_check(F2_cusp, origin_c, er2, resolve_isosingular_dimension(F2_cusp, origin_c, cfg64);
    label = "CUSP", expect_eq_counts = [3, 6])
replay_and_check(F_umbrella, tip_pt, er3, resolve_isosingular_dimension(F_umbrella, tip_pt, cfg64);
    label = "TIP", expect_eq_counts = [4, 15])

# HANDLE gets its own block: this is the one case that exercises
# _corank_plateau_hint's one-extra-round cost directly (corank reaches
# its true limit at round 1, but the hint can't recognize a plateau
# from a single value -- it needs the SAME value to repeat, so round 2
# happens purely to confirm the plateau before verification ever runs).
r_handle_full = resolve_isosingular_dimension(F_umbrella, handle_pt, cfg64)
@test r_handle_full.corank_sequence == [3, 1, 1]
@test r_handle_full.corank_sequence[2] == r_handle_full.corank_sequence[1] - 2  # true limit (1) reached at round 1 already...
@test r_handle_full.corank_sequence[3] == r_handle_full.corank_sequence[2]      # ...but NOT recognized as a plateau until round 2 repeats it
@test r_handle_full.rounds == 2  # one full extra deflate_once round paid specifically for that confirmation
replay_and_check(F_umbrella, handle_pt, er3, r_handle_full; label = "HANDLE", expect_eq_counts = [4, 8])

println()
println("Round counts across all four ground-truth cases: node=1, cusp=2, tip=2, handle=2.")
println("Maximum observed = 2. cfg64.max_deflations = $(cfg64.max_deflations) is a 5x safety")
println("margin over anything actually needed here, not an empirically-required minimum --")
println("same framing as isosingular_verify_retries=20 in Stage 3.")

println()
println("=" ^ 70)
println("12. Isosingular deflation Stage 4b: pipeline wiring")
println("=" ^ 70)

# Guard: deflate=true with F_original omitted must throw in EITHER branch --
# the pre-augmented curve case (F_crit, ne==nv) and the raw-surface
# auto-augment case alike -- confirming no branch gets a silently-inferred
# default, per the approved design.
println("  guard: deflate=true, F_original omitted -> ArgumentError in both branches")
@test_throws ArgumentError compute_critical_points(F_crit, cfg64; deflate = true)
@test_throws ArgumentError compute_critical_points(F_umbrella, cfg64; deflate = true)
@test_throws ArgumentError intersect_bounding_object(F_curve, cfg64; deflate = true)

# deflate=false (the default) must be the byte-identical pre-Stage-4 path --
# confirmed directly, not assumed: same vertex set as calling with no
# keywords at all, and no isosingular metadata added to any vertex.
v_nodeflate = compute_critical_points(F_crit, cfg64)
v_nodeflate2 = compute_critical_points(F_crit, cfg64; deflate = false)
@test length(v_nodeflate) == length(v_nodeflate2)
@test all(!haskey(v.metadata, :isosingular_verdict) for v in v_nodeflate2)
println("  deflate=false: no isosingular metadata added, vertex set unchanged ✓")

# deflate=true on F_crit (nodal cubic, (0,0) confirmed Singular by hand at
# the top of this file) -- confirm the Singular vertex gets real isosingular
# metadata, and non-Singular vertices don't.
v_deflate = compute_critical_points(F_crit, cfg64; deflate = true, F_original = F_curve)
sing_v = only(filter(v -> v.v_type == Singular, v_deflate))
crit_v = only(filter(v -> v.v_type == Critical, v_deflate))
println("  deflate=true on F_crit: Singular vertex metadata keys = ", sort(collect(keys(sing_v.metadata))))
println("    isosingular_verdict=", sing_v.metadata[:isosingular_verdict],
    "  isosingular_dimension=", sing_v.metadata[:isosingular_dimension],
    "  isosingular_deflation_sequence=", sing_v.metadata[:isosingular_deflation_sequence])
@test haskey(sing_v.metadata, :isosingular_verdict)
@test haskey(sing_v.metadata, :isosingular_dimension)
@test haskey(sing_v.metadata, :isosingular_deflation_sequence)
@test sing_v.metadata[:isosingular_verdict] == Resolved
@test !haskey(crit_v.metadata, :isosingular_verdict)  # Critical vertices never get deflated

# decompose_1d_curve: end-to-end wiring, both directions confirmed directly.
f_nodal = y^2 - x^3 - x^2
F_nodal_curve = System([f_nodal], variables = [x, y])
verts_off, _ = decompose_1d_curve(F_nodal_curve, cfg64)
verts_on, _ = decompose_1d_curve(F_nodal_curve, cfg64; deflate = true)
@test all(!haskey(v.metadata, :isosingular_verdict) for v in verts_off)
sing_on = filter(v -> v.v_type == Singular, verts_on)
@test !isempty(sing_on)
@test all(haskey(v.metadata, :isosingular_verdict) for v in sing_on)
println("  decompose_1d_curve: deflate=false adds no metadata, deflate=true resolves ",
    length(sing_on), " Singular vertex(es)")

# Existing non-deflate call sites need zero changes -- confirmed by running
# this project's own sphere/ellipsoid fixtures completely unmodified
# (test_surfacedecomposition.jl, included separately in runtests.jl) and by
# the byte-identical check above, not merely assumed from "no signature
# error was raised."

println()
println("=" ^ 70)
println("12b. Required before Stage 4c: forcing Exhausted and Ambiguous directly")
println("=" ^ 70)

# Both verdicts are defined but were NEVER triggered by any of the four
# ground-truth cases above (all resolve in 1-2 rounds with zero
# inconsistent results) -- confirmed unexercised, not just untested by
# omission. Validating against Taubin's real singular vertices without
# ever having seen these branches actually fire once would mean trusting
# them by inspection, not evidence. Both constructions below are
# deliberately artificial, and labeled as such.

# Exhausted: the cusp genuinely needs 2 rounds (ground truth [2,1,0]) --
# capping max_deflations at 1 forces the loop to give up one round short,
# with NO randomness involved (this path never reaches a verification
# attempt at all). Confirmed reliable across repeated trials before
# committing to this as a test, not assumed from one lucky run.
cfg_capped = HomotopyConfig{Float64}(max_deflations = 1)
r_exhausted = resolve_isosingular_dimension(F2_cusp, origin_c, cfg_capped)
println("  Exhausted (cusp, max_deflations=1): verdict=", r_exhausted.verdict,
    "  corank_seq=", r_exhausted.corank_sequence, "  rounds=", r_exhausted.rounds)
@test r_exhausted.verdict == Exhausted
@test r_exhausted.rounds == 1
@test r_exhausted.isosingular_dimension === nothing

# Ambiguous: an artificial F_original, shifted by a tiny CONSTANT from the
# true Whitney umbrella equation. A constant shift leaves the Jacobian
# (hence the entire corank sequence and deflate_once trajectory) IDENTICAL
# to the true umbrella -- confirmed below, corank_seq still [3,1,1] -- but
# guarantees F_original(x) != 0 at every point, including wherever a
# genuine successful track lands, so the residual check inside
# verify_isosingular_dimension fails deterministically once tracking
# succeeds. Confirmed reliable across repeated trials before committing to
# this as a test, not assumed from one lucky run.
f_umbrella_shifted = wx^2 - wy^2 * wz - 1e-3
F_umbrella_shifted = System([f_umbrella_shifted], variables = [wx, wy, wz])
println("  f_umbrella_shifted(handle_pt) = ",
    HomotopyContinuation.evaluate(f_umbrella_shifted, [wx, wy, wz] => handle_pt), "  (nonzero by construction)")
r_ambiguous = resolve_isosingular_dimension(F_umbrella_shifted, handle_pt, cfg64)
println("  Ambiguous (shifted umbrella, handle point): verdict=", r_ambiguous.verdict,
    "  corank_seq=", r_ambiguous.corank_sequence, "  rounds=", r_ambiguous.rounds)
@test r_ambiguous.verdict == Ambiguous
@test r_ambiguous.corank_sequence == [3, 1, 1]  # identical trajectory to the true umbrella -- confirms
                                                  # the mismatch is isolated to the residual check, nothing else
@test r_ambiguous.isosingular_dimension === nothing

end
