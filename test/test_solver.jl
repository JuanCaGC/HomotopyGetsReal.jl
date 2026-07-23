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

c_node_origin = estimate_corank(F_node, ComplexF64[0, 0], cfg64)
c_node_smooth = estimate_corank(F_node, ComplexF64[1, 1], cfg64)
c_cusp_origin = estimate_corank(F_cusp, ComplexF64[0, 0], cfg64)
c_cusp_smooth = estimate_corank(F_cusp, ComplexF64[1, 1], cfg64)

println("  node  corank at (0,0) = $c_node_origin  (hand-computed: 1)")
println("  node  corank at (1,1) = $c_node_smooth  (hand-computed: 0)")
println("  cusp  corank at (0,0) = $c_cusp_origin  (hand-computed: 1)")
println("  cusp  corank at (1,1) = $c_cusp_smooth  (hand-computed: 0)")

@test c_node_origin == 1
@test c_node_smooth == 0
@test c_cusp_origin == 1
@test c_cusp_smooth == 0

println("  deflation_stabilized([1,1,0])  = ", deflation_stabilized([1, 1, 0]))
println("  deflation_stabilized([1,0])    = ", deflation_stabilized([1, 0]))
println("  deflation_stabilized([])       = ", deflation_stabilized(Int[]))
println("  deflation_stabilized([1,1,1])  = ", deflation_stabilized([1, 1, 1]))

@test deflation_stabilized([1, 1, 0]) == true
@test deflation_stabilized([1, 0]) == true
@test deflation_stabilized(Int[]) == false
@test deflation_stabilized([1, 1, 1]) == false
@test_throws ArgumentError deflation_stabilized([1, 2])

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

println("  node deflated system: ", F_node_defl)
println("  node corank sequence: [", estimate_corank(F_node_witness, origin, cfg64), ", ", c_node_new, "]  (hand: [1,0])")
println("  cusp deflated system: ", F_cusp_defl)
println("  cusp corank sequence: [", estimate_corank(F_cusp_witness, origin, cfg64), ", ", c_cusp_new, "]  (hand: [1,0])")

@test length(F_node_defl.expressions) == 3
@test length(F_cusp_defl.expressions) == 3
@test c_node_new == 0
@test c_cusp_new == 0
@test deflation_stabilized([estimate_corank(F_node_witness, origin, cfg64), c_node_new]) == true
@test deflation_stabilized([estimate_corank(F_cusp_witness, origin, cfg64), c_cusp_new]) == true

end
