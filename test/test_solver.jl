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

end
