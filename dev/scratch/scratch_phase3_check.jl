# scratch_phase3_check.jl
#
# Standalone Phase 3 sanity check (NOT part of test/ yet), continuing
# directly from Phase 2's nodal cubic so the two phases compose
# end-to-end on a hand-verified example.
#
# Test curve: the nodal cubic f(x,y) = y^2 - x^3 - x^2 (same as
# scratch_phase2_check.jl). By hand (see that script's comment block):
#   - compute_critical_points finds (-1,0) [Critical] and (0,0)
#     [Singular, the node].
#   - intersect_bounding_object finds two boundary vertices at
#     (x0, -4) and (x0, +4) for the same x0 (root of x^3+x^2-16=0,
#     x0 ~ 2.318 -- since f depends on y only through y^2, a box
#     crossing at y=+c always mirrors one at y=-c at the same x).
#
# So decompose_1d_curve should produce 4 final vertices and 4 edges:
#   - 2 edges between Critical@(-1,0) and Singular@(0,0) (the curve has
#     two y-branches, y = +-sqrt(x^2*(x+1)), for x in (-1,0)).
#   - 1 edge from Singular@(0,0) to Boundary@(x0,+4).
#   - 1 edge from Singular@(0,0) to Boundary@(x0,-4).
# with NO spurious edge directly between the two Boundary vertices --
# this is exactly the naive-adjacency counterexample from the Phase 3
# architecture discussion: cluster_scalars must collapse the two
# boundary vertices' near-identical x into a single slot before
# interval selection, and connect_the_dots!'s full (x,y) matching must
# tell the two boundary vertices apart from each other despite sharing
# an x-slot.
#
# Run with:
#   julia --project=. scratch_phase3_check.jl

using Test
using HomotopyContinuation
include(joinpath(@__DIR__, "src", "HomotopyGetsReal.jl"))
using .HomotopyGetsReal

println("=" ^ 70)
println("Setup: nodal cubic f(x,y) = y^2 - x^3 - x^2")
println("=" ^ 70)

@var x y
f = y^2 - x^3 - x^2
F_curve = System([f], variables = [x, y])

cfg64 = HomotopyConfig{Float64}()

println()
println("=" ^ 70)
println("1. decompose_1d_curve with default HomotopyConfig{Float64}()")
println("=" ^ 70)

vertices, edges = decompose_1d_curve(F_curve, cfg64)

println("\nVertices (", length(vertices), "):")
for v in vertices
    println("  id=$(v.id)  coords=$(v.coordinates)  v_type=$(v.v_type)")
end

println("\nEdges (", length(edges), "):")
for e in edges
    println(
        "  id=$(e.id)  left=$(e.left_vertex_id)  right=$(e.right_vertex_id)  " *
        "is_singular=$(e.is_singular)  n_points=$(length(e.sampled_points))",
    )
end

@test length(vertices) == 4
@test count(v -> v.v_type == Critical, vertices) == 1
@test count(v -> v.v_type == Singular, vertices) == 1
@test count(v -> v.v_type == Boundary, vertices) == 2
@test any(v -> v.v_type == Critical && isapprox(real(v.coordinates[1]), -1.0; atol = 1e-3) &&
              isapprox(real(v.coordinates[2]), 0.0; atol = 1e-3), vertices)
@test any(v -> v.v_type == Singular && isapprox(real(v.coordinates[1]), 0.0; atol = 1e-3) &&
              isapprox(real(v.coordinates[2]), 0.0; atol = 1e-3), vertices)
@test any(v -> v.v_type == Boundary && isapprox(real(v.coordinates[2]), 4.0; atol = 1e-3), vertices)
@test any(v -> v.v_type == Boundary && isapprox(real(v.coordinates[2]), -4.0; atol = 1e-3), vertices)
println("Vertex checks passed: exactly 1 Critical, 1 Singular, 2 Boundary, at the expected coordinates.")

sing_vertex = only(filter(v -> v.v_type == Singular, vertices))
crit_vertex = only(filter(v -> v.v_type == Critical, vertices))
bnd_plus = only(filter(v -> v.v_type == Boundary && real(v.coordinates[2]) > 0, vertices))
bnd_minus = only(filter(v -> v.v_type == Boundary && real(v.coordinates[2]) < 0, vertices))

function endpoints_of(e)
    Set((e.left_vertex_id, e.right_vertex_id))
end

@test length(edges) == 4
crit_sing_edges = filter(e -> endpoints_of(e) == Set((crit_vertex.id, sing_vertex.id)), edges)
sing_bndplus_edges = filter(e -> endpoints_of(e) == Set((sing_vertex.id, bnd_plus.id)), edges)
sing_bndminus_edges = filter(e -> endpoints_of(e) == Set((sing_vertex.id, bnd_minus.id)), edges)
bnd_bnd_edges = filter(e -> endpoints_of(e) == Set((bnd_plus.id, bnd_minus.id)), edges)

@test length(crit_sing_edges) == 2
@test length(sing_bndplus_edges) == 1
@test length(sing_bndminus_edges) == 1
@test isempty(bnd_bnd_edges)
@test all(e -> e.is_singular, edges) # every edge touches the Singular node at one end
println("Edge topology checks passed: node has degree 4, no spurious Boundary-Boundary edge.")

println()
println("=" ^ 70)
println("2. Clustering.cluster_scalars in isolation")
println("=" ^ 70)

xs = [real(v.coordinates[1]) for v in vertices]
println("  raw x-coordinates: ", xs)
distinct_xs = cluster_scalars(xs, cfg64.vertex_match_tol)
println("  cluster_scalars(xs, vertex_match_tol) -> ", distinct_xs)
@test length(distinct_xs) == 3 # -1, 0, and the (collapsed) boundary x
@test issorted(distinct_xs)
println("cluster_scalars checks passed: the two near-identical boundary x's collapsed into one slot.")

println()
println("=" ^ 70)
println("3. sample_edge point-count / decoupling checks")
println("=" ^ 70)

for e in edges
    @test length(e.sampled_points) == cfg64.edge_sample_density
end
println("Every edge has exactly cfg.edge_sample_density = $(cfg64.edge_sample_density) points.")

cfg_dense = HomotopyConfig{Float64}(edge_sample_density = 200)
_, edges_dense = decompose_1d_curve(F_curve, cfg_dense)
@test length(edges_dense) == length(edges) # connectivity unchanged
@test all(e -> length(e.sampled_points) == 200, edges_dense)
println("Varying edge_sample_density (50 -> 200) changed only the resampled count, not connectivity ($(length(edges_dense)) edges either way).")

cfg_fewsteps = HomotopyConfig{Float64}(max_path_steps = 20)
vertices_fewsteps, edges_fewsteps = decompose_1d_curve(F_curve, cfg_fewsteps)
println(
    "With max_path_steps=20: ", length(vertices_fewsteps), " vertices, ",
    length(edges_fewsteps), " edges (tracking robustness/step-budget only, not final point count).",
)
@test all(e -> length(e.sampled_points) == cfg_fewsteps.edge_sample_density, edges_fewsteps)

println()
println("=" ^ 70)
println("4. HomotopyConfig{BigFloat}() T-genericity check (non-tracking parts)")
println("=" ^ 70)

cfg_big = HomotopyConfig{BigFloat}()
vertices_big, edges_big = decompose_1d_curve(F_curve, cfg_big)
println("Vertices (BigFloat, ", length(vertices_big), "):")
for v in vertices_big
    println("  id=$(v.id)  coords=$(v.coordinates)  v_type=$(v.v_type)")
end
@test eltype(vertices_big) == NativeVertex{BigFloat}
@test eltype(edges_big) == Edge{BigFloat}
@test length(vertices_big) == 4
@test count(v -> v.v_type == Critical, vertices_big) == 1
@test count(v -> v.v_type == Singular, vertices_big) == 1
@test count(v -> v.v_type == Boundary, vertices_big) == 2
@test eltype(vertices_big[1].coordinates) == Complex{BigFloat}
@test eltype(edges_big[1].sampled_points[1]) == BigFloat
println("BigFloat run confirmed T-generic vertex coordinates and classification.")

println()
println("=" ^ 70)
println("5. @inferred type-stability checks")
println("=" ^ 70)

infer_scalars(xs_, tol) = cluster_scalars(xs_, tol)
infer_midslice(F_, xl, xr, cfg_) = compute_midslice(F_, xl, xr, cfg_)
infer_connect(F_, xl, xm, xr, ym, eid, verts, cfg_) = connect_the_dots!(F_, xl, xm, xr, ym, eid, verts, cfg_)
infer_sample(e, cfg_) = sample_edge(e, cfg_)
infer_decompose(F_, cfg_) = decompose_1d_curve(F_, cfg_)

r1 = @inferred infer_scalars(xs, cfg64.vertex_match_tol)
println("@inferred cluster_scalars       -> ", typeof(r1), "  OK")

x_left, x_right = distinct_xs[2], distinct_xs[3] # the (0, boundary_x) interval
r2 = @inferred infer_midslice(F_curve, x_left, x_right, cfg64)
println("@inferred compute_midslice      -> ", typeof(r2), "  OK")

x_mid = (x_left + x_right) / 2
y_mid = r2[1]
verts_for_infer = deepcopy(vertices)
r3 = @inferred infer_connect(F_curve, x_left, x_mid, x_right, y_mid, 999, verts_for_infer, cfg64)
println("@inferred connect_the_dots!     -> ", typeof(r3), "  OK")

r4 = @inferred infer_sample(edges[1], cfg64)
println("@inferred sample_edge           -> ", typeof(r4), "  OK")

r5 = @inferred infer_decompose(F_curve, cfg64)
println("@inferred decompose_1d_curve    -> ", typeof(r5), "  OK")

@test r1 isa Vector{Float64}
@test r2 isa Vector{ComplexF64}
@test r3 isa Edge{Float64}
@test r4 isa Edge{Float64}
@test r5 isa Tuple{Vector{NativeVertex{Float64}},Vector{Edge{Float64}}}

println()
println("=" ^ 70)
println("All Phase 3 sanity checks PASSED.")
println("=" ^ 70)
