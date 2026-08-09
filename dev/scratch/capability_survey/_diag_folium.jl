using HomotopyContinuation
using HomotopyGetsReal

@var x y
F = System([x^3 + y^3 - 3 * x * y], variables = [x, y])
cfg = HomotopyConfig{Float64}(edge_sample_density = 6)

vertices, edges = decompose_1d_curve(F, cfg)

println("Total vertices: ", length(vertices))
for v in vertices
    println("  id=$(v.id) type=$(v.v_type) coords=$(v.coordinates) metadata=$(v.metadata)")
end
println("Total edges: ", length(edges))
for e in edges
    println("  id=$(e.id) left=$(e.left_vertex_id) right=$(e.right_vertex_id) n_samples=$(length(e.sampled_points)) is_singular=$(e.is_singular)")
end
