using HomotopyContinuation
using HomotopyGetsReal
using HomotopyContinuation: evaluate

@var x y
F = System([x * (x - y) * (x + y)], variables = [x, y])
cfg = HomotopyConfig{Float64}(edge_sample_density = 6)

residual_at(F, point) = abs(only(evaluate(F.expressions, F.variables => point)))

vertices, edges = decompose_1d_curve(F, cfg)

println("Total vertices: ", length(vertices))
referenced_ids = Set{Int}()
for e in edges
    push!(referenced_ids, e.left_vertex_id)
    push!(referenced_ids, e.right_vertex_id)
end
for v in vertices
    r = residual_at(F, real.(v.coordinates))
    referenced = v.id in referenced_ids
    println("  id=$(v.id) type=$(v.v_type) coords=$(real.(v.coordinates)) residual=$(r) referenced_by_edge=$(referenced) metadata=$(v.metadata)")
end
println("Total edges: ", length(edges))
for e in edges
    println("  id=$(e.id) left=$(e.left_vertex_id) right=$(e.right_vertex_id) n_samples=$(length(e.sampled_points))")
end
