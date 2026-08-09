using HomotopyContinuation
using HomotopyGetsReal
using HomotopyContinuation: evaluate

@var x y z
F = System([x^4 + y^4 + z^4 - 1], variables = [x, y, z])
cfg = HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8)

residual_at(F, point) = abs(only(evaluate(F.expressions, F.variables => point)))

vertices, edges, faces, mesh = decompose_3d_surface(F, cfg)

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
println("Total faces: ", length(faces))
println("mesh: $(length(GeometryBasics.coordinates(mesh))) verts, $(length(GeometryBasics.faces(mesh))) tris")
