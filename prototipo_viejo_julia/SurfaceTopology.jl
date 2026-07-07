module SurfaceTopology

using HomotopyContinuation, LinearAlgebra, GeometryBasics
using ..Types
using ..Solver
using ..Topology

export generate_3d_topology, weld_faces_to_mesh

# Helper to convert raw points into NativeVertices for the 1D engine
function build_2d_vertices(f_2d, x_var, y_var)
    # 1. Find Critical Points of the 2D slice
    crit_sys = System([f_2d, differentiate(f_2d, y_var)], variables=[x_var, y_var])
    raw_crits = find_critical_points(crit_sys)
    
    # 2. Add some boundary points (Bounding box: x in [-4, 4], y in [-4, 4])
    raw_bounds = find_boundary_points(f_2d, [x_var, y_var], [-4.0, 4.0], [-4.0, 4.0])
    
    vertices = NativeVertex[]
    v_id = 1
    
    for pt in raw_crits
        push!(vertices, NativeVertex(v_id, [ComplexF64(pt[1]), ComplexF64(pt[2])], Critical))
        v_id += 1
    end
    for pt in raw_bounds
        push!(vertices, NativeVertex(v_id, [ComplexF64(pt[1]), ComplexF64(pt[2])], Boundary))
        v_id += 1
    end
    
    return vertices
end

# The Parameter Homotopy Sweeper - Final "Welded" Version
function track_face(f_expr, vars, edge_1d, z_start, z_target, face_id; steps=300)
    x_var, y_var, z_var = vars
    @var x0 y0
    
    radial_eq = x_var * y0 - y_var * x0
    H_sys = System([f_expr, radial_eq], variables=[x_var, y_var], parameters=[z_var, x0, y0])
    
    ph = ParameterHomotopy(H_sys; start_parameters=[z_start, 0.0, 0.0], target_parameters=[z_target, 0.0, 0.0])
    tracker = Tracker(ph)
    
    # Safety margin for the numerical solver
    ε = 1e-4
    direction = sign(z_target - z_start)
    safe_z_target = z_target - (direction * ε)
    
    z_vals = range(z_start, safe_z_target, length=steps)
    mesh_vertices = Vector{Float64}[]
    
    active_vertices = trues(length(edge_1d.path))
    current_coords = [Float64[pt[1], pt[2], z_start] for pt in edge_1d.path]
    z_prev = z_start
    
    # --- MAIN TRACKING LOOP ---
    for z_curr in z_vals
        for (i, pt) in enumerate(edge_1d.path)
            x_anchor, y_anchor = pt[1], pt[2]
            x_old, y_old, z_old = current_coords[i][1], current_coords[i][2], current_coords[i][3]
            
            if !active_vertices[i]
                push!(mesh_vertices, [x_old, y_old, z_old])
                continue
            end
            
            start_parameters!(ph, [z_prev, x_anchor, y_anchor])
            target_parameters!(ph, [z_curr, x_anchor, y_anchor])
            
            res = track(tracker, [ComplexF64(x_old), ComplexF64(y_old)], 1.0, 0.0)
            
            if is_success(res)
                sol = solution(res)
                x_new, y_new = real(sol[1]), real(sol[2])
                
                # Bounding Box Clamp (Prevents holes and wild divergence)
                if abs(x_new) > 2.0 || abs(y_new) > 2.0
                    active_vertices[i] = false
                    push!(mesh_vertices, [x_old, y_old, z_old])
                else
                    current_coords[i] = [x_new, y_new, z_curr]
                    push!(mesh_vertices, [x_new, y_new, z_curr])
                end
            else
                active_vertices[i] = false
                push!(mesh_vertices, [x_old, y_old, z_old])
            end
        end
        z_prev = z_curr 
    end
    
    # --- THE WELDING ROW (The Gap Closer) ---
    # We add one final row that bridges the 1e-4 gap to the actual z_target
    for (i, pt) in enumerate(edge_1d.path)
        x_last, y_last, z_last = current_coords[i]
        
        if active_vertices[i]
            # Stretch the successful points to the absolute boundary
            push!(mesh_vertices, [x_last, y_last, z_target])
        else
            # Keep dead points frozen to prevent horns
            push!(mesh_vertices, [x_last, y_last, z_last])
        end
    end
    
    mesh_matrix = reduce(vcat, transpose.(mesh_vertices))
    
    # Important: total rows is now steps + 1
    dim_matrix = [(steps + 1) length(edge_1d.path)]
    
    return Face(face_id, z_start, [edge_1d.id], mesh_matrix, dim_matrix)
end

function generate_3d_topology(f_expr, vars; z_range=(-3.5, 3.5))
    x, y, z = vars
    faces = Face[]
    face_id_counter = 1
    
    println("--- Starting 3D Decomposition ---")
    
    # 1. Find Critical Z-Slices (The Slabs)
    z_crits = find_critical_z_slices(f_expr, vars)
    z_bounds = sort(unique([z_range[1]; z_crits; z_range[2]]))
    
    println("Z-Milestones identified at: ", round.(z_bounds, digits=3))
    
    # 2. Loop through the Slabs
    for i in 1:(length(z_bounds)-1)
        z_bottom = z_bounds[i]
        z_top = z_bounds[i+1]
        z_mid = z_bottom + 0.4137 * (z_top - z_bottom)
        
        println("\nProcessing Slab $(i): Z in [$(round(z_bottom, digits=2)), $(round(z_top, digits=2))]")
        println("  -> Slicing at mid-point: z = $(round(z_mid, digits=2))")
        
        # 3. Slicing & 1D Extraction
        f_2d = slice_and_extract_2d(f_expr, vars, z_mid)
        vertices_2d = build_2d_vertices(f_2d, x, y)
        
        println("  -> Handing 2D slice to 1D Topology Engine...")
        edges_2d = generate_1d_topology(vertices_2d, f_2d, x, y)
        println("  -> 1D Engine found $(length(edges_2d)) distinct curves at this slice.")
        
        # 4. Sweep edges to build 3D Faces
        for edge in edges_2d
            # Sweep UP
            face_up = track_face(f_expr, vars, edge, z_mid, z_top, face_id_counter)
            push!(faces, face_up)
            face_id_counter += 1
            
            # Sweep DOWN
            face_down = track_face(f_expr, vars, edge, z_mid, z_bottom, face_id_counter)
            push!(faces, face_down)
            face_id_counter += 1
        end
    end
    
    println("\n--- 3D Decomposition Complete! ---")
    println("Total 3D Faces generated: ", length(faces))
    return faces
end

"""
    weld_faces_to_mesh(faces; tol=1e-5)

Takes a list of Face objects, merges coincident vertices within `tol`, 
and returns a single watertight HomogeneousMesh.
"""
function weld_faces_to_mesh(faces; tol=1e-5)
    all_vertices = Point3f[]
    all_triangles = TriangleFace{Int}[]
    
    # This dictionary maps a unique 3D coordinate to a single vertex index
    # We round to the tolerance to ensure nearly identical points merge
    vertex_map = Dict{Vector{Float32}, Int}()
    
    function get_unique_idx(pt)
        # 1. Round to handle floating point noise
        # 2. Convert to Float32 for the GPU-friendly Point3f
        pt_f32 = Float32.(pt)
        rounded = round.(pt_f32, digits=5) 
        
        if haskey(vertex_map, rounded)
            return vertex_map[rounded]
        else
            push!(all_vertices, Point3f(pt_f32...))
            new_idx = length(all_vertices)
            vertex_map[rounded] = new_idx
            return new_idx
        end
    end

    for f in faces
        rows, cols = Int(f.dim_matrix[1]), Int(f.dim_matrix[2])
        m = f.mesh_matrix
        
        # Step A: Map all points in this face to global unique IDs
        grid_ids = Matrix{Int}(undef, rows, cols)
        for r in 1:rows
            for c in 1:cols
                # Calculate index in the flat mesh_matrix
                idx_in_matrix = (r - 1) * cols + c
                grid_ids[r, c] = get_unique_idx(m[idx_in_matrix, :])
            end
        end
        
        # Step B: Create triangles (triangulate the grid)
        for r in 1:(rows - 1)
            for c in 1:(cols - 1)
                v1 = grid_ids[r,   c]
                v2 = grid_ids[r+1, c]
                v3 = grid_ids[r+1, c+1]
                v4 = grid_ids[r,   c+1]
                
                # We check if the triangle has "area" (v1 != v2 != v3)
                # This prevents glitches at the sharp tips of the heart
                if v1 != v2 && v2 != v3 && v3 != v1
                    push!(all_triangles, TriangleFace(v1, v2, v3))
                end
                if v1 != v3 && v3 != v4 && v4 != v1
                    push!(all_triangles, TriangleFace(v1, v3, v4))
                end
            end
        end
    end
    
    return Mesh(all_vertices, all_triangles)
end

end