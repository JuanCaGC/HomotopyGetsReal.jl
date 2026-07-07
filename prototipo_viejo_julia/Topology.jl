module Topology

using HomotopyContinuation, ..Types, LinearAlgebra

export generate_1d_topology


# Tolerance tuned to 0.15 to bridge the steep Y-drop near vertical tangents
function match_vertex(x_target, y_target, vertices; tol=0.15)
    y_real = real(y_target)
    
    distances = [sqrt((x_target - real(v.coordinates[1]))^2 + 
                      (y_real - real(v.coordinates[2]))^2) for v in vertices]
    
    min_dist, idx = findmin(distances)
    
    if min_dist < tol
        return vertices[idx].id
    else
        return nothing
    end
end

function generate_1d_topology(vertices::Vector{NativeVertex}, f_expr, x_var, y_var; n_samples=80)
    edges = Edge[]
    edge_id_counter = 1 
    
    v_sorted = sort(vertices, by=v -> real(v.coordinates[1]))
    x_coords = [real(v.coordinates[1]) for v in v_sorted]
    
    H_sys = System([f_expr], variables=[y_var], parameters=[x_var])
    ph = ParameterHomotopy(H_sys; start_parameters=[0.0], target_parameters=[0.0])
    path_tracker = Tracker(ph)

    for i in 1:(length(x_coords)-1)
        x_left, x_right = x_coords[i], x_coords[i+1]
        
        # Skip identical or extremely close X coordinates
        if abs(x_left - x_right) < 1e-5
            continue
        end

        x_mid = (x_left + x_right) / 2.0
        sys_mid = System([subs(f_expr, x_var => x_mid)], variables=[y_var])
        res_mid = solve(sys_mid; show_progress=false)
        y_sols = real_solutions(res_mid)

        for current_y in y_sols
            # --- INCREMENTAL LEFT TRACKING ---
            left_samples = range(x_mid, x_left + 1e-3, length=n_samples)
            segment_left = Vector{Float64}[]
            temp_y_left = copy(current_y) # Buffer to hold the previous Y
            prev_x_left = x_mid           # Buffer to hold the previous X
            
            for x_val in left_samples
                start_parameters!(ph, [prev_x_left]) # Start from the LAST known point
                target_parameters!(ph, [x_val])      # Step forward slightly
                res = track(path_tracker, temp_y_left, 1.0, 0.0)
                
                if is_success(res)
                    temp_y_left = solution(res) # Update the Y buffer for the next step
                end
                
                push!(segment_left, [x_val, real(temp_y_left[1])])
                prev_x_left = x_val # Update the X buffer
            end
            
            # --- INCREMENTAL RIGHT TRACKING ---
            right_samples = range(x_mid, x_right - 1e-3, length=n_samples)
            segment_right = Vector{Float64}[]
            temp_y_right = copy(current_y) # Buffer to hold the previous Y
            prev_x_right = x_mid           # Buffer to hold the previous X
            
            for x_val in right_samples
                start_parameters!(ph, [prev_x_right]) # Start from the LAST known point
                target_parameters!(ph, [x_val])       # Step forward slightly
                res = track(path_tracker, temp_y_right, 1.0, 0.0)
                
                if is_success(res)
                    temp_y_right = solution(res) # Update the Y buffer for the next step
                end
                
                push!(segment_right, [x_val, real(temp_y_right[1])])
                prev_x_right = x_val # Update the X buffer
            end

            # Stitch the two halves together
            combined_path = vcat(reverse(segment_left), segment_right)

            # Match the endpoints to our exact NativeVertices
            id_l = match_vertex(x_left, combined_path[1][2], vertices; tol=0.15)
            id_r = match_vertex(x_right, combined_path[end][2], vertices; tol=0.15)

            if id_l !== nothing && id_r !== nothing
                # Pull the exact coordinates from the vertex struct to prevent gaps
                v_l = vertices[findfirst(v->v.id == id_l, vertices)]
                v_r = vertices[findfirst(v->v.id == id_r, vertices)]
                
                # Push the exact critical points to the ends of our tracked path
                pushfirst!(combined_path, [real(v_l.coordinates[1]), real(v_l.coordinates[2])])
                push!(combined_path, [real(v_r.coordinates[1]), real(v_r.coordinates[2])])

                # Save the complete edge
                push!(edges, Edge(edge_id_counter, id_l, combined_path, id_r))
                edge_id_counter += 1
            end
        end
    end
    
    return edges
end

end