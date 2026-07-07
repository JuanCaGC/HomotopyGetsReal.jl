module Solver
    using HomotopyContinuation, LinearAlgebra, ..Types

    export find_critical_points, find_boundary_points, check_singularity
    export find_critical_z_slices, slice_and_extract_2d

    # Existing critical point finder
    function find_critical_points(poly_system)
        result = solve(poly_system; seed = 0x12345, compile = true, show_progress = false)
        raw_real = real_solutions(result)
        return cluster_points(raw_real, 1e-2) 
    end

    # Bertini_real Machinery for Bounding Boxes
    function find_boundary_points(f_expr, vars, x_range, y_range)
        x_var, y_var = vars
        boundary_pts = Vector{Vector{Float64}}()

        # Solve for x fixed at boundaries: find y intersections
        for val in x_range
            sys = System([subs(f_expr, x_var => val)], variables=[y_var])
            sols = real_solutions(solve(sys; show_progress=false))
            for s in sols
                if y_range[1] <= s[1] <= y_range[2]
                    push!(boundary_pts, [val, s[1]])
                end
            end
        end

        # Solve for y fixed at boundaries: find x intersections
        for val in y_range
            sys = System([subs(f_expr, y_var => val)], variables=[x_var])
            sols = real_solutions(solve(sys; show_progress=false))
            for s in sols
                if x_range[1] <= s[1] <= x_range[2]
                    push!(boundary_pts, [s[1], val])
                end
            end
        end
        return cluster_points(boundary_pts, 1e-3)
    end

    function cluster_points(points, tol)
        unique_pts = Vector{Vector{Float64}}()
        for pt in points
            if isempty(unique_pts) || all(u -> norm(pt - u) > tol, unique_pts)
                push!(unique_pts, pt)
            end
        end
        sort!(unique_pts, by = x -> (round(x[1], digits=6), round(x[2], digits=6)))
        return unique_pts
    end

    function check_singularity(f_expr, vars, point; tol=1e-3)
        # 1. Create a System object
        F = System([f_expr], variables=vars)
        
        # 2. Evaluate the Jacobian numerically at the point
        J_eval = jacobian(F, point)
        
        # 3. Now svdvals will work because J_eval contains numbers, not symbols
        s = svdvals(J_eval)
        
        # Numerical Rank Check: 
        # If all singular values are smaller than our tolerance, the rank is 0.
        return all(val < tol for val in s)
    end

    # --- 3D SURFACE MACHINERY ---

    function find_critical_z_slices(f_expr, vars)
        x, y, z = vars
        # The system for critical points in 3D
        # This finds points where the surface cannot be expressed as z = g(x,y)
        crit_sys = System([
            f_expr, 
            differentiate(f_expr, x), 
            differentiate(f_expr, y)
        ], variables=[x, y, z])
        
        result = solve(crit_sys; show_progress=false)
        sols = real_solutions(result)
        
        # Extract unique Z-coordinates
        z_values = [s[3] for s in sols]
        return sort(cluster_points_1d(z_values, 1e-5))
    end

    function cluster_points_1d(points, tol)
        if isempty(points) return Float64[] end
        sorted = sort(points)
        clusters = [sorted[1]]
        for i in 2:length(sorted)
            if sorted[i] - clusters[end] > tol
                push!(clusters, sorted[i])
            end
        end
        return clusters
    end

    # The bridge that passes a 2D slice to the 1D Engine
    function slice_and_extract_2d(f_expr, vars, z_val)
        x, y, z = vars
        
        # Substitute z with our slice height to get g(x,y)
        f_2d = subs(f_expr, z => z_val)
        
        return f_2d
    end
end