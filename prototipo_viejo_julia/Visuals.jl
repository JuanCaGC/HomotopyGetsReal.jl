module Visuals

using GLMakie
using ..Types

export plot_topology, plot_surface

# Your existing 1D Plotting Function
function plot_topology(vertices::Vector{NativeVertex}, edges::Vector{Edge})
    fig = Figure(size = (800, 800))
    ax = Axis(fig[1, 1], aspect = DataAspect(), title="1D Algebraic Topology")
    
    # Plot edges
    for edge in edges
        x_vals = [pt[1] for pt in edge.path]
        y_vals = [pt[2] for pt in edge.path]
        lines!(ax, x_vals, y_vals, linewidth=3)
    end
    
    # Plot vertices
    if !isempty(vertices)
        vx = [real(v.coordinates[1]) for v in vertices]
        vy = [real(v.coordinates[2]) for v in vertices]
        v_colors = [v.v_type == Singular ? :red : (v.v_type == Boundary ? :gray : :black) for v in vertices]
        scatter!(ax, vx, vy, color=v_colors, markersize=12)
    end
    
    display(fig)
    return fig
end

# NEW: The 3D Surface Plotter
function plot_surface(faces::Vector{Face})
    println("Rendering Solid 3D Surface Decomposition...")
    
    fig = Figure(size = (1000, 800))
    ax = Axis3(fig[1, 1], 
               aspect = :data, 
               elevation = pi/8, 
               azimuth = pi/4,
               title = "Solid 3D Cellular Algebraic Decomposition")
    
    # We can add a deep crimson :red for the heart here!
    cell_colors = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :gold]
    
    for (i, face) in enumerate(faces)
        X_flat = face.mesh_vertices[:, 1]
        Y_flat = face.mesh_vertices[:, 2]
        Z_flat = face.mesh_vertices[:, 3]
        
        # THE DIMENSION FIX: Read the exact grid dimensions from the engine
        if size(face.mesh_faces) == (1, 2)
            steps = face.mesh_faces[1, 1]
            num_path_points = face.mesh_faces[1, 2]
        else
            # Fallback for old saved runs
            z_start = Z_flat[1]
            num_path_points = count(z -> isapprox(z, z_start, atol=1e-8), Z_flat)
            steps = div(length(X_flat), num_path_points)
        end
        
        # ... (keep the grid reshaping code above exactly the same) ...
        X_grid = reshape(X_flat, num_path_points, steps)
        Y_grid = reshape(Y_flat, num_path_points, steps)
        Z_grid = reshape(Z_flat, num_path_points, steps)
        
        c_idx = mod1(i, length(cell_colors))
        
        # THE FIX: Create a 2D matrix of solid colors matching the grid size
        color_matrix = fill(cell_colors[c_idx], size(Z_grid))
        
        surface!(ax, X_grid, Y_grid, Z_grid, color=color_matrix, transparency=false)
    end
    
    display(fig)
    return fig
end

end