module HomotopyGetsReal

using HomotopyContinuation, LinearAlgebra

# Include the machinery modules we've refined
include("Types.jl")
include("Solver.jl")
include("Topology.jl")
include("SurfaceTopology.jl")
include("Visuals.jl")

using .Types, .Solver, .Topology, .SurfaceTopology, .Visuals

# NEW: Export the torus example
export run_astroid_example, run_unbounded_example, run_torus_example, run_sphere_example, run_heart_example, plot_surface

"""
    run_astroid_example()
Demonstrates the machinery on a bounded curve with singular cusps.
"""
function run_astroid_example()
    @var x y
    f = (x^2 + y^2 - 1)^3 + 27*x^2*y^2
    
    # 1. Define the Machinery's Bounding Box (The viewing window)
    x_range = (-1.5, 1.5)
    y_range = (-1.5, 1.5)

    println("Phase 1: Finding Critical Points...")
    sys = System([f, differentiate(f, y)], variables=[x, y])
    crit_pts = find_critical_points(sys)
    
    println("Phase 2: Finding Boundary Intersections...")
    bound_pts = find_boundary_points(f, [x, y], x_range, y_range)
    
    # 3. Classify and Create the Vertex Set
    vertices = NativeVertex[]
    id_counter = 1
    
    # Process Critical/Singular Points
    for p in crit_pts
        is_sing = check_singularity(f, [x, y], p)
        v_type = is_sing ? Singular : Critical
        push!(vertices, NativeVertex(id_counter, ComplexF64.(p), v_type))
        id_counter += 1
    end
    
    # Process Boundary Points
    for p in bound_pts
        push!(vertices, NativeVertex(id_counter, ComplexF64.(p), Boundary))
        id_counter += 1
    end

    println("Phase 3: Establishing Topology (Sweep-line)...")
    # n_samples=30 for extremely high-fidelity Bertini-style arcs
    edges = generate_1d_topology(vertices, f, x, y; n_samples=30)

    println("\nMachinery Summary:")
    println(" - Vertices: $(length(vertices)) ($(count(v->v.v_type==Singular, vertices)) Singular)")
    println(" - Edges:    $(length(edges))")

    # 4. Visual Analysis
    plot_topology(vertices, edges)
    
    return vertices, edges
end

"""
    run_unbounded_example()
Demonstrates how the Bounding Box 'clips' an infinite hyperbola.
"""
function run_unbounded_example()
    @var x y
    f = x^2 - y^2 - 1  # Hyperbola
    x_range = (-3.0, 3.0)
    y_range = (-3.0, 3.0)

    # Simplified pipeline for demonstration
    sys = System([f, differentiate(f, y)], variables=[x, y])
    crit_pts = find_critical_points(sys)
    bound_pts = find_boundary_points(f, [x, y], x_range, y_range)
    
    vertices = NativeVertex[]
    id_counter = 1
    for p in crit_pts
        push!(vertices, NativeVertex(id_counter, ComplexF64.(p), Critical))
        id_counter += 1
    end
    for p in bound_pts
        push!(vertices, NativeVertex(id_counter, ComplexF64.(p), Boundary))
        id_counter += 1
    end

    edges = generate_1d_topology(vertices, f, x, y; n_samples=20)
    plot_topology(vertices, edges)
end

"""
Demonstrates the 3D Cellular Algebraic Decomposition machinery on a vertical torus.
"""

function run_torus_example()
    @var x y z
    
    # The mathematical definition of the vertical torus
    f = (x^2 + y^2 + z^2 + 3)^2 - 16*(x^2 + z^2)

    println("Phase 1 & 2: Sweeping 1D Slices into 3D Faces...")
    
    # Call YOUR exact function signature!
    # We pass the expanded z_range as a keyword argument to capture the caps.
    faces = generate_3d_topology(f, [x, y, z]; z_range=(-3.5, 3.5))
    plot_surface(faces)
    return faces
end

"""
    run_sphere_example()
Demonstrates the 3D CAD machinery on a unit sphere to verify surface tracking.
"""
function run_sphere_example()
    @var x y z
    
    # 1. The mathematical definition of the unit sphere
    f = x^2 + y^2 + z^2 - 1
    
    println("--- Testing Pipeline with Unit Sphere ---")
    
    # 2. Call the 3D topology generator. 
    # Bounding box is set slightly wider than the radius of 1.
    faces = generate_3d_topology(f, [x, y, z]; z_range=(-1.5, 1.5))
    plot_surface(faces)
    return faces
end

"""
    run_heart_example()
Demonstrates the 3D CAD machinery on the Taubin Heart Surface.
Tests sharp cusps and deep topological cleavages.
"""
function run_heart_example()
    @var x y z
    
    # The mathematical definition of the Taubin Heart Surface
    f = (x^2 + (9/4)*y^2 + z^2 - 1)^3 - x^2 * z^3 - (9/80) * y^2 * z^3

    println("--- Testing Pipeline with the Taubin Heart ---")
    
    # The heart fits nicely within Z = [-1.5, 1.5]
    faces = generate_3d_topology(f, [x, y, z]; z_range=(-3, 3))
    plot_surface(faces)
    return faces
end

end # module