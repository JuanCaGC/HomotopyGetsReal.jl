# End-to-end smoke test: Oscar -> oscar_ideal_to_system -> real
# decompose_3d_surface call, residuals checked. Deferred from the
# previous round (blocked by display sleep); screen confirmed active now.

using Oscar, HomotopyContinuation, HomotopyGetsReal, OscarHomotopyContinuation
using GeometryBasics
using Statistics
using CairoMakie
using HomotopyContinuation: evaluate

include("/Users/juancagc/HomotopyGetsReal/examples/oscar_integration.jl")

println("="^70)
println("Building ellipsoid via Oscar, ring order [c_var, a_var, b_var]")
println("="^70)
R, (c_var, a_var, b_var) = Oscar.polynomial_ring(Oscar.QQ, ["c_var", "a_var", "b_var"])
f = a_var^2 + 4 * b_var^2 + 9 * c_var^2 - 1
I = Oscar.ideal(R, [f])
F = oscar_ideal_to_system(I)

println("F.variables (right before decompose call): ", F.variables)
println("gens(base_ring(I)):                          ", Oscar.gens(Oscar.base_ring(I)))
order_ok = string.(F.variables) == string.(Oscar.gens(Oscar.base_ring(I)))
println("Matches expected [c_var, a_var, b_var] order? ", order_ok)

println("\n" * "="^70)
println("Running decompose_3d_surface(F, cfg)")
println("="^70)
cfg = HomotopyConfig{Float64}()

ran_ok = true
local vertices, edges, faces, mesh
try
    global vertices, edges, faces, mesh = decompose_3d_surface(F, cfg)
    println("Ran without error.")
catch e
    global ran_ok = false
    println("FAILED:")
    showerror(stdout, e, catch_backtrace())
    println()
end

if ran_ok
    println("\n" * "="^70)
    println("Counts")
    println("="^70)
    println("vertices=$(length(vertices))  edges=$(length(edges))  faces=$(length(faces))")
    pts = GeometryBasics.coordinates(mesh)
    tris = GeometryBasics.faces(mesh)
    println("mesh_vertices=$(length(pts))  mesh_triangles=$(length(tris))")

    println("\n" * "="^70)
    println("Residuals: |F(p)| over every mesh point, F.variables order used directly")
    println("="^70)
    resids = Float64[abs(only(evaluate(F.expressions, F.variables => Float64.(p)))) for p in pts]
    println("n=$(length(resids))")
    println("mean=$(Statistics.mean(resids))")
    println("median=$(Statistics.median(resids))")
    println("p90=$(Statistics.quantile(resids, 0.9))")
    println("max=$(maximum(resids))")

    println("\n" * "="^70)
    println("Saving render")
    println("="^70)
    fig = plot_surface_decomposition(mesh; color_by = :z, cfg = cfg)
    render_path = "/Users/juancagc/HomotopyGetsReal/dev/scratch/oscar_investigation/end_to_end_ellipsoid.png"
    CairoMakie.save(render_path, fig)
    println("Saved: $render_path")

    println("\nRAN_OK=true  ORDER_OK=$order_ok  MEDIAN_RESIDUAL=$(Statistics.median(resids))")
else
    println("\nRAN_OK=false")
end

println("\nSMOKE TEST DONE")
