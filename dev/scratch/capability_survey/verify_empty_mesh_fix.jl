# Verify the plot_surface_decomposition empty-mesh fix against the three
# fixtures that actually crashed this way in the capability survey.

using HomotopyContinuation
using HomotopyGetsReal
using GeometryBasics
using CairoMakie

const RENDER_DIR = "/Users/juancagc/HomotopyGetsReal/dev/scratch/capability_survey/renders"

@var x y z

fixtures = [
    ("cone", System([x^2 + y^2 - z^2], variables = [x, y, z])),
    ("horn_torus", System([(x^2 + y^2 + z^2)^2 - 4 * (x^2 + y^2)], variables = [x, y, z])),
    ("empty_surface", System([x^2 + y^2 + z^2 + 1], variables = [x, y, z])),
]

cfg = HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8)

for (name, F) in fixtures
    println("\n=== $name ===")
    v, e, f, m = decompose_3d_surface(F, cfg)
    println("  vertices=$(length(v))  edges=$(length(e))  faces=$(length(f))  mesh_pts=$(length(GeometryBasics.coordinates(m)))  mesh_tris=$(length(GeometryBasics.faces(m)))")
    try
        fig = plot_surface_decomposition(m; color_by = :z, cfg = cfg)
        path = joinpath(RENDER_DIR, "$(name)_POST_FIX.png")
        CairoMakie.save(path, fig)
        println("  PASS: rendered without error -> $path")
    catch e
        println("  FAIL: still throws:")
        showerror(stdout, e, catch_backtrace())
        println()
    end
end

println("\nALL DONE")
