# Targeted follow-up #2: does taubin_heart_incidence.png's visible
# top-lobe seam persist at production density / the paper's own actual
# config, or is it a coarse-density-only artifact?

using HomotopyContinuation
using HomotopyGetsReal
using GeometryBasics
using CairoMakie

const RENDER_DIR = "/Users/juancagc/HomotopyGetsReal/dev/scratch/capability_survey/renders"

naked_count(mesh) = length(HomotopyGetsReal._naked_mesh_edges(mesh))

function report(label, v, e, f, m, wall)
    nk = naked_count(m)
    println("$label: wall=$(round(wall,digits=2))s  vertices=$(length(v))  edges=$(length(e))  faces=$(length(f))  mesh_verts=$(length(GeometryBasics.coordinates(m)))  mesh_tris=$(length(GeometryBasics.faces(m)))  naked_edges=$nk")
    return nk
end

@var x y z
Ftaubin = System([(x^2 + (1.2 * y)^2 + z^2 - 1)^3 - x^2 * z^3 - 0.1 * (1.2 * y)^2 * z^3], variables = [x, y, z])

# --- Run A: literal HomotopyConfig{Float64}() full defaults (as literally
#     requested), incidence=true ---
cfg_defaults = HomotopyConfig{Float64}()
println("=== TAUBIN full-defaults config: edge_sample_density=$(cfg_defaults.edge_sample_density), midslice_sample_density=$(cfg_defaults.midslice_sample_density), bbox_x=$(cfg_defaults.bbox_x), bbox_y=$(cfg_defaults.bbox_y), bbox_z=$(cfg_defaults.bbox_z) ===")
t0 = time()
v, e, f, m = decompose_3d_surface(Ftaubin, cfg_defaults; incidence = true)
nk_defaults = report("taubin_full_defaults_incidence", v, e, f, m, time() - t0)
fig = plot_surface_decomposition(m; color_by = :z, cfg = cfg_defaults)
CairoMakie.save(joinpath(RENDER_DIR, "taubin_heart_incidence_production_density.png"), fig)

# --- Run B: the ACTUAL paper-figure config, from
#     paper_artifacts/taubin_singular_structure_example.jl (read directly,
#     not assumed): edge_sample_density=8, midslice_sample_density=8,
#     tightened bbox. This is the one that's actually comparable to the
#     real existing PDF. ---
cfg_paper = HomotopyConfig{Float64}(
    bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3),
    edge_sample_density = 8, midslice_sample_density = 8,
)
println("\n=== TAUBIN actual-paper-figure config: edge_sample_density=8, midslice_sample_density=8, bbox=(-1.5,1.5)/(-1.5,1.5)/(-1.3,1.3) ===")
t0 = time()
v, e, f, m = decompose_3d_surface(Ftaubin, cfg_paper; incidence = true)
nk_paper = report("taubin_paper_config_incidence", v, e, f, m, time() - t0)
fig = plot_surface_decomposition(m; color_by = :z, cfg = cfg_paper)
CairoMakie.save(joinpath(RENDER_DIR, "taubin_heart_incidence_paper_config.png"), fig)

println("\n=== SUMMARY ===")
println("coarse survey config (edge=6, midslice=8, default bbox), incidence=true: naked_edges=29 (measured in the prior follow-up)")
println("full HomotopyConfig{Float64}() defaults (edge=50, midslice=100, default bbox), incidence=true: naked_edges=$nk_defaults")
println("actual paper-figure config (edge=8, midslice=8, tightened bbox), incidence=true: naked_edges=$nk_paper")
println("ALL DONE")
