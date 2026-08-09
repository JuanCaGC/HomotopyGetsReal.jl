# Follow-up #3: the paper-config vs production-density comparison in
# investigate_taubin_seam.jl was confounded by different plot axis
# limits (plot_surface_decomposition sets xlims!/ylims!/zlims! from
# cfg.bbox_x/y/z). Re-render both meshes through an IDENTICAL, tight
# plotting cfg so zoom level cannot be doing the work of making one
# defect look bigger than the other.

using HomotopyContinuation
using HomotopyGetsReal
using GeometryBasics
using CairoMakie

const RENDER_DIR = "/Users/juancagc/HomotopyGetsReal/dev/scratch/capability_survey/renders"

naked_count(mesh) = length(HomotopyGetsReal._naked_mesh_edges(mesh))

@var x y z
Ftaubin = System([(x^2 + (1.2 * y)^2 + z^2 - 1)^3 - x^2 * z^3 - 0.1 * (1.2 * y)^2 * z^3], variables = [x, y, z])

# Identical, tight plotting cfg for BOTH renders -- controls axis limits
# only (plot_surface_decomposition uses cfg.bbox_x/y/z for xlims!/ylims!/
# zlims!; it does not re-derive anything about the mesh itself, which is
# already computed). Matches the actual paper figure's own bbox exactly.
zoom_cfg = HomotopyConfig{Float64}(bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3))

# --- Paper-config mesh (own decomposition cfg = its own tight bbox, same
#     as before) ---
cfg_paper = HomotopyConfig{Float64}(
    bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3),
    edge_sample_density = 8, midslice_sample_density = 8,
)
println("=== Recomputing paper-config mesh (edge=8, midslice=8, tight bbox) ===")
t0 = time()
v, e, f, m_paper = decompose_3d_surface(Ftaubin, cfg_paper; incidence = true)
t_paper = time() - t0
nk_paper = naked_count(m_paper)
println("paper_config: wall=$(round(t_paper,digits=2))s  mesh_verts=$(length(GeometryBasics.coordinates(m_paper)))  mesh_tris=$(length(GeometryBasics.faces(m_paper)))  naked_edges=$nk_paper")

fig_paper = plot_surface_decomposition(m_paper; color_by = :z, cfg = zoom_cfg)
CairoMakie.save(joinpath(RENDER_DIR, "taubin_heart_paper_config_MATCHED_ZOOM.png"), fig_paper)
println("Saved: taubin_heart_paper_config_MATCHED_ZOOM.png")

# --- Production-density mesh (own decomposition cfg = its own wide
#     default bbox, same as before -- decomposition itself is unaffected
#     by what we later plot with) ---
cfg_defaults = HomotopyConfig{Float64}()
println("\n=== Recomputing production-density mesh (edge=50, midslice=100, default wide bbox) ===")
t0 = time()
v, e, f, m_prod = decompose_3d_surface(Ftaubin, cfg_defaults; incidence = true)
t_prod = time() - t0
nk_prod = naked_count(m_prod)
println("production_density: wall=$(round(t_prod,digits=2))s  mesh_verts=$(length(GeometryBasics.coordinates(m_prod)))  mesh_tris=$(length(GeometryBasics.faces(m_prod)))  naked_edges=$nk_prod")

fig_prod = plot_surface_decomposition(m_prod; color_by = :z, cfg = zoom_cfg)
CairoMakie.save(joinpath(RENDER_DIR, "taubin_heart_production_density_MATCHED_ZOOM.png"), fig_prod)
println("Saved: taubin_heart_production_density_MATCHED_ZOOM.png")

println("\n=== SUMMARY (both rendered through the SAME zoom_cfg axis limits: bbox_x=bbox_y=(-1.5,1.5), bbox_z=(-1.3,1.3)) ===")
println("paper_config (edge=8,midslice=8):        naked_edges=$nk_paper  wall=$(round(t_paper,digits=2))s")
println("production_density (edge=50,midslice=100): naked_edges=$nk_prod  wall=$(round(t_prod,digits=2))s")
println("ALL DONE")
