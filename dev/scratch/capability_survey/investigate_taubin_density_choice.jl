# Follow-up #4: test an intermediate density (edge=20, midslice=25) for
# the taubin_singular_structure.pdf regeneration decision, and get REAL
# (not estimated) PDF file sizes for intermediate + production density so
# the three-way comparison against the current edge=8 figure is rigorous.

using HomotopyContinuation
using HomotopyGetsReal
using GeometryBasics
using CairoMakie

const RENDER_DIR = "/Users/juancagc/HomotopyGetsReal/dev/scratch/capability_survey/renders"

naked_count(mesh) = length(HomotopyGetsReal._naked_mesh_edges(mesh))

@var x y z
Ftaubin = System([(x^2 + (1.2 * y)^2 + z^2 - 1)^3 - x^2 * z^3 - 0.1 * (1.2 * y)^2 * z^3], variables = [x, y, z])

# Same matched-zoom plotting cfg as the previous follow-up, for direct
# visual comparability against the two already-produced matched-zoom
# renders.
zoom_cfg = HomotopyConfig{Float64}(bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3))

function run_density(label, edge_density, midslice_density, save_pdf_stem)
    cfg = HomotopyConfig{Float64}(
        bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3),
        edge_sample_density = edge_density, midslice_sample_density = midslice_density,
    )
    println("\n=== $label (edge=$edge_density, midslice=$midslice_density) ===")
    t0 = time()
    v, e, f, m = decompose_3d_surface(Ftaubin, cfg; incidence = true)
    wall = time() - t0
    nk = naked_count(m)
    ntris = length(GeometryBasics.faces(m))
    nverts = length(GeometryBasics.coordinates(m))
    println("$label: wall=$(round(wall,digits=2))s  mesh_verts=$nverts  mesh_tris=$ntris  naked_edges=$nk")

    fig = plot_surface_decomposition(m; color_by = :z, cfg = zoom_cfg)
    png_path = joinpath(RENDER_DIR, "$(save_pdf_stem)_MATCHED_ZOOM.png")
    CairoMakie.save(png_path, fig)
    println("Saved PNG: $png_path")

    pdf_path = joinpath(RENDER_DIR, "$(save_pdf_stem)_MATCHED_ZOOM.pdf")
    CairoMakie.save(pdf_path, fig; pt_per_unit = 1)
    pdf_bytes = filesize(pdf_path)
    println("Saved PDF: $pdf_path  ($(round(pdf_bytes/1024/1024, digits=3)) MB, $pdf_bytes bytes)")

    return (wall = wall, naked_edges = nk, mesh_tris = ntris, mesh_verts = nverts, pdf_bytes = pdf_bytes, pdf_path = pdf_path, png_path = png_path)
end

r_intermediate = run_density("INTERMEDIATE", 20, 25, "taubin_heart_intermediate_density")
r_production = run_density("PRODUCTION (re-measured for a real PDF size)", 50, 100, "taubin_heart_production_density")

println("\n=== THREE-WAY SUMMARY (current edge=8 figure checked separately via `ls`) ===")
println("intermediate (edge=20,midslice=25): naked_edges=$(r_intermediate.naked_edges)  tris=$(r_intermediate.mesh_tris)  pdf=$(round(r_intermediate.pdf_bytes/1024/1024,digits=3))MB  wall=$(round(r_intermediate.wall,digits=1))s")
println("production   (edge=50,midslice=100): naked_edges=$(r_production.naked_edges)  tris=$(r_production.mesh_tris)  pdf=$(round(r_production.pdf_bytes/1024/1024,digits=3))MB  wall=$(round(r_production.wall,digits=1))s")
println("ALL DONE")
