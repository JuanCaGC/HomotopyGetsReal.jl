# Candidate replacement for paper_artifacts/figures/paper/surfaces/taubin_singular_structure.pdf
# at the recommended intermediate density (edge_sample_density=20,
# midslice_sample_density=25), mirroring
# paper_artifacts/scripts/taubin_singular_structure_example.jl's own Part 1 logic
# (same camera pin, same vertex-type overlay, same composition) EXACTLY,
# except: (1) the new density, (2) output path is a new CANDIDATE file,
# never overwriting the original, (3) does NOT touch paper_artifacts/data/results.json
# -- this is a candidate for review, not a finalized regeneration.

using HomotopyContinuation
using HomotopyGetsReal
using GeometryBasics
using CairoMakie

const FIGURES_DIR = "/Users/juancagc/HomotopyGetsReal/paper_artifacts/figures/archive"

function vertex_type_counts(vertices)
    counts = Dict{String,Int}("Critical" => 0, "Boundary" => 0, "Singular" => 0, "Artificial" => 0)
    for v in vertices
        counts[string(v.v_type)] += 1
    end
    return counts
end

@var x y z
f = (x^2 + (1.2 * y)^2 + z^2 - 1)^3 - x^2 * z^3 - 0.1 * (1.2 * y)^2 * z^3
F = System([f], variables = [x, y, z])

# Recommended intermediate density; same tight bbox as the original figure.
cfg = HomotopyConfig{Float64}(
    bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3),
    edge_sample_density = 20, midslice_sample_density = 25,
)

println("="^70)
println("decompose_3d_surface(F, cfg; incidence=true) -- candidate at edge=20, midslice=25")
println("="^70)
t1 = @elapsed (vertices, edges, faces, mesh, inc) = decompose_3d_surface(F, cfg; incidence = true)
counts = vertex_type_counts(vertices)
counts_anchor = vertex_type_counts(inc.critical_vertices)
overlay_vertices = vcat(vertices, inc.critical_vertices)
counts_overlay = vertex_type_counts(overlay_vertices)
println("  decompose_3d_surface(incidence=true): $(round(t1; digits = 2))s")
println("  main vertices:       $(length(vertices)) total  $(counts)")
println("  fold-anchor vertices: $(length(inc.critical_vertices)) total  $(counts_anchor)")
println("  combined overlay:    $(length(overlay_vertices)) total  $(counts_overlay)")
println("  naked_edges: $(length(HomotopyGetsReal._naked_mesh_edges(mesh)))")
println("  mesh_verts=$(length(GeometryBasics.coordinates(mesh)))  mesh_tris=$(length(GeometryBasics.faces(mesh)))")
flush(stdout)

CairoMakie.activate!()

println()
println("="^70)
println("Candidate figure: mesh + combined VertexType overlay, same composition as taubin_singular_structure.pdf")
println("="^70)
fig1 = plot_surface_decomposition(mesh; color_by = :z, cfg = cfg, vertices = overlay_vertices)
ax1 = fig1.content[1]
ax1.title = "Taubin heart: welded mesh + VertexType overlay"
# Identical camera pin to the original script (kept as a literal there too).
ax1.azimuth[] = 1.275 * pi
ax1.elevation[] = pi / 8
ax1.viewmode[] = :fit

path1 = joinpath(FIGURES_DIR, "taubin_singular_structure_CANDIDATE.pdf")
CairoMakie.save(path1, fig1; pt_per_unit = 1)
png1 = replace(path1, ".pdf" => "_preview.png")
CairoMakie.save(png1, fig1; px_per_unit = 2)
pdf_bytes = filesize(path1)
println("  wrote $path1 ($(round(pdf_bytes/1024/1024, digits=3)) MB, $pdf_bytes bytes)")
println("  wrote $png1")
println("NOTE: paper_artifacts/data/results.json was NOT touched. The original")
println("taubin_singular_structure.pdf was NOT overwritten.")
println("ALL DONE")
