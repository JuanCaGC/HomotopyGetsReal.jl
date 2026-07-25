# dev/scratch/scratch_paper_final_renders.jl
#
# Polished, publication-quality renders of Curves 1-3 for the paper
# (Part 1 of a two-part rendering request, 2026-07). Tight auto-fit axis
# limits (not the full configured bbox), DataAspect (already
# plot_curve_decomposition's default), edge_color_by comparison for
# Curves 1-2, an inset panel for Curve 2's separate oval component, and a
# two-view treatment of Curve 3's confirmed unbounded structure (zoomed
# finite structure + wide ray view). High-DPI PNG via GLMakie's
# px_per_unit (no CairoMakie dependency added just for vector export --
# not in this project's Project.toml, and pulling it in for a one-off
# figure-export task felt like the wrong tradeoff against a well-
# established high-DPI PNG fallback the user's own instructions allowed for).

using HomotopyContinuation
using LinearAlgebra
using GeometryBasics
using GLMakie
using HomotopyGetsReal

const OUTDIR = joinpath(@__DIR__, "renders_paper_final")
mkpath(OUTDIR)

const TITLESIZE = 18
const LABELSIZE = 14
const MARKERSIZE = 14
const LINEWIDTH = 2.5

@var x y
f1 = (x^2 + y^2 - 1)^3 + 27 * x^2 * y^2
f2 = (x^3 - x * y^2 + y + 1)^2 * (x^2 + y^2 - 1) + y^2 - 5
f3 = x^8 - 28 * x^6 * y^2 + 70 * x^4 * y^4 - 28 * x^2 * y^6 + y^8 + 15 * x^4 * y^2 - 15 * x^2 * y^4

cfg = HomotopyConfig{Float64}()
F1 = System([f1], variables = [x, y])
F2 = System([f2], variables = [x, y])
F3 = System([f3], variables = [x, y])

v1, e1 = decompose_1d_curve(F1, cfg)
v2, e2 = decompose_1d_curve(F2, cfg)
v3, e3 = decompose_1d_curve(F3, cfg)                      # default bbox (-4,4): finite/near-origin structure
cfg_wide = HomotopyConfig{Float64}(bbox_x = (-12.0, 12.0), bbox_y = (-12.0, 12.0))
v3w, e3w = decompose_1d_curve(F3, cfg_wide)               # wide bbox: shows the ray artifact

"""
    draw_curve!(ax, vertices, edges; kwargs...)

Mirrors plot_curve_decomposition's own internals (same edge-coloring
scheme, same HomotopyGetsReal._plot_vertices_by_type! vertex styling) but
draws into a caller-supplied Axis instead of building its own Figure --
needed for multi-panel figures (Curve 2's inset) that
plot_curve_decomposition itself doesn't support.
"""
function draw_curve!(
    ax, vertices::Vector{<:NativeVertex}, edges::Vector{<:Edge};
    edge_color_by::Symbol = :cell, edge_color = :steelblue, colormap = :viridis,
    show_vertices::Bool = true, show_legend::Bool = true, legend_position = :rb,
)
    n_edges = length(edges)
    edge_colors = edge_color_by == :cell ? Makie.cgrad(colormap, max(n_edges, 2); categorical = true) : nothing
    for (i, e) in enumerate(edges)
        isempty(e.sampled_points) && continue
        xs = [p[1] for p in e.sampled_points]
        ys = [p[2] for p in e.sampled_points]
        c = edge_color_by == :cell ? edge_colors[i] : edge_color
        lines!(ax, xs, ys; color = c, linewidth = LINEWIDTH)
    end
    if show_vertices
        handles = HomotopyGetsReal._plot_vertices_by_type!(ax, vertices)
        show_legend && !isempty(handles) && axislegend(ax; position = legend_position, labelsize = LABELSIZE)
    end
end

"""
    tight_lims(vertices, edges; margin_frac=0.12) -> (xlims, ylims)

Auto-fit axis limits from the ACTUAL plotted data (vertex + edge sample
coordinates), not the configured bbox -- so each curve fills its panel.
"""
function tight_lims(vertices, edges; margin_frac = 0.12)
    xs, ys = Float64[], Float64[]
    for v in vertices
        push!(xs, real(v.coordinates[1])); push!(ys, real(v.coordinates[2]))
    end
    for e in edges, p in e.sampled_points
        push!(xs, p[1]); push!(ys, p[2])
    end
    xlo, xhi = extrema(xs); ylo, yhi = extrema(ys)
    m = margin_frac * max(xhi - xlo, yhi - ylo)
    m = m == 0 ? 0.5 : m
    return (xlo - m, xhi + m), (ylo - m, yhi + m)
end

"""
    connected_components(vertices, edges) -> Dict{Int,Vector{Int}}

Union-find over the vertex/edge graph -- used to isolate Curve 2's
separate small oval component from its main loop for the inset panel.
"""
function connected_components(vertices, edges)
    parent = Dict(v.id => v.id for v in vertices)
    function find(a)
        while parent[a] != a
            parent[a] = parent[parent[a]]
            a = parent[a]
        end
        return a
    end
    for e in edges
        ra, rb = find(e.left_vertex_id), find(e.right_vertex_id)
        ra != rb && (parent[ra] = rb)
    end
    comps = Dict{Int,Vector{Int}}()
    for v in vertices
        push!(get!(comps, find(v.id), Int[]), v.id)
    end
    return comps
end

println("=" ^ 70)
println("Curve 1: edge_color_by comparison (:cell vs :mono)")
println("=" ^ 70)
for (label, ecb) in (("cell", :cell), ("mono", :mono))
    fig = Figure(size = (750, 750))
    ax = Axis(fig[1, 1]; aspect = DataAspect(), xlabel = "x", ylabel = "y",
        title = "Curve 1: (x²+y²-1)³ + 27x²y² = 0", titlesize = TITLESIZE,
        xlabelsize = LABELSIZE, ylabelsize = LABELSIZE)
    draw_curve!(ax, v1, e1; edge_color_by = ecb, legend_position = :rt)
    xl, yl = tight_lims(v1, e1)
    xlims!(ax, xl...); ylims!(ax, yl...)
    path = joinpath(OUTDIR, "curve1_compare_$(label).png")
    GLMakie.save(path, fig; px_per_unit = 4)
    println("  saved $path")
end

println()
println("=" ^ 70)
println("Curve 2: edge_color_by comparison + component census (for inset placement)")
println("=" ^ 70)
comps2 = connected_components(v2, e2)
comp_sizes = sort(collect(comps2); by = kv -> length(kv[2]))
println("  connected components: ", [(k, length(ids)) for (k, ids) in comp_sizes])
small_root = first(comp_sizes)[1]
small_ids = Set(comp_sizes[1][2])
v2_small = filter(v -> v.id in small_ids, v2)
e2_small = filter(e -> (e.left_vertex_id in small_ids) && (e.right_vertex_id in small_ids), e2)
println("  small (oval) component: $(length(v2_small)) vertices, $(length(e2_small)) edges")
for v in v2_small
    println("    id=$(v.id) coords=$(round.(real.(v.coordinates); digits=3)) type=$(v.v_type)")
end

for (label, ecb) in (("cell", :cell), ("mono", :mono))
    fig = Figure(size = (750, 750))
    ax = Axis(fig[1, 1]; aspect = DataAspect(), xlabel = "x", ylabel = "y",
        title = "Curve 2: (x³-xy²+y+1)²(x²+y²-1) + y² - 5 = 0", titlesize = TITLESIZE,
        xlabelsize = LABELSIZE, ylabelsize = LABELSIZE)
    draw_curve!(ax, v2, e2; edge_color_by = ecb, legend_position = :rt)
    xl, yl = tight_lims(v2, e2)
    xlims!(ax, xl...); ylims!(ax, yl...)
    path = joinpath(OUTDIR, "curve2_compare_$(label)_noinset.png")
    GLMakie.save(path, fig; px_per_unit = 4)
    println("  saved $path")
end

println()
println("=" ^ 70)
println("Curve 3: zoomed (finite structure) vs wide (ray artifact) views")
println("=" ^ 70)
for (label, ecb) in (("cell", :cell), ("mono", :mono))
    fig = Figure(size = (750, 750))
    ax = Axis(fig[1, 1]; aspect = DataAspect(), xlabel = "x", ylabel = "y",
        title = "Curve 3, zoomed near origin (finite critical/singular structure)", titlesize = TITLESIZE - 2,
        xlabelsize = LABELSIZE, ylabelsize = LABELSIZE)
    draw_curve!(ax, v3, e3; edge_color_by = ecb, legend_position = :rt)
    xlims!(ax, -2.2, 2.2); ylims!(ax, -2.2, 2.2)
    path = joinpath(OUTDIR, "curve3_zoomed_$(label).png")
    GLMakie.save(path, fig; px_per_unit = 4)
    println("  saved $path")
end

fig3w = Figure(size = (750, 750))
axw = Axis(fig3w[1, 1]; aspect = DataAspect(), xlabel = "x", ylabel = "y",
    title = "Curve 3, wide view (unbounded real structure: 16 asymptotic rays)", titlesize = TITLESIZE - 2,
    xlabelsize = LABELSIZE, ylabelsize = LABELSIZE)
draw_curve!(axw, v3w, e3w; edge_color_by = :mono, legend_position = :rt)
xl, yl = tight_lims(v3w, e3w; margin_frac = 0.05)
xlims!(axw, xl...); ylims!(axw, yl...)
path = joinpath(OUTDIR, "curve3_wide_rays.png")
GLMakie.save(path, fig3w; px_per_unit = 4)
println("  saved $path")

println()
println("All comparison renders saved to $OUTDIR")
println("Next: inspect renders, pick edge_color_by per curve, finalize Curve 2's inset placement.")
