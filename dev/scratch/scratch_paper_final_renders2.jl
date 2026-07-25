# dev/scratch/scratch_paper_final_renders2.jl
#
# FINAL versions after comparing :cell vs :mono on the first pass
# (scratch_paper_final_renders.jl): Curve 1 -> :mono (reads as ONE closed
# curve, which is the point -- :cell fragments it into 12 arbitrary
# colors that fight the "closed astroid" gestalt). Curve 2 -> :cell (the
# curve is NOT a single smooth loop -- 2 disconnected components with
# sharp turns -- so per-edge color usefully marks where the topology
# actually changes, and doesn't fight any "one closed shape" reading
# since there isn't one). Curve 3 wide view: the plain auto-tight-bbox
# pass came out absurdly tall/narrow because a few :endpoint_fallback
# Artificial vertices (the fixed-x-axis-sweep artifact investigated
# separately) land at |y| ~70, far past where the 16-ray structure itself
# needs to be shown -- switched to a manually chosen square window that
# shows the genuine unbounded-ray structure cleanly without that
# artifact dominating the frame.

using HomotopyContinuation
using LinearAlgebra
using GeometryBasics
using GLMakie
using HomotopyGetsReal

const OUTDIR = joinpath(@__DIR__, "renders_paper_final")
mkpath(OUTDIR)

const TITLESIZE = 18
const LABELSIZE = 14
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
v3, e3 = decompose_1d_curve(F3, cfg)
cfg_wide = HomotopyConfig{Float64}(bbox_x = (-12.0, 12.0), bbox_y = (-12.0, 12.0))
v3w, e3w = decompose_1d_curve(F3, cfg_wide)

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
println("Curve 1 FINAL (:mono)")
println("=" ^ 70)
fig1 = Figure(size = (750, 750))
ax1 = Axis(fig1[1, 1]; aspect = DataAspect(), xlabel = "x", ylabel = "y",
    title = "Curve 1: (x²+y²-1)³ + 27x²y² = 0", titlesize = TITLESIZE,
    xlabelsize = LABELSIZE, ylabelsize = LABELSIZE)
draw_curve!(ax1, v1, e1; edge_color_by = :mono, legend_position = :rt)
xl, yl = tight_lims(v1, e1)
xlims!(ax1, xl...); ylims!(ax1, yl...)
GLMakie.save(joinpath(OUTDIR, "curve1_final.png"), fig1; px_per_unit = 4)
println("  saved curve1_final.png")

println()
println("=" ^ 70)
println("Curve 2 FINAL (:cell + inset on the small oval component)")
println("=" ^ 70)
comps2 = connected_components(v2, e2)
comp_sizes = sort(collect(comps2); by = kv -> length(kv[2]))
small_ids = Set(comp_sizes[1][2])
v2_small = filter(v -> v.id in small_ids, v2)
e2_small = filter(e -> (e.left_vertex_id in small_ids) && (e.right_vertex_id in small_ids), e2)
xl_s, yl_s = tight_lims(v2_small, e2_small; margin_frac = 0.3)
println("  oval component bbox for inset: x=$xl_s y=$yl_s")

fig2 = Figure(size = (800, 750))
ax2 = Axis(fig2[1, 1]; aspect = DataAspect(), xlabel = "x", ylabel = "y",
    title = "Curve 2: (x³-xy²+y+1)²(x²+y²-1) + y² - 5 = 0", titlesize = TITLESIZE,
    xlabelsize = LABELSIZE, ylabelsize = LABELSIZE)
draw_curve!(ax2, v2, e2; edge_color_by = :cell, legend_position = :lb)
xl, yl = tight_lims(v2, e2)
xlims!(ax2, xl...); ylims!(ax2, yl...)

# Inset axis in the same grid cell, positioned in empty space near plot center
# (found by inspecting the main render: the region around (0,0) has no curve).
ax2_inset = Axis(fig2[1, 1];
    width = Relative(0.36), height = Relative(0.36),
    halign = 0.40, valign = 0.42,
    aspect = DataAspect(), backgroundcolor = :white,
    xticklabelsize = 10, yticklabelsize = 10, title = "oval component (zoom)", titlesize = 12)
draw_curve!(ax2_inset, v2_small, e2_small; edge_color_by = :cell, show_legend = false)
xlims!(ax2_inset, xl_s...); ylims!(ax2_inset, yl_s...)
translate!(ax2_inset.scene, 0, 0, 100)
translate!(ax2_inset.blockscene, 0, 0, 101)
GLMakie.save(joinpath(OUTDIR, "curve2_final.png"), fig2; px_per_unit = 4)
println("  saved curve2_final.png")

println()
println("=" ^ 70)
println("Curve 3 FINAL: zoomed (Boundary type filtered out -- off-frame at this")
println("zoom, would otherwise show an empty, confusing legend entry) + wide")
println("(manually squared window, avoiding the fixed-axis-sweep artifact outliers)")
println("=" ^ 70)
v3_nozero = filter(v -> v.v_type != Boundary, v3)
fig3z = Figure(size = (750, 750))
ax3z = Axis(fig3z[1, 1]; aspect = DataAspect(), xlabel = "x", ylabel = "y",
    title = "Curve 3, zoomed near origin: finite critical/singular structure", titlesize = TITLESIZE - 2,
    xlabelsize = LABELSIZE, ylabelsize = LABELSIZE)
draw_curve!(ax3z, v3_nozero, e3; edge_color_by = :mono, legend_position = :rt)
xlims!(ax3z, -2.2, 2.2); ylims!(ax3z, -2.2, 2.2)
GLMakie.save(joinpath(OUTDIR, "curve3_zoomed_final.png"), fig3z; px_per_unit = 4)
println("  saved curve3_zoomed_final.png")

fig3w = Figure(size = (750, 750))
ax3w = Axis(fig3w[1, 1]; aspect = DataAspect(), xlabel = "x", ylabel = "y",
    title = "Curve 3, wide view: unbounded real structure (16 asymptotic rays)", titlesize = TITLESIZE - 2,
    xlabelsize = LABELSIZE, ylabelsize = LABELSIZE)
draw_curve!(ax3w, v3w, e3w; edge_color_by = :mono, legend_position = :rt)
xlims!(ax3w, -13.5, 13.5); ylims!(ax3w, -13.5, 13.5)
GLMakie.save(joinpath(OUTDIR, "curve3_wide_final.png"), fig3w; px_per_unit = 4)
println("  saved curve3_wide_final.png")

println()
println("All final renders in $OUTDIR")
