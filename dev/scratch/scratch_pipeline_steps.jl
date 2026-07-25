# dev/scratch/scratch_pipeline_steps.jl
#
# Pedagogical step-by-step illustration of decompose_1d_curve's pipeline
# (Part 2 of the paper-figures request, 2026-07), using Curve 1 (the
# astroid) as the running example -- cleanest topology, and already the
# subject of Part 1's final render for direct before/after comparison.
#
# Proposed 6-step sequence (as suggested by the user, adopted as-is: it
# maps cleanly onto compute_critical_points / intersect_bounding_object /
# compute_midslice / connect_the_dots! / decompose_1d_curve, each already
# a real, independently-callable pipeline function -- no awkward
# reach-into-internals needed for any step):
#   1. Raw curve: dense implicit contour of f=0, nothing computed.
#   2. + critical points (compute_critical_points on the y-derivative-
#      augmented system, exactly what decompose_1d_curve itself builds).
#   3. + boundary points (intersect_bounding_object).
#   4. One interval's midslice witness point(s) (compute_midslice),
#      isolated to a single representative [x_left, x_right] interval.
#   5. That same interval's tracked edge(s) (connect_the_dots!), shown
#      against the still-visible raw curve so the reader can see the
#      tracked polyline actually follows it.
#   6. The final assembled decomposition (decompose_1d_curve, full run) --
#      directly comparable to step 1's raw curve and to Part 1's polished
#      curve1_final.png.

using HomotopyContinuation
using LinearAlgebra
using GeometryBasics
using GLMakie
using HomotopyGetsReal

const OUTDIR = joinpath(@__DIR__, "renders_pipeline_steps")
mkpath(OUTDIR)

const TITLESIZE = 18
const LABELSIZE = 14
const LINEWIDTH = 2.5
const XLIM = (-2.9, 2.9)
const YLIM = (-2.9, 2.9)

@var xv yv
f1 = (xv^2 + yv^2 - 1)^3 + 27 * xv^2 * yv^2
cfg = HomotopyConfig{Float64}()
F1 = System([f1], variables = [xv, yv])

function new_axis(title)
    fig = Figure(size = (750, 750))
    ax = Axis(fig[1, 1]; aspect = DataAspect(), xlabel = "x", ylabel = "y",
        title = title, titlesize = TITLESIZE, xlabelsize = LABELSIZE, ylabelsize = LABELSIZE)
    xlims!(ax, XLIM...); ylims!(ax, YLIM...)
    return fig, ax
end

function draw_implicit_curve!(ax, fexpr, xvar, yvar; n = 400, color = :gray40, linewidth = LINEWIDTH)
    xs = range(XLIM...; length = n)
    ys = range(YLIM...; length = n)
    Z = [real(HomotopyContinuation.evaluate(fexpr, [xvar, yvar] => [xg, yg])) for xg in xs, yg in ys]
    contour!(ax, xs, ys, Z; levels = [0.0], color = color, linewidth = linewidth)
end

function draw_curve!(ax, vertices, edges; edge_color = :steelblue, show_legend = true, legend_position = :rt)
    for e in edges
        isempty(e.sampled_points) && continue
        xs = [p[1] for p in e.sampled_points]
        ys = [p[2] for p in e.sampled_points]
        lines!(ax, xs, ys; color = edge_color, linewidth = LINEWIDTH)
    end
    handles = HomotopyGetsReal._plot_vertices_by_type!(ax, vertices)
    show_legend && !isempty(handles) && axislegend(ax; position = legend_position, labelsize = LABELSIZE)
end

println("=" ^ 70)
println("Step 1: raw curve, nothing computed")
println("=" ^ 70)
fig1, ax1 = new_axis("Step 1: the raw curve f(x,y) = 0\n(x²+y²-1)³ + 27x²y² = 0, nothing computed yet")
draw_implicit_curve!(ax1, f1, xv, yv)
GLMakie.save(joinpath(OUTDIR, "01_raw_curve.png"), fig1; px_per_unit = 4)
println("  saved 01_raw_curve.png")

println()
println("=" ^ 70)
println("Step 2: critical points (compute_critical_points on the y-derivative-augmented system)")
println("=" ^ 70)
F_aug = System([f1, differentiate(f1, yv)], variables = [xv, yv])
crit_vertices = HomotopyGetsReal.compute_critical_points(F_aug, cfg)
println("  found $(length(crit_vertices)) critical vertices")
for v in crit_vertices
    println("    id=$(v.id) coords=$(round.(real.(v.coordinates);digits=4)) type=$(v.v_type)")
end
fig2, ax2 = new_axis("Step 2: critical points\ncompute_critical_points(F_aug, cfg) -- vertical-tangent + node points")
draw_implicit_curve!(ax2, f1, xv, yv)
draw_curve!(ax2, crit_vertices, HomotopyGetsReal.Edge{Float64}[]; legend_position = :rt)
GLMakie.save(joinpath(OUTDIR, "02_critical_points.png"), fig2; px_per_unit = 4)
println("  saved 02_critical_points.png")

println()
println("=" ^ 70)
println("Step 3: boundary points (intersect_bounding_object)")
println("=" ^ 70)
bnd_vertices = HomotopyGetsReal.intersect_bounding_object(F1, cfg)
println("  found $(length(bnd_vertices)) boundary vertices (bbox = $(cfg.bbox_x) x $(cfg.bbox_y))")
fig3, ax3 = new_axis("Step 3: boundary points\nintersect_bounding_object(F, cfg) -- curve ∩ bbox edges")
draw_implicit_curve!(ax3, f1, xv, yv)
draw_curve!(ax3, bnd_vertices, HomotopyGetsReal.Edge{Float64}[]; legend_position = :rt)
println("  (Curve 1 is fully contained in the default bbox (-4,4) -- expect 0 boundary points;")
println("   this step is included for pipeline completeness even though it's empty here.)")
GLMakie.save(joinpath(OUTDIR, "03_boundary_points.png"), fig3; px_per_unit = 4)
println("  saved 03_boundary_points.png")

println()
println("=" ^ 70)
println("Step 4: one interval's midslice witness point(s)")
println("=" ^ 70)
all_vertices_raw = vcat(crit_vertices, bnd_vertices)
all_vertices = HomotopyGetsReal.cluster_vertices(all_vertices_raw, cfg.vertex_match_tol)
xs_all = Float64[real(v.coordinates[1]) for v in all_vertices]
distinct_xs = HomotopyGetsReal.cluster_scalars(xs_all, cfg.vertex_match_tol)
sort!(distinct_xs)
println("  distinct critical/boundary x-values: ", round.(distinct_xs; digits = 4))
# Pick the interval strictly between x=0 and the next positive critical x
# (the upper-right petal's inner arc) -- an interior interval, not touching
# either end of the curve's x-range, for a representative illustration.
mid_interval_idx = findfirst(i -> distinct_xs[i] >= 0 && distinct_xs[i+1] > distinct_xs[i], 1:(length(distinct_xs)-1))
x_left, x_right = distinct_xs[mid_interval_idx], distinct_xs[mid_interval_idx+1]
x_mid = (x_left + x_right) / 2
println("  chosen illustration interval: x_left=$x_left, x_right=$x_right, x_mid=$x_mid")
y_mids = HomotopyGetsReal.compute_midslice(F1, x_left, x_right, cfg)
println("  compute_midslice found $(length(y_mids)) witness y-value(s) at x=$x_mid: ", y_mids)

fig4, ax4 = new_axis("Step 4: one interval's midslice witness\ncompute_midslice(F, x_left, x_right, cfg), x ∈ [$(round(x_left;digits=2)), $(round(x_right;digits=2))]")
draw_implicit_curve!(ax4, f1, xv, yv)
vlines!(ax4, [x_left, x_right]; color = :gray70, linestyle = :dash, linewidth = 1.5)
vlines!(ax4, [x_mid]; color = :black, linestyle = :dot, linewidth = 1.5)
scatter!(ax4, fill(x_mid, length(y_mids)), real.(y_mids); color = :purple, marker = :star5, markersize = 22, label = "midslice witness")
axislegend(ax4; position = :rt, labelsize = LABELSIZE)
GLMakie.save(joinpath(OUTDIR, "04_midslice_witness.png"), fig4; px_per_unit = 4)
println("  saved 04_midslice_witness.png")

println()
println("=" ^ 70)
println("Step 5: connect-the-dots for that interval (tracked edge reaching outward)")
println("=" ^ 70)
tracked_edges = HomotopyGetsReal.Edge{Float64}[]
for (k, y_mid) in enumerate(y_mids)
    e = HomotopyGetsReal.connect_the_dots!(F1, x_left, x_mid, x_right, y_mid, k, all_vertices, cfg)
    push!(tracked_edges, e)
    println("  tracked edge $k: left_vertex_id=$(e.left_vertex_id) right_vertex_id=$(e.right_vertex_id) $(length(e.sampled_points)) raw points")
end
fig5, ax5 = new_axis("Step 5: connect-the-dots for that interval\nconnect_the_dots! tracks bidirectionally from the midslice witness")
draw_implicit_curve!(ax5, f1, xv, yv; color = (:gray70, 0.6))
for e in tracked_edges
    xs = [p[1] for p in e.sampled_points]
    ys = [p[2] for p in e.sampled_points]
    lines!(ax5, xs, ys; color = :crimson, linewidth = LINEWIDTH + 1.5)
end
scatter!(ax5, fill(x_mid, length(y_mids)), real.(y_mids); color = :purple, marker = :star5, markersize = 22, label = "midslice witness")
relevant_vertex_ids = Set(vcat([e.left_vertex_id for e in tracked_edges], [e.right_vertex_id for e in tracked_edges]))
draw_curve!(ax5, filter(v -> v.id in relevant_vertex_ids, all_vertices), HomotopyGetsReal.Edge{Float64}[]; show_legend = false)
axislegend(ax5; position = :rt, labelsize = LABELSIZE)
GLMakie.save(joinpath(OUTDIR, "05_connect_the_dots.png"), fig5; px_per_unit = 4)
println("  saved 05_connect_the_dots.png")

println()
println("=" ^ 70)
println("Step 6: final assembled decomposition (decompose_1d_curve, full run)")
println("=" ^ 70)
v_final, e_final = decompose_1d_curve(F1, cfg)
println("  $(length(v_final)) vertices, $(length(e_final)) edges")
fig6, ax6 = new_axis("Step 6: final assembled decomposition\ndecompose_1d_curve(F, cfg) -- compare against Step 1's raw curve")
draw_curve!(ax6, v_final, e_final; legend_position = :rt)
GLMakie.save(joinpath(OUTDIR, "06_final_decomposition.png"), fig6; px_per_unit = 4)
println("  saved 06_final_decomposition.png")

println()
println("All pipeline-step renders saved to $OUTDIR")
