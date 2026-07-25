# dev/scratch/scratch_pipeline_steps_verify.jl
#
# Verification re-run of scratch_pipeline_steps.jl (2026-07): prints every
# piece of raw data requested for each of the 6 pipeline-step figures
# (not just the images), plus a density-convergence comparison for Part 2
# of the follow-up request. Reuses the exact same code path as the
# original script -- not a re-description, a live re-run.

using HomotopyContinuation
using LinearAlgebra
using GeometryBasics
using GLMakie
using HomotopyGetsReal

const OUTDIR = joinpath(@__DIR__, "renders_pipeline_steps")
mkpath(OUTDIR)
const XLIM = (-2.9, 2.9)
const YLIM = (-2.9, 2.9)

@var xv yv
f1 = (xv^2 + yv^2 - 1)^3 + 27 * xv^2 * yv^2
cfg = HomotopyConfig{Float64}()
F1 = System([f1], variables = [xv, yv])

println("=" ^ 70)
println("cfg fields actually in force for this whole script (HomotopyConfig{Float64}(), no overrides)")
println("=" ^ 70)
println("  edge_sample_density = ", cfg.edge_sample_density)
println("  midslice_sample_density = ", cfg.midslice_sample_density)
println("  bbox_x = ", cfg.bbox_x, "   bbox_y = ", cfg.bbox_y)
println("  vertex_match_tol = ", cfg.vertex_match_tol)

println()
println("=" ^ 70)
println("STEP 1: raw curve rendering -- exact method, no pipeline function")
println("=" ^ 70)
n_grid = 400
xs_grid = range(XLIM...; length = n_grid)
ys_grid = range(YLIM...; length = n_grid)
println("  grid: $(n_grid) x $(n_grid) = $(n_grid^2) points over x in $(XLIM), y in $(YLIM)")
println("  each grid point: Z[i,j] = real(HomotopyContinuation.evaluate(f1, [xv,yv] => [xs[i],ys[j]]))")
println("  contour method: Makie's contour!(ax, xs, ys, Z; levels=[0.0]) -- standard marching-squares")
println("  contouring on the SIGN of f1 alone. No NativeVertex/Edge/decompose_1d_curve/compute_* call")
println("  anywhere in this step -- confirmed by inspection: the only HC call is `evaluate`.")
Z_sample = [real(HomotopyContinuation.evaluate(f1, [xv, yv] => [xs_grid[1], ys_grid[1]]))]
println("  spot check, f1 at grid corner (x=$(xs_grid[1]), y=$(ys_grid[1])) = $(Z_sample[1])")

println()
println("=" ^ 70)
println("STEP 2: compute_critical_points(F_aug, cfg) -- full raw output")
println("=" ^ 70)
F_aug = System([f1, differentiate(f1, yv)], variables = [xv, yv])
crit_vertices = compute_critical_points(F_aug, cfg)
println("  F_aug = ", F_aug)
println("  returned $(length(crit_vertices)) NativeVertex objects:")
for v in sort(crit_vertices; by = x -> x.id)
    println("    id=$(v.id)  coordinates=$(v.coordinates)  v_type=$(v.v_type)")
end

println()
println("=" ^ 70)
println("STEP 3: intersect_bounding_object(F1, cfg) -- full raw output")
println("=" ^ 70)
bnd_vertices = intersect_bounding_object(F1, cfg)
println("  returned $(length(bnd_vertices)) NativeVertex objects (bbox_x=$(cfg.bbox_x), bbox_y=$(cfg.bbox_y))")
for v in bnd_vertices
    println("    id=$(v.id)  coordinates=$(v.coordinates)  v_type=$(v.v_type)")
end
println("  isempty(bnd_vertices) = ", isempty(bnd_vertices), "  (programmatic confirmation, not assumed)")

println()
println("=" ^ 70)
println("STEP 4: compute_midslice -- interval selection rule + full raw output")
println("=" ^ 70)
all_vertices_raw = vcat(crit_vertices, bnd_vertices)
all_vertices = cluster_vertices(all_vertices_raw, cfg.vertex_match_tol)
xs_all = Float64[real(v.coordinates[1]) for v in all_vertices]
distinct_xs = sort(cluster_scalars(xs_all, cfg.vertex_match_tol))
println("  distinct critical/boundary x-values (cluster_scalars output): ", round.(distinct_xs; digits = 6))
println("  selection rule (literal code): findfirst(i -> distinct_xs[i] >= 0 && distinct_xs[i+1] > distinct_xs[i], 1:(length(distinct_xs)-1))")
println("  i.e. \"the first interval whose left endpoint is >= 0\" -- a fixed, mechanical, curve-independent")
println("  rule, evaluated against whatever distinct_xs THIS curve happens to produce. Not hand-picked")
println("  per-curve to look good; for this specific curve's x-values it lands on [0,1].")
mid_interval_idx = findfirst(i -> distinct_xs[i] >= 0 && distinct_xs[i+1] > distinct_xs[i], 1:(length(distinct_xs)-1))
x_left, x_right = distinct_xs[mid_interval_idx], distinct_xs[mid_interval_idx+1]
x_mid = (x_left + x_right) / 2
println("  => x_left=$x_left, x_right=$x_right, x_mid=$x_mid")
y_mids = compute_midslice(F1, x_left, x_right, cfg)
println("  compute_midslice(F1, $x_left, $x_right, cfg) returned $(length(y_mids)) witness value(s):")
for y in y_mids
    println("    $y   (Float64 real part: $(real(y)))")
end

println()
println("=" ^ 70)
println("STEP 5: connect_the_dots! -- full raw Edge output for this interval")
println("=" ^ 70)
tracked_edges = Edge{Float64}[]
for (k, y_mid) in enumerate(y_mids)
    e = connect_the_dots!(F1, x_left, x_mid, x_right, y_mid, k, all_vertices, cfg)
    push!(tracked_edges, e)
end
for e in tracked_edges
    println("  Edge: id=$(e.id)  left_vertex_id=$(e.left_vertex_id)  right_vertex_id=$(e.right_vertex_id)  is_singular=$(e.is_singular)")
    println("    sampled_points: $(length(e.sampled_points)) points")
    println("    first point: $(e.sampled_points[1])")
    println("    last  point: $(e.sampled_points[end])")
    lv = only(filter(v -> v.id == e.left_vertex_id, all_vertices))
    rv = only(filter(v -> v.id == e.right_vertex_id, all_vertices))
    println("    left_vertex ($(e.left_vertex_id)) coords = $(round.(real.(lv.coordinates);digits=4)), type=$(lv.v_type)")
    println("    right_vertex ($(e.right_vertex_id)) coords = $(round.(real.(rv.coordinates);digits=4)), type=$(rv.v_type)")
end

println()
println("=" ^ 70)
println("STEP 6: decompose_1d_curve(F1, cfg) -- final summary + per-edge sample counts")
println("=" ^ 70)
v_final, e_final = decompose_1d_curve(F1, cfg)
n_crit = count(v -> v.v_type == Critical, v_final)
n_sing = count(v -> v.v_type == Singular, v_final)
n_bnd = count(v -> v.v_type == Boundary, v_final)
n_art = count(v -> v.v_type == Artificial, v_final)
println("  vertices: $(length(v_final)) total  (Critical=$n_crit, Boundary=$n_bnd, Singular=$n_sing, Artificial=$n_art)")
println("  edges: $(length(e_final))")
println("  critical points (Critical+Singular, this session's earlier convention) = $(n_crit + n_sing)")
matches = (n_crit + n_sing == 7) && (length(v_final) == 12) && (length(e_final) == 12)
println("  matches this session's earlier Curve 1 measurement (7 crit, 12 vertices, 12 edges)? ", matches)
println("  per-edge sampled_points counts (should each equal cfg.edge_sample_density=$(cfg.edge_sample_density)):")
for e in e_final
    flag = length(e.sampled_points) == cfg.edge_sample_density ? "" : "  <-- DOES NOT MATCH cfg.edge_sample_density!"
    println("    edge id=$(e.id): $(length(e.sampled_points)) points$flag")
end

println()
println("=" ^ 70)
println("PART 2 (follow-up): honest density comparison -- coarse(8) vs actual")
println("default(50, == Step 6 exactly) vs raw ground truth")
println("=" ^ 70)

const TITLESIZE = 16
const LABELSIZE = 12
const LINEWIDTH = 2.2

function draw_curve_only!(ax, vertices, edges; edge_color = :steelblue)
    for e in edges
        isempty(e.sampled_points) && continue
        xs = [p[1] for p in e.sampled_points]
        ys = [p[2] for p in e.sampled_points]
        lines!(ax, xs, ys; color = edge_color, linewidth = LINEWIDTH)
    end
end

cfg8 = HomotopyConfig{Float64}(edge_sample_density = 8)
v8, e8 = decompose_1d_curve(F1, cfg8)
println("  edge_sample_density=8:  $(length(v8)) vertices, $(length(e8)) edges, ",
        "per-edge samples: ", [length(e.sampled_points) for e in e8])

# density=50 case: reuse v_final/e_final already computed above (literally Step 6's own data)
println("  edge_sample_density=50: $(length(v_final)) vertices, $(length(e_final)) edges (this IS Step 6)")

fig = Figure(size = (2250, 750))
ax_a = Axis(fig[1, 1]; aspect = DataAspect(), title = "edge_sample_density = 8\n(genuinely coarse, for contrast)",
    titlesize = TITLESIZE, xlabelsize = LABELSIZE, ylabelsize = LABELSIZE, xlabel = "x", ylabel = "y")
draw_curve_only!(ax_a, v8, e8)
xlims!(ax_a, XLIM...); ylims!(ax_a, YLIM...)

ax_b = Axis(fig[1, 2]; aspect = DataAspect(), title = "edge_sample_density = 50 (actual default)\n= Step 6, exactly, unchanged",
    titlesize = TITLESIZE, xlabelsize = LABELSIZE, ylabelsize = LABELSIZE, xlabel = "x", ylabel = "y")
draw_curve_only!(ax_b, v_final, e_final)
xlims!(ax_b, XLIM...); ylims!(ax_b, YLIM...)

ax_c = Axis(fig[1, 3]; aspect = DataAspect(), title = "raw implicit curve\n(Step 1, ground truth)",
    titlesize = TITLESIZE, xlabelsize = LABELSIZE, ylabelsize = LABELSIZE, xlabel = "x", ylabel = "y")
xs_grid2 = range(XLIM...; length = n_grid)
ys_grid2 = range(YLIM...; length = n_grid)
Z2 = [real(HomotopyContinuation.evaluate(f1, [xv, yv] => [xg, yg])) for xg in xs_grid2, yg in ys_grid2]
contour!(ax_c, xs_grid2, ys_grid2, Z2; levels = [0.0], color = :gray30, linewidth = LINEWIDTH)
xlims!(ax_c, XLIM...); ylims!(ax_c, YLIM...)

path = joinpath(OUTDIR, "07_density_comparison.png")
GLMakie.save(path, fig; px_per_unit = 4)
println("  saved $path")

println()
println("=" ^ 70)
println("PART 2b: does a genuinely higher density close the gap to the raw curve?")
println("=" ^ 70)
cfg200 = HomotopyConfig{Float64}(edge_sample_density = 200)
v200, e200 = decompose_1d_curve(F1, cfg200)
println("  edge_sample_density=200: $(length(v200)) vertices, $(length(e200)) edges, ",
        "per-edge samples: ", [length(e.sampled_points) for e in e200])

fig2 = Figure(size = (1500, 750))
ax_d = Axis(fig2[1, 1]; aspect = DataAspect(), title = "edge_sample_density = 50 (actual default) = Step 6",
    titlesize = TITLESIZE, xlabelsize = LABELSIZE, ylabelsize = LABELSIZE, xlabel = "x", ylabel = "y")
draw_curve_only!(ax_d, v_final, e_final)
xlims!(ax_d, XLIM...); ylims!(ax_d, YLIM...)

ax_e = Axis(fig2[1, 2]; aspect = DataAspect(), title = "edge_sample_density = 200",
    titlesize = TITLESIZE, xlabelsize = LABELSIZE, ylabelsize = LABELSIZE, xlabel = "x", ylabel = "y")
draw_curve_only!(ax_e, v200, e200)
xlims!(ax_e, XLIM...); ylims!(ax_e, YLIM...)

path2 = joinpath(OUTDIR, "08_density_50_vs_200.png")
GLMakie.save(path2, fig2; px_per_unit = 4)
println("  saved $path2")

println()
println("=" ^ 70)
println("PART 2c: are the visible 'kinks' actually AT vertices (tangent mismatch")
println("between independently-tracked edges), not a within-edge sampling artifact?")
println("=" ^ 70)
println("  v200 vertex list:")
for v in sort(v200; by=x->x.id)
    println("    id=$(v.id) coords=$(round.(real.(v.coordinates);digits=4)) type=$(v.v_type)")
end
# Zoom on the artificial vertex near (1, 2.28) at high density to see if the
# kink survives right up close.
fig3 = Figure(size = (900, 900))
axz = Axis(fig3[1,1]; aspect=DataAspect(), title="edge_sample_density=200, zoomed on an Artificial-vertex junction",
    xlabel="x", ylabel="y")
draw_curve_only!(axz, v200, e200)
HomotopyGetsReal._plot_vertices_by_type!(axz, v200)
axislegend(axz; position=:rb)
xlims!(axz, 0.7, 1.3); ylims!(axz, 2.0, 2.5)
path3 = joinpath(OUTDIR, "09_kink_zoom.png")
GLMakie.save(path3, fig3; px_per_unit=4)
println("  saved $path3")
