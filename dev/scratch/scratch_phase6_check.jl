# scratch_phase6_check.jl
#
# Phase 6 (Visuals.jl) sanity check. Per the Phase 6 prompt's own
# testing-plan section: proportionally lighter than Phases 2-5 (no
# tracking, no tolerances, no correctness-critical math here -- this
# renders already-validated Phase 1-5 data). No @inferred section: a
# Makie.Figure's return type isn't a numerically meaningful thing to
# assert type-stability on, unlike Phases 2-5's actual tracked/computed
# data.
#
# Reuses the sphere/ellipsoid/Taubin heart fixtures already validated in
# scratch_phase5_check.jl / scratch_phase5_taubin_check.jl -- no new
# numerical fixtures needed, since Phase 6 doesn't compute anything new.

using Pkg
Pkg.activate(".")

include(joinpath(pwd(), "src", "HomotopyGetsReal.jl"))
using .HomotopyGetsReal
using .HomotopyGetsReal: Critical, Boundary, Singular, Artificial
using HomotopyContinuation
using GeometryBasics
using GLMakie
using Test

mkpath("phase6_renders")

println("=" ^ 70)
println("0. GLMakie headless-render capability (already confirmed empirically")
println("   before implementation; re-confirmed here as part of the regression)")
println("=" ^ 70)
let
    fig = Figure()
    ax = Axis(fig[1, 1])
    scatter!(ax, [0.0], [0.0])
    GLMakie.save("phase6_renders/00_smoke_test.png", fig)
    println("  OK: Figure/Axis/scatter!/save round-trip succeeded.")
end

println()
println("=" ^ 70)
println("1. plot_curve_decomposition: unit sphere's equatorial slice at z=0")
println("=" ^ 70)

@var x y z
f_sphere = x^2 + y^2 + z^2 - 1
F_sphere = System([f_sphere], variables = [x, y, z])
cfg = HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8)

vertices2d, edges2d = slice_at_z(F_sphere, 0.0, cfg)
println("  slice_at_z(F_sphere, 0.0, cfg) -> $(length(vertices2d)) vertices, $(length(edges2d)) edges")

fig1 = plot_curve_decomposition(vertices2d, edges2d; cfg = cfg, show_labels = true)
@test fig1 isa Makie.Figure
GLMakie.save("phase6_renders/01_sphere_equator_curve.png", fig1)
println("  plot_curve_decomposition (cfg-bounded, labeled) -> Figure, saved.")

fig1b = plot_curve_decomposition(vertices2d, edges2d; edge_color_by = :mono, show_vertices = false)
@test fig1b isa Makie.Figure
GLMakie.save("phase6_renders/01b_sphere_equator_curve_mono_noverts.png", fig1b)
println("  plot_curve_decomposition (edge_color_by=:mono, show_vertices=false) -> Figure, saved.")

println()
println("=" ^ 70)
println("2. plot_surface_decomposition(mesh): unit sphere, full pipeline")
println("=" ^ 70)

all_vertices, all_edges, all_faces, mesh = decompose_3d_surface(F_sphere, cfg)
println("  decompose_3d_surface -> $(length(all_vertices)) vertices, $(length(all_faces)) faces, ",
        "$(length(GeometryBasics.coordinates(mesh))) mesh vertices, $(length(GeometryBasics.faces(mesh))) triangles")

# color_by=:z has a meaningful range (-1 to 1) -- must NOT trip the new
# near-constant-colorrange warning (checked before radial_fn's near-constant
# case below ever has a chance to latch the one-shot Ref).
z_warns, fig2 = Test.collect_test_logs() do
    plot_surface_decomposition(mesh; color_by = :z, show_wireframe = true, cfg = cfg, vertices = all_vertices)
end
@test fig2 isa Makie.Figure
@test isempty([l for l in z_warns if l.level == Base.CoreLogging.Warn])
GLMakie.save("phase6_renders/02_sphere_mesh_by_z.png", fig2)
println("  plot_surface_decomposition(mesh; color_by=:z, wireframe, vertex overlay) -> Figure, saved, ",
        "0 near-constant-color warnings (meaningfully-varying range).")

# radial_fn is ~1.0 everywhere on a unit sphere (up to Float32 round-off) --
# must trip the near-constant-colorrange gate exactly once (first call ever
# to hit it in this process).
radial_fn(px, py, pz) = sqrt(px^2 + py^2 + pz^2)
radial_warns, fig2b = Test.collect_test_logs() do
    plot_surface_decomposition(mesh; color_by = radial_fn, show_colorbar = true)
end
@test fig2b isa Makie.Figure
radial_warn_logs = [l for l in radial_warns if l.level == Base.CoreLogging.Warn]
@test length(radial_warn_logs) == 1
@test occursin("near-constant", radial_warn_logs[1].message)
GLMakie.save("phase6_renders/02b_sphere_mesh_by_function.png", fig2b)
println("  plot_surface_decomposition(mesh; color_by=Function) on a near-constant function -> Figure, ",
        "saved, exactly 1 near-constant-color warning fired.")

# Second near-constant call -- must NOT warn again (one-shot latch).
radial_warns2, fig2b2 = Test.collect_test_logs() do
    plot_surface_decomposition(mesh; color_by = radial_fn)
end
@test isempty([l for l in radial_warns2 if l.level == Base.CoreLogging.Warn])
println("  second near-constant color_by call: 0 warnings (one-shot latch confirmed).")

@test_throws ArgumentError plot_surface_decomposition(mesh; color_by = :bogus)
println("  plot_surface_decomposition(mesh; color_by=:bogus) correctly throws ArgumentError.")

println()
println("=" ^ 70)
println("3. plot_surface_decomposition(faces): one-shot winding warning + per-cell coloring")
println("=" ^ 70)

# First call in this session -- must warn exactly once.
local_warns = Test.collect_test_logs() do
    global fig3 = plot_surface_decomposition(all_faces; show_wireframe = true, cfg = cfg)
end
warn_logs1 = [l for l in local_warns[1] if l.level == Base.CoreLogging.Warn]
@test length(warn_logs1) == 1
@test occursin("winding correction", warn_logs1[1].message)
@test fig3 isa Makie.Figure
GLMakie.save("phase6_renders/03_sphere_faces_by_cell.png", fig3)
println("  first call: exactly 1 warning emitted, message mentions winding correction. Figure saved.")

# Second call -- must NOT warn again.
local_warns2 = Test.collect_test_logs() do
    global fig3b = plot_surface_decomposition(all_faces)
end
warn_logs2 = [l for l in local_warns2[1] if l.level == Base.CoreLogging.Warn]
@test isempty(warn_logs2)
println("  second call: 0 warnings emitted (one-shot latch confirmed).")

println()
println("=" ^ 70)
println("4. interactive_3d_viewer: forces display, reuses plot_surface_decomposition(mesh)")
println("=" ^ 70)

fig4 = interactive_3d_viewer(mesh; color_by = :z, show_wireframe = false, cfg = cfg)
@test fig4 isa Makie.Figure
GLMakie.save("phase6_renders/04_sphere_interactive.png", fig4)
println("  interactive_3d_viewer -> Figure, saved (display() call did not error headlessly).")

println()
println("=" ^ 70)
println("5. Asymmetric ellipsoid (catches x/y/z mixups a symmetric sphere can't)")
println("=" ^ 70)

@var xe ye ze
f_ell = xe^2 + 4 * ye^2 + 9 * ze^2 - 1
F_ell = System([f_ell], variables = [xe, ye, ze])
ev, ee, ef, emesh = decompose_3d_surface(F_ell, cfg)
fig5 = plot_surface_decomposition(emesh; color_by = :z, cfg = cfg)
@test fig5 isa Makie.Figure
GLMakie.save("phase6_renders/05_ellipsoid_mesh.png", fig5)
println("  ellipsoid: decompose_3d_surface -> $(length(ef)) faces; plot_surface_decomposition -> Figure, saved.")

println()
println("=" ^ 70)
println("6. Taubin heart: multi-branch curve (edge_color_by=:cell distinguishes branches)")
println("   and full welded mesh (the most topologically complex Phase 5 test case)")
println("=" ^ 70)

@var xh yh zh
f_heart = (xh^2 + (1.2 * yh)^2 + zh^2 - 1)^3 - xh^2 * zh^3 - 0.1 * (1.2 * yh)^2 * zh^3
F_heart = System([f_heart], variables = [xh, yh, zh])
cfg_heart = HomotopyConfig{Float64}(
    bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3),
    edge_sample_density = 8, midslice_sample_density = 8,
)

# the narrow, multi-branch slab from the Phase 5 investigation, derived
# live from compute_critical_z_slices (not hardcoded) as the midpoint of
# the two closest consecutive critical z-values -- the same "narrowest
# slab" that stressed decompose_3d_surface's re-anchoring in Phase 5.
z_crits_heart = compute_critical_z_slices(F_heart, cfg_heart)
sort!(z_crits_heart)
gaps = diff(z_crits_heart)
narrow_i = argmin(gaps)
z_narrow_mid = (z_crits_heart[narrow_i] + z_crits_heart[narrow_i + 1]) / 2
v_narrow, e_narrow = slice_at_z(F_heart, z_narrow_mid, cfg_heart)
println("  narrowest slab [$(z_crits_heart[narrow_i]), $(z_crits_heart[narrow_i + 1])], ",
        "z_mid=$(z_narrow_mid): $(length(e_narrow)) edges")
fig6 = plot_curve_decomposition(v_narrow, e_narrow; cfg = cfg_heart, edge_color_by = :cell)
@test fig6 isa Makie.Figure
GLMakie.save("phase6_renders/06_heart_narrow_slab_curve.png", fig6)
println("  plot_curve_decomposition (8-edge multi-branch slab, edge_color_by=:cell) -> Figure, saved.")

t_heart = @elapsed (hv, he, hf, hmesh) = decompose_3d_surface(F_heart, cfg_heart)
println("  decompose_3d_surface (Taubin heart): $(length(hf)) faces, ",
        "$(length(GeometryBasics.coordinates(hmesh))) mesh vertices, ",
        "$(length(GeometryBasics.faces(hmesh))) triangles, elapsed=$(round(t_heart; digits=2))s")
fig6b = plot_surface_decomposition(hmesh; color_by = :z, show_wireframe = false, cfg = cfg_heart)
@test fig6b isa Makie.Figure
GLMakie.save("phase6_renders/06b_heart_mesh_by_z.png", fig6b)
println("  plot_surface_decomposition(mesh; color_by=:z) on full Taubin heart -> Figure, saved.")

fig6c = plot_surface_decomposition(hf; cfg = cfg_heart)  # already warned once above; expect none here
@test fig6c isa Makie.Figure
GLMakie.save("phase6_renders/06c_heart_faces_by_cell.png", fig6c)
println("  plot_surface_decomposition(faces) on full Taubin heart (14 faces, by-cell) -> Figure, saved.")

println()
println("=" ^ 70)
println("All Phase 6 sanity checks PASSED. Rendered PNGs saved under phase6_renders/")
println("for manual visual review (not a substitute for the automated checks above,")
println("a supplement -- see the Phase 6 prompt's testing-plan section).")
println("=" ^ 70)
