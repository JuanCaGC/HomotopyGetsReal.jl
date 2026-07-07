# scratch_phase6_investigate.jl
#
# Follow-up investigation requested after reviewing phase6_renders/:
#   1. Is 02b's patchy dark colorbar-collapsed-to-1.0 look just
#      auto-scaling noise on a near-constant radial_fn, or a mesh defect?
#   2. Does interactive_3d_viewer (04) show the same dark-patch artifact
#      as the KNOWN-uncorrected faces path (03), and if so, is that a
#      call-path bug (wrong method internally) or a genuine winding
#      defect in weld_mesh's "corrected" mesh that 02's z-coloring
#      (plus wireframe overlay) happened to visually mask?
#   3. Re-render sphere mesh/faces/interactive at realistic resolution
#      (cfg defaults: edge_sample_density=50, midslice_sample_density=100)
#      to separate "low-res faceting" from the patchy-dark-region issue.

using Pkg
Pkg.activate(".")

include(joinpath(pwd(), "src", "HomotopyGetsReal.jl"))
using .HomotopyGetsReal
using HomotopyContinuation
using GeometryBasics
using GLMakie
using LinearAlgebra
using Statistics
using InteractiveUtils

mkpath("phase6_renders")

@var x y z
f_sphere = x^2 + y^2 + z^2 - 1
F_sphere = System([f_sphere], variables = [x, y, z])
cfg = HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8)

println("=" ^ 70)
println("ISSUE 1: radial_fn range on mesh coordinates (low-res cfg, same as 02b)")
println("=" ^ 70)

_, _, faces, mesh = decompose_3d_surface(F_sphere, cfg)
pts = GeometryBasics.coordinates(mesh)
radial_fn(px, py, pz) = sqrt(px^2 + py^2 + pz^2)
radvals = [radial_fn(p[1], p[2], p[3]) for p in pts]
println("  n points = ", length(radvals))
println("  min(radial_fn) = ", minimum(radvals))
println("  max(radial_fn) = ", maximum(radvals))
println("  range (max-min) = ", maximum(radvals) - minimum(radvals))
println("  mean = ", mean(radvals), "  std = ", std(radvals))
# Show what plain Float32 mesh-point round-trip alone contributes, isolated
# from any residual solver error: points are stored as Point3f (Float32).
f32_only = [radial_fn(Float32(p[1]), Float32(p[2]), Float32(p[3])) for p in pts]
println("  (for comparison) range if computed in Float32 from otherwise-exact unit-sphere points: ",
        "not meaningful here since points already come from the solve, skip")

println()
println("=" ^ 70)
println("ISSUE 2a: is interactive_3d_viewer calling plot_surface_decomposition correctly?")
println("=" ^ 70)
println("Source of interactive_3d_viewer:")
println(@which interactive_3d_viewer(mesh))
m = methods(interactive_3d_viewer)
for mm in m
    println("  ", mm)
end

println()
println("=" ^ 70)
println("ISSUE 2b: direct comparison -- plot_surface_decomposition(mesh; show_wireframe=false)")
println("          called directly (NOT via interactive_3d_viewer) vs via interactive_3d_viewer")
println("=" ^ 70)

fig_direct = plot_surface_decomposition(mesh; color_by = :z, show_wireframe = false, cfg = cfg)
GLMakie.save("phase6_renders/investigate_direct_nowireframe.png", fig_direct)
println("  Saved investigate_direct_nowireframe.png (direct call, no interactive_3d_viewer).")

fig_viewer = interactive_3d_viewer(mesh; color_by = :z, show_wireframe = false, cfg = cfg)
GLMakie.save("phase6_renders/investigate_via_viewer.png", fig_viewer)
println("  Saved investigate_via_viewer.png (via interactive_3d_viewer).")

println()
println("=" ^ 70)
println("ISSUE 2c: genuine winding check on weld_mesh's output -- for EVERY triangle,")
println("          is its normal (from vertex winding) pointing outward (same sign as")
println("          the position vector at its centroid, since a sphere centered at 0")
println("          has outward normal == outward radial direction)?")
println("=" ^ 70)

mesh_pts = GeometryBasics.coordinates(mesh)
mesh_tris = GeometryBasics.faces(mesh)
n_inward = 0
n_outward = 0
n_degenerate = 0
worst_dot = 1.0
for t in mesh_tris
    i1, i2, i3 = Int(t[1]), Int(t[2]), Int(t[3])
    p1, p2, p3 = mesh_pts[i1], mesh_pts[i2], mesh_pts[i3]
    v1 = p2 .- p1
    v2 = p3 .- p1
    n = cross(Vector(v1), Vector(v2))
    nnorm = norm(n)
    centroid = (Vector(p1) .+ Vector(p2) .+ Vector(p3)) ./ 3
    if nnorm < 1e-12
        global n_degenerate += 1
        continue
    end
    n_hat = n ./ nnorm
    radial_hat = centroid ./ norm(centroid)
    d = dot(n_hat, radial_hat)
    global worst_dot = min(worst_dot, d)
    if d > 0
        global n_outward += 1
    else
        global n_inward += 1
    end
end
println("  total triangles = ", length(mesh_tris))
println("  outward-facing (dot(normal, radial) > 0) = ", n_outward)
println("  INWARD-facing  (dot(normal, radial) < 0) = ", n_inward)
println("  degenerate (zero-area) = ", n_degenerate)
println("  worst (most negative or smallest) dot(normal, radial) = ", worst_dot)

println()
println("=" ^ 70)
println("ISSUE 2d: does the SAME inward-normal defect exist in 02's wireframe+vertex")
println("          render? Re-render it explicitly to check whether wireframe visually")
println("          masks the same patches (same mesh, so same defect if any is present).")
println("=" ^ 70)
_, _, _, mesh2 = decompose_3d_surface(F_sphere, cfg)  # re-derive, should be identical
println("  mesh === mesh2 (same object)? ", mesh === mesh2, "  (expect false, re-solved; check below if numerically identical)")
pts1 = GeometryBasics.coordinates(mesh)
pts2 = GeometryBasics.coordinates(mesh2)
println("  same number of points: ", length(pts1) == length(pts2))
if length(pts1) == length(pts2)
    maxdiff = maximum(norm(Vector(pts1[i]) .- Vector(pts2[i])) for i in 1:length(pts1))
    println("  max per-point difference across two independent decompose_3d_surface calls: ", maxdiff)
end

println()
println("=" ^ 70)
println("3. Re-render sphere at realistic resolution (edge_sample_density=40, midslice_sample_density=60)")
println("   to separate low-res faceting from the patchy-dark-region issue.")
println("=" ^ 70)

cfg_hi = HomotopyConfig{Float64}(edge_sample_density = 40, midslice_sample_density = 60)
t_hi = @elapsed (hi_vertices, hi_edges, hi_faces, hi_mesh) = decompose_3d_surface(F_sphere, cfg_hi)
println("  decompose_3d_surface (hi-res): $(length(hi_faces)) faces, ",
        "$(length(GeometryBasics.coordinates(hi_mesh))) mesh vertices, ",
        "$(length(GeometryBasics.faces(hi_mesh))) triangles, elapsed=$(round(t_hi; digits=2))s")

fig_hi_mesh = plot_surface_decomposition(hi_mesh; color_by = :z, show_wireframe = false, cfg = cfg_hi)
GLMakie.save("phase6_renders/hires_02_sphere_mesh_by_z.png", fig_hi_mesh)
println("  Saved hires_02_sphere_mesh_by_z.png")

fig_hi_faces = plot_surface_decomposition(hi_faces; cfg = cfg_hi)
GLMakie.save("phase6_renders/hires_03_sphere_faces_by_cell.png", fig_hi_faces)
println("  Saved hires_03_sphere_faces_by_cell.png")

fig_hi_viewer = interactive_3d_viewer(hi_mesh; color_by = :z, show_wireframe = false, cfg = cfg_hi)
GLMakie.save("phase6_renders/hires_04_sphere_interactive.png", fig_hi_viewer)
println("  Saved hires_04_sphere_interactive.png")

# Also redo the winding check at hi-res to see if it's resolution-dependent.
hi_pts = GeometryBasics.coordinates(hi_mesh)
hi_tris = GeometryBasics.faces(hi_mesh)
n_inward_hi = 0
n_outward_hi = 0
for t in hi_tris
    i1, i2, i3 = Int(t[1]), Int(t[2]), Int(t[3])
    p1, p2, p3 = hi_pts[i1], hi_pts[i2], hi_pts[i3]
    v1 = p2 .- p1
    v2 = p3 .- p1
    n = cross(Vector(v1), Vector(v2))
    nnorm = norm(n)
    nnorm < 1e-12 && continue
    n_hat = n ./ nnorm
    centroid = (Vector(p1) .+ Vector(p2) .+ Vector(p3)) ./ 3
    radial_hat = centroid ./ norm(centroid)
    d = dot(n_hat, radial_hat)
    if d > 0
        global n_outward_hi += 1
    else
        global n_inward_hi += 1
    end
end
println("  hi-res: outward = $(n_outward_hi), inward = $(n_inward_hi), total = $(length(hi_tris))")

println()
println("Investigation complete.")
