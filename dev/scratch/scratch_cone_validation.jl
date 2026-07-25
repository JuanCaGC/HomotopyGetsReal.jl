# dev/scratch/scratch_cone_validation.jl
#
# Issue 3 follow-up (2026-07): trial run of a cone (x^2+y^2-z^2=0) as the
# "intermediate complexity between sphere and Taubin heart" fixture --
# proposed as a lower-risk alternative to the Whitney umbrella, which
# was confirmed (dev/scratch/scratch_torus_validation.jl's sibling
# investigation, same session) to crash decompose_3d_surface outright
# and separately hit the same reducible-curve critical-point-detection
# gap found on the node curve. The cone has exactly ONE isolated apex
# singularity at z=0 (not an entire singular line), and its z!=0 slices
# are ordinary circles -- the same degeneracy CLASS (single degenerate
# z, surrounded by well-conditioned slices) that Taubin's own [-1,1]
# slab already handles via _robust_slice_at_z's retry mechanism.
# NOT run before this script -- reporting exactly what happens, live.

using HomotopyContinuation
using LinearAlgebra
using GeometryBasics
using GLMakie
using HomotopyGetsReal

@var x y z
f_cone = x^2 + y^2 - z^2
F_cone = System([f_cone], variables = [x, y, z])
f_eval(px, py, pz) = px^2 + py^2 - pz^2

cfg = HomotopyConfig{Float64}(bbox_x = (-2.0, 2.0), bbox_y = (-2.0, 2.0), bbox_z = (-2.0, 2.0))
println("cfg bbox: x=", cfg.bbox_x, " y=", cfg.bbox_y, " z=", cfg.bbox_z)
println("(bbox chosen so the cone's radius stays comparable to its z-extent within the box,")
println(" unlike default (-4,4)^3 where an unbounded cone would run far past a useful window)")

println()
println("="^70)
println("1. compute_critical_z_slices -- the apex is a single isolated point, not a smooth fold")
println("="^70)
try
    t = @elapsed zc = compute_critical_z_slices(F_cone, cfg)
    println("  elapsed: $(round(t; digits=3))s")
    println("  z_crits found: ", sort(zc))
catch e
    println("  ERROR: ", sprint(showerror, e))
end

println()
println("="^70)
println("2. slice_at_z at a few representative z values (apex at z=0, ordinary circles elsewhere)")
println("="^70)
for zval in [-1.0, -0.001, 0.0, 0.001, 1.0]
    print("  z=$zval: ")
    try
        v3, e3 = slice_at_z(F_cone, zval, cfg)
        println("OK, ", length(v3), " vertices, ", length(e3), " edges")
    catch e
        println("ERROR: ", sprint(showerror, e))
    end
end

println()
println("="^70)
println("3. Full decompose_3d_surface(F_cone, cfg) -- no incidence")
println("="^70)
try
    t0 = @elapsed (v1, e1, f1, mesh1) = decompose_3d_surface(F_cone, cfg)
    println("  elapsed: $(round(t0; digits=2))s")
    counts1 = Dict(t2 => count(vv -> vv.v_type == t2, v1) for t2 in (Critical, Boundary, Singular, Artificial))
    println("  vertices: $(length(v1)) total  $(counts1)")
    println("  edges: $(length(e1))  faces: $(length(f1))")

    mesh_pts = GeometryBasics.coordinates(mesh1)
    mesh_tris = GeometryBasics.faces(mesh1)
    resids = [abs(f_eval(Float64(p[1]), Float64(p[2]), Float64(p[3]))) for p in mesh_pts]
    println("  mesh: $(length(mesh_pts)) vertices, $(length(mesh_tris)) triangles")
    println("  |f| residuals: min=$(minimum(resids)) mean=$(sum(resids)/length(resids)) max=$(maximum(resids))")
    n_degenerate = count(tri -> length(unique((tri[1], tri[2], tri[3]))) < 3, mesh_tris)
    println("  degenerate triangles: ", n_degenerate)

    global mesh_ok = true
    global v1g, mesh1g = v1, mesh1
catch e
    println("  ERROR: ", sprint(showerror, e))
    global mesh_ok = false
end

println()
println("="^70)
println("4. Full decompose_3d_surface(F_cone, cfg; incidence=true) -- crit-slice census")
println("="^70)
try
    t2 = @elapsed (v2, e2, f2, mesh2, inc) = decompose_3d_surface(F_cone, cfg; incidence = true)
    println("  elapsed: $(round(t2; digits=2))s")
    println("  crit_slices: $(length(inc.crit_slices))")
    for cs in inc.crit_slices
        println("    j=$(cs.boundary_index)  z=$(round(cs.z; digits=4))  n_edges=$(length(cs.edges))  n_vertices=$(length(cs.vertices))  is_degenerate=$(cs.is_degenerate)")
    end
    naked = HomotopyGetsReal._naked_mesh_edges(mesh2)
    println("  naked edges (stitched mesh): ", length(naked))
    global v2g, mesh2g = v2, mesh2
    global incidence_ok = true
catch e
    println("  ERROR: ", sprint(showerror, e))
    global incidence_ok = false
end

println()
println("="^70)
println("SUMMARY")
println("="^70)
println("  decompose_3d_surface (no incidence) succeeded: ", mesh_ok)
println("  decompose_3d_surface (incidence=true) succeeded: ", incidence_ok)

if mesh_ok
    OUTDIR = joinpath(@__DIR__, "renders_cone_validation")
    mkpath(OUTDIR)
    GLMakie.activate!()
    fig = plot_surface_decomposition(incidence_ok ? mesh2g : mesh1g; color_by = :z, cfg = cfg, vertices = incidence_ok ? v2g : v1g)
    GLMakie.save(joinpath(OUTDIR, "cone.png"), fig; px_per_unit = 3)
    println("\n  wrote ", joinpath(OUTDIR, "cone.png"))
end
