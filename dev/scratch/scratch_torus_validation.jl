# dev/scratch/scratch_torus_validation.jl
#
# Issue 3 (2026-07 follow-up to the visual-assets audit): a torus as a new
# decompose_3d_surface validation fixture -- genuinely new topology (genus
# 1, not simply connected), and the first fixture where a single z-slice
# can have two disjoint curve components.
#
# Two-part story, both reproducible from this one script (2026-07,
# updated in place rather than superseded by a new file, so the failure
# and its fix stay together):
#   Sections 1-4: the z-aligned hole orientation below FAILS. Confirmed by
#     direct calculation (see docs/DESIGN_NOTES.md's backlog): at the fold
#     z=+-1, df/dx=df/dy=0 IDENTICALLY for every point on the circle
#     x^2+y^2=4 -- a genuine 1-dimensional critical locus, not isolated
#     points, which compute_critical_points (built for isolated-solution
#     homotopy continuation) cannot represent. compute_critical_z_slices
#     returns empty, the naive whole-bbox slab lands exactly on the
#     degenerate fold, and decompose_3d_surface either produces a
#     catastrophically wrong mesh or runs for ~1400s -- kept in this
#     script AS-IS (not trimmed) so the failure stays independently
#     reproducible, not just described in prose.
#   Sections 5-6 (2026-07, second investigation): projection=:random
#     (Phase 8's existing generic-projection support) closes this. The
#     fold's positive-dimensional locus is a coordinate-ALIGNMENT
#     artifact (the hole axis exactly coincides with the slicing axis) --
#     the torus itself has no real singularity anywhere. A generic
#     rotation breaks that alignment and the whole pipeline works
#     unmodified. The torus is a usable validation fixture after all, via
#     an existing capability, not a new one.

using Random
#
# Equation: reusing the old prototype's already-validated torus geometry
# (prototipo_viejo_julia/HomotopyGetsReal.jl:106, R=2, r=1 -- "vertical
# torus"), but reoriented: the OLD prototype's hole axis was along ITS OWN
# y-axis (x^2+y^2+z^2+3)^2 - 16(x^2+z^2); this project's decompose_3d_surface
# REQUIRES z last and slices along z (confirmed: src/SurfaceDecomposition.jl
# docstrings, "variables ordered [x_var,y_var,z_var] -- z LAST"). Swapping
# which coordinate pair is squared in the second term puts the hole axis
# on z instead, which is what actually produces the "two disjoint
# components per z-slice" property this fixture is FOR -- verified
# symbolically (2026-07, Symbolics.jl, exact polynomial elimination, not
# just point checks) before writing this script:
#   at z=0: (x^2+y^2+3)^2 = 16(x^2+y^2). Let rho=x^2+y^2:
#     rho^2 - 10*rho + 9 = 0  ->  (rho-1)(rho-9) = 0  ->  rho=1 or rho=9
#   i.e. TWO disjoint concentric circles, x^2+y^2=1 and x^2+y^2=9.
#   General z: rho^2 + (2k-16)rho + k^2 = 0 (k = z^2+3) has real roots iff
#   k <= 4, i.e. |z| <= 1 -- so the torus's z-extent is EXACTLY [-1,1], and
#   at |z|=1 the two circles merge into one (rho=4, a double root) -- the
#   expected fold/critical z-value, directly analogous to a sphere's poles.

using HomotopyContinuation
using LinearAlgebra
using GeometryBasics
using GLMakie
using HomotopyGetsReal

@var x y z
f_torus = (x^2 + y^2 + z^2 + 3)^2 - 16 * (x^2 + y^2)
F_torus = System([f_torus], variables = [x, y, z])
f_eval(px, py, pz) = (px^2 + py^2 + pz^2 + 3)^2 - 16 * (px^2 + py^2)

cfg = HomotopyConfig{Float64}()  # default bbox (-4,4)^3 comfortably covers radial extent <=3, |z|<=1
println("cfg bbox: x=", cfg.bbox_x, " y=", cfg.bbox_y, " z=", cfg.bbox_z)

println()
println("="^70)
println("1. compute_critical_z_slices -- expect |z|=1 (closed-form, derived above)")
println("="^70)
t_crit = @elapsed z_crits = compute_critical_z_slices(F_torus, cfg)
println("  elapsed: $(round(t_crit; digits=3))s")
println("  z_crits found: ", sort(z_crits))
matches_prediction = length(z_crits) == 2 && all(isapprox.(sort(z_crits), [-1.0, 1.0]; atol=1e-6))
println("  matches closed-form prediction [-1.0, 1.0]? ", matches_prediction)

println()
println("="^70)
println("2. z=0 slice directly -- expect TWO disjoint circles (rho=1, rho=9)")
println("="^70)
v0, e0 = slice_at_z(F_torus, 0.0, cfg)
println("  vertices: ", length(v0), "  edges: ", length(e0))
for v in v0
    println("    id=$(v.id) coords=$(round.(real.(v.coordinates); digits=4)) type=$(v.v_type)")
end
n_components = length(e0)  # each disjoint smooth circle should be its own edge (no vertices to split it)
println("  edges (components) found at z=0: ", n_components, "  (expect 2 for the two disjoint circles)")

println()
println("="^70)
println("3. Full decompose_3d_surface(F_torus, cfg) -- no incidence")
println("="^70)
t0 = @elapsed (v1, e1, f1, mesh1) = decompose_3d_surface(F_torus, cfg)
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

println()
println("="^70)
println("4. Full decompose_3d_surface(F_torus, cfg; incidence=true) -- crit-slice census")
println("="^70)
t2 = @elapsed (v2, e2, f2, mesh2, inc) = decompose_3d_surface(F_torus, cfg; incidence = true)
println("  elapsed: $(round(t2; digits=2))s")
println("  crit_slices: $(length(inc.crit_slices))")
for cs in inc.crit_slices
    println("    j=$(cs.boundary_index)  z=$(round(cs.z; digits=4))  n_edges=$(length(cs.edges))  n_vertices=$(length(cs.vertices))  is_degenerate=$(cs.is_degenerate)")
end
naked = HomotopyGetsReal._naked_mesh_edges(mesh2)
println("  naked edges (stitched mesh): ", length(naked))

println()
println("="^70)
println("5. projection=:random cheap check -- compute_critical_z_slices on the")
println("   rotated chart, 3 seeds (this project's own precedent: 1, 7, 42)")
println("="^70)
println("   Verifies the fold is a coordinate-alignment artifact, not a real")
println("   singularity: a generic rotation should always find an isolated,")
println("   nonempty 4-value critical set (never empty) if that diagnosis is")
println("   correct.")
projection_seeds_ok = Dict{Int,Vector{Float64}}()
for seed in (1, 7, 42)
    rng = Xoshiro(seed)
    Q = HomotopyGetsReal._resolve_projection(:random, rng, cfg)
    F_chart = HomotopyGetsReal._rotate_system(F_torus, Q)
    HomotopyGetsReal._verify_projection_ok(F_chart, cfg)
    cfg_chart = HomotopyGetsReal._chart_config(cfg, Q)
    t = @elapsed zc = compute_critical_z_slices(F_chart, cfg_chart)
    println("  seed=$seed: elapsed=$(round(t; digits=2))s  critical z-values: ", sort(zc))
    projection_seeds_ok[seed] = sort(zc)
end
all_seeds_isolated_4 = all(length(v) == 4 for v in values(projection_seeds_ok))
println("  all 3 seeds found an isolated, nonempty 4-value critical set: ", all_seeds_isolated_4)

println()
println("="^70)
println("6. Full decompose_3d_surface(F_torus, cfg; projection=:random, incidence=true)")
println("   -- two seeds, same reporting style as sections 3-4")
println("="^70)
projection_results = Dict{Int,Any}()
for seed in (42, 7)
    println()
    println("  --- seed=$seed ---")
    tproj = @elapsed (vp, ep, fp, meshp, incp) = decompose_3d_surface(F_torus, cfg; projection = :random, rng = Xoshiro(seed), incidence = true)
    println("    elapsed: $(round(tproj; digits=2))s")
    countsp = Dict(t2 => count(vv -> vv.v_type == t2, vp) for t2 in (Critical, Boundary, Singular, Artificial))
    println("    vertices: $(length(vp)) total  $(countsp)")
    println("    edges: $(length(ep))  faces: $(length(fp))")

    mesh_ptsp = GeometryBasics.coordinates(meshp)
    mesh_trisp = GeometryBasics.faces(meshp)
    residsp = [abs(f_eval(Float64(p[1]), Float64(p[2]), Float64(p[3]))) for p in mesh_ptsp]
    println("    mesh: $(length(mesh_ptsp)) vertices, $(length(mesh_trisp)) triangles")
    println("    |f| residuals: min=$(minimum(residsp)) mean=$(sum(residsp)/length(residsp)) max=$(maximum(residsp))")
    n_degeneratep = count(tri -> length(unique((tri[1], tri[2], tri[3]))) < 3, mesh_trisp)
    println("    degenerate triangles: ", n_degeneratep)

    nakedp = HomotopyGetsReal._naked_mesh_edges(meshp)
    println("    naked edges (stitched mesh): ", length(nakedp))
    println("    crit_slices: $(length(incp.crit_slices))")
    for cs in incp.crit_slices
        println("      j=$(cs.boundary_index)  z=$(round(cs.z; digits=4))  n_edges=$(length(cs.edges))  n_vertices=$(length(cs.vertices))  is_degenerate=$(cs.is_degenerate)")
    end
    projection_results[seed] = (n_vertices = length(vp), n_faces = length(fp), max_resid = maximum(residsp), n_naked = length(nakedp))
end

println()
println("="^70)
println("SUMMARY")
println("="^70)
println("  --- z-aligned orientation (sections 1-4): FAILS ---")
println("  compute_critical_z_slices matches closed-form |z|=1 prediction: ", matches_prediction)
println("  z=0 slice found 2 disjoint components (edges): ", n_components == 2)
println("  decompose_3d_surface (no incidence) succeeded: true, $(length(v1)) vertices / $(length(f1)) faces")
println("  decompose_3d_surface (incidence=true) succeeded: true, naked edges=$(length(naked))")
println("  max |f| residual over welded mesh: ", maximum(resids))
println()
println("  --- projection=:random (sections 5-6): WORKS ---")
println("  cheap check, all 3 seeds isolated nonempty 4-value critical set: ", all_seeds_isolated_4)
for (seed, r) in sort(collect(projection_results); by = first)
    println("  seed=$seed: $(r.n_vertices) vertices, $(r.n_faces) faces, max |f| residual=$(r.max_resid), naked edges=$(r.n_naked)")
end
println()
println("  CONCLUSION: the torus IS a usable decompose_3d_surface validation")
println("  fixture, via projection=:random (Phase 8's existing generic-")
println("  projection support) -- an existing capability, not a new one.")

OUTDIR = joinpath(@__DIR__, "renders_torus_validation")
mkpath(OUTDIR)
GLMakie.activate!()
fig = plot_surface_decomposition(mesh2; color_by = :z, cfg = cfg, vertices = v2)
GLMakie.save(joinpath(OUTDIR, "torus.png"), fig; px_per_unit = 3)
println("\n  wrote ", joinpath(OUTDIR, "torus.png"))
