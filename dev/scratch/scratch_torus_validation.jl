# dev/scratch/scratch_torus_validation.jl
#
# Issue 3 (2026-07 follow-up to the visual-assets audit): a torus as a new
# decompose_3d_surface validation fixture -- genuinely new topology (genus
# 1, not simply connected), and the first fixture where a single z-slice
# can have two disjoint curve components.
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
println("SUMMARY")
println("="^70)
println("  compute_critical_z_slices matches closed-form |z|=1 prediction: ", matches_prediction)
println("  z=0 slice found 2 disjoint components (edges): ", n_components == 2)
println("  decompose_3d_surface (no incidence) succeeded: true, $(length(v1)) vertices / $(length(f1)) faces")
println("  decompose_3d_surface (incidence=true) succeeded: true, naked edges=$(length(naked))")
println("  max |f| residual over welded mesh: ", maximum(resids))

OUTDIR = joinpath(@__DIR__, "renders_torus_validation")
mkpath(OUTDIR)
GLMakie.activate!()
fig = plot_surface_decomposition(mesh2; color_by = :z, cfg = cfg, vertices = v2)
GLMakie.save(joinpath(OUTDIR, "torus.png"), fig; px_per_unit = 3)
println("\n  wrote ", joinpath(OUTDIR, "torus.png"))
