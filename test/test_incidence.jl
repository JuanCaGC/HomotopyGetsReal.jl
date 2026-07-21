@testset "Incidence (Phase 9a)" begin
    using LinearAlgebra
    using Random

println("=" ^ 70)
println("1. sphere: 5-tuple arity, degenerate crit-slices, fold anchors at the poles")
println("=" ^ 70)

@var x y z
F_sph = System([x^2 + y^2 + z^2 - 1], variables = [x, y, z])
cfg = HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8)

res5 = decompose_3d_surface(F_sph, cfg; incidence = true)
@test length(res5) == 5
vs, es, fs, ms, inc = res5
@test inc isa SurfaceIncidence{Float64}

# z_bounds = [-4, -1, 1, 4]: two interior boundaries (the poles), both of which
# decompose to EMPTY slices (measured: fold-type point-slices fail the imaginary
# filter) -> is_degenerate, no edges, fold anchors do the work.
@test length(inc.crit_slices) == 2
@test all(cs -> cs.is_degenerate, inc.crit_slices)
@test all(cs -> isempty(cs.edges), inc.crit_slices)
@test isempty(inc.edge_faces)

for f in fs
    lb = inc.column_landings_bottom[f.id]
    lt = inc.column_landings_top[f.id]
    @test all(l -> l.kind == :critical_point, lb)  # south pole anchor
    @test all(l -> l.kind == :critical_point, lt)  # north pole anchor
    # mechanism check on a KNOWN fixture (not a 9b threshold): columns converge
    # to the poles within ~1e-2 at these densities, so the assigned distance to
    # the pole critical vertex must be small.
    @test all(l -> l.dist < 0.05, lb)
    @test all(l -> l.dist < 0.05, lt)
    @test inc.face_bottom_anchor[f.id] != 0
    @test inc.face_top_anchor[f.id] != 0
    @test isempty(inc.face_bottom_edges[f.id])
    @test isempty(inc.face_top_edges[f.id])
end
println("  sphere: both pole boundaries degenerate; every column of every face anchored to a pole")
println("  critical vertex (max assigned dist over all faces/columns: ",
        maximum(l.dist for d in (inc.column_landings_bottom, inc.column_landings_top) for ls in values(d) for l in ls), ")")

res4 = decompose_3d_surface(F_sph, cfg)
@test length(res4) == 4
println("  incidence = false still returns the 4-tuple.")

println()
println("=" ^ 70)
println("2. ellipsoid: same structure, asymmetric control")
println("=" ^ 70)

@var xe ye ze
F_ell = System([xe^2 + 4 * ye^2 + 9 * ze^2 - 1], variables = [xe, ye, ze])
_, _, fe, _, ince = decompose_3d_surface(F_ell, cfg; incidence = true)
@test length(ince.crit_slices) == 2
@test all(cs -> cs.is_degenerate, ince.crit_slices)
@test all(f -> ince.face_bottom_anchor[f.id] != 0 && ince.face_top_anchor[f.id] != 0, fe)
println("  ellipsoid: both fold boundaries anchored (dist max: ",
        maximum(l.dist for d in (ince.column_landings_bottom, ince.column_landings_top) for ls in values(d) for l in ls), ")")

println()
println("=" ^ 70)
println("3. z - x^2 under projection = :random + incidence: bbox-only boundaries, Phase 8 interplay")
println("=" ^ 70)

F_par = System([z - x^2], variables = [x, y, z])
vp, ep, fp, mp, incp = decompose_3d_surface(F_par, cfg; projection = :random, rng = Xoshiro(42), incidence = true)
@test incp isa SurfaceIncidence{Float64}
@test isempty(incp.crit_slices)   # no critical z at all -> no interior boundaries
all_landings = vcat(collect(values(incp.column_landings_bottom))..., collect(values(incp.column_landings_top))...)
@test !isempty(all_landings)
@test all(l -> l.kind == :bbox, all_landings)
println("  rotated parabolic sheet: no interior boundaries, all $(length(all_landings)) landings :bbox,")
println("  incidence world-mapped through the projection recursion without error.")

println()
println("=" ^ 70)
println("4. id-namespace uniqueness across main return and crit-slice cells")
println("=" ^ 70)

# Taubin-lite check is in the slow suite; here assert on the sphere run.
all_vids = vcat([v.id for v in vs], [v.id for cs in inc.crit_slices for v in cs.vertices])
all_eids = vcat([e.id for e in es], [e.id for cs in inc.crit_slices for e in cs.edges])
@test length(unique(all_vids)) == length(all_vids)
@test length(unique(all_eids)) == length(all_eids)
println("  no vertex/edge id collisions between slab cells and crit-slice cells.")

println()
println("=" ^ 70)
println("5. _naked_mesh_edges unit test on a hand-built open mesh")
println("=" ^ 70)

pts_unit = [GeometryBasics.Point3f(0, 0, 0), GeometryBasics.Point3f(1, 0, 0),
            GeometryBasics.Point3f(1, 1, 0), GeometryBasics.Point3f(0, 1, 0)]
tris_unit = [GeometryBasics.TriangleFace{Int}(1, 2, 3), GeometryBasics.TriangleFace{Int}(1, 3, 4)]
open_mesh = GeometryBasics.Mesh(pts_unit, tris_unit)
naked_unit = HomotopyGetsReal._naked_mesh_edges(open_mesh)
@test length(naked_unit) == 4                       # square perimeter
@test (1, 3) ∉ naked_unit                           # shared diagonal has 2 triangles
@test Set(naked_unit) == Set([(1, 2), (2, 3), (3, 4), (1, 4)])
println("  square-with-diagonal: 4 naked perimeter edges, shared diagonal excluded. ✓")

end
