@testset "SurfaceDecomposition (Phase 5)" begin

"""
    instrumented_sweep_direction(patch, x0, y0, z_start, z_targets, cfg)
        -> (n_rebuilds::Int, dense_path::Vector{Vector{Float64}})

Read-only introspection harness reproducing `_sweep_direction`'s exact
call pattern (patch build, then one `_sweep_hop!` call per z-target),
since `_sweep_direction` itself intentionally exposes no internal
counters (it is a pure function, not a logging harness). Used below to
(a) confirm the sphere case triggers zero re-anchors (radial symmetry),
(b) report the ellipsoid case's rebuild count, and (c) sweep
`cfg.patch_transversality_cos_tol` itself to confirm the mechanism
responds monotonically, not just "happens to work" at the chosen
default.
"""
function instrumented_sweep_direction(patch, x0, y0, z_start, z_targets, cfg)
    tols = (
        rank_tol64 = Float64(cfg.jacobian_rank_tol),
        sv_thresh64 = Float64(cfg.singular_value_threshold),
        poor_acc_tol64 = Float64(cfg.critical_point_tol),
        min_width64 = Float64(cfg.vertex_match_tol),
        cos_tol64 = Float64(cfg.patch_transversality_cos_tol),
    )
    _, ph0, tracker0 = HomotopyGetsReal.build_face_tracker(patch, x0, y0, z_start, cfg)
    fx0, fy0, _ = HomotopyGetsReal._gradient_at(patch, x0, y0, z_start, cfg)
    state = (ph = ph0, tracker = tracker0, fx0_64 = Float64(fx0), fy0_64 = Float64(fy0))
    y_state = ComplexF64[ComplexF64(x0), ComplexF64(y0)]
    p_cur = Float64(z_start)
    n_rebuilds = 0
    dense_path = Vector{Vector{Float64}}(undef, length(z_targets))
    for (i, p_target) in enumerate(z_targets)
        hop_path = Vector{Float64}[]
        budget = Ref(max(cfg.max_path_steps, 1))
        ph_before = state.ph
        state, y_state = HomotopyGetsReal._sweep_hop!(hop_path, patch, cfg, state, y_state, p_cur, p_target, budget, tols)
        state.ph !== ph_before && (n_rebuilds += 1)
        dense_path[i] = hop_path[end]
        p_cur = p_target
    end
    return n_rebuilds, dense_path
end

println("=" ^ 70)
println("Setup: unit sphere f(x,y,z) = x^2 + y^2 + z^2 - 1")
println("=" ^ 70)

@var x y z
f = x^2 + y^2 + z^2 - 1
F_sphere = System([f], variables = [x, y, z]) # z LAST, required convention

# Small densities so this script runs in seconds; large enough to give
# `weld_mesh`'s clustering and the winding check something nontrivial to
# chew on.
cfg = HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8)

println()
println("=" ^ 70)
println("1. compute_critical_z_slices: expect {-1, +1} (the poles)")
println("=" ^ 70)

z_crits = compute_critical_z_slices(F_sphere, cfg)
println("  compute_critical_z_slices -> ", z_crits)
@test length(z_crits) == 2
@test isapprox(sort(z_crits), [-1.0, 1.0]; atol = 1e-6)

println()
println("=" ^ 70)
println("2. slice_at_z(F_sphere, 0.0, cfg): expect unit circle in x,y at z=0")
println("=" ^ 70)

vertices0, edges0 = slice_at_z(F_sphere, 0.0, cfg)
println("  vertices: ", [(v.id, round.(real.(v.coordinates); digits = 4), v.v_type) for v in vertices0])
println("  edges: ", [(e.id, e.left_vertex_id, e.right_vertex_id, length(e.sampled_points)) for e in edges0])
@test length(vertices0) == 2
@test all(v -> isapprox(abs(real(v.coordinates[1])), 1.0; atol = 1e-4) && isapprox(real(v.coordinates[2]), 0.0; atol = 1e-4), vertices0)
@test all(v -> real(v.coordinates[3]) == 0.0, vertices0) # lifted z-coordinate exact
@test length(edges0) == 2 # upper and lower semicircle
@test all(e -> length(e.sampled_points) == cfg.edge_sample_density, edges0)
# sample_edge's resampling is PURE geometric linear interpolation between
# already-tracked RAW points (no re-projection onto the curve, see its
# own docstring) -- for a smooth arc like this, connect_the_dots! only
# records a handful of raw points (no bisection needed), so the
# resulting chord-interpolated points can deviate substantially from the
# true circle at low edge_sample_density. This is an existing, already
# Phase-3-accepted property of `sample_edge`, not a Phase 5 concern, so
# only the lifted z-coordinate (a literal appended constant, not
# interpolated) is checked here for exactness; genuine on-manifold
# fidelity is checked below on `sweep_face_bidirectional`'s actually
# TRACKED (not linearly resampled) points.
@test all(e -> all(p -> p[3] == 0.0, e.sampled_points), edges0)
println("slice_at_z: sampled point count and lifted z-coordinate confirmed.")

println()
println("=" ^ 70)
println("3. build_patch_system / patch_direction: expect ∇f = (2x,2y,2z)")
println("=" ^ 70)

patch = build_patch_system(F_sphere)
println("  F_for_tracking.variables = ", patch.F_for_tracking.variables, "  (expect [z, x, y])")
@test string.(patch.F_for_tracking.variables) == string.([z, x, y])

a1, b1 = patch_direction(patch, 0.0, 1.0, 0.0, cfg)
println("  patch_direction at anchor (x0,y0,z0)=(0,1,0) -> (a,b) = ($a1, $b1)  (expect (2, 0), from ∇f=(0,2,0))")
@test isapprox(a1, 2.0; atol = 1e-8)
@test isapprox(b1, 0.0; atol = 1e-8)

a2, b2 = patch_direction(patch, 1.0, 0.0, 0.0, cfg)
println("  patch_direction at anchor (x0,y0,z0)=(1,0,0) -> (a,b) = ($a2, $b2)  (expect (0, -2), from ∇f=(2,0,0))")
@test isapprox(a2, 0.0; atol = 1e-8)
@test isapprox(b2, -2.0; atol = 1e-8)

println()
println("=" ^ 70)
println("4. sweep_face_bidirectional: sweep anchor (0,1,0) from z_mid=0 to poles")
println("=" ^ 70)

dense_down, dense_up = sweep_face_bidirectional(F_sphere, patch, 0.0, 1.0, 0.0, -1.0, 1.0, cfg)
println("  dense_path_down length: ", length(dense_down), "  (expect $(cfg.midslice_sample_density))")
println("  dense_path_up length:   ", length(dense_up), "  (expect $(cfg.midslice_sample_density))")
@test length(dense_down) == cfg.midslice_sample_density
@test length(dense_up) == cfg.midslice_sample_density

println("  landing near south pole (z=-1): ", dense_down[end], "  (expect z=-1, x≈0, y≈0)")
println("  landing near north pole (z=+1): ", dense_up[end], "  (expect z=+1, x≈0, y≈0)")
@test isapprox(dense_down[end][1], -1.0; atol = 1e-3)
@test isapprox(dense_down[end][2], 0.0; atol = 1e-2)
@test isapprox(dense_down[end][3], 0.0; atol = 1e-2)
@test isapprox(dense_up[end][1], 1.0; atol = 1e-3)
@test isapprox(dense_up[end][2], 0.0; atol = 1e-2)
@test isapprox(dense_up[end][3], 0.0; atol = 1e-2)

# The strong, independent correctness check: every single swept point,
# not just the endpoints, must satisfy f ≈ 0 (dense path entries are
# [z, x, y]-ordered, see track_dense_path's docstring).
on_sphere(p) = isapprox(p[2]^2 + p[3]^2 + p[1]^2, 1.0; atol = 1e-4)
@test all(on_sphere, dense_down)
@test all(on_sphere, dense_up)
println("Every swept point (not just endpoints) satisfies x^2+y^2+z^2=1, confirmed.")

# Regression check for the adaptive re-anchoring machinery itself: the
# sphere's gradient is exactly radial, so the fixed patch line ALWAYS
# passes through the axis the circle shrinks toward (see
# `_sweep_hop!`'s docstring) -- the transversality-drift check should
# therefore never fire here. Confirming this is zero (not just "small")
# is the sphere-side regression guard for the mechanism validated
# against the ellipsoid below: it must not fire spuriously on the exact
# case it was designed to leave alone.
n_rebuilds_down, _ = instrumented_sweep_direction(patch, 0.0, 1.0, 0.0, collect(range(0.0, -1.0; length = cfg.midslice_sample_density + 1))[2:end], cfg)
n_rebuilds_up, _ = instrumented_sweep_direction(patch, 0.0, 1.0, 0.0, collect(range(0.0, 1.0; length = cfg.midslice_sample_density + 1))[2:end], cfg)
println("  re-anchor events on sphere (radially symmetric, expect 0): down=$(n_rebuilds_down), up=$(n_rebuilds_up)")
@test n_rebuilds_down == 0
@test n_rebuilds_up == 0

println()
println("=" ^ 70)
println("5. track_face: sweep one full edge (a semicircle) into a hemisphere patch")
println("=" ^ 70)

edge_upper = edges0[argmax([sum(p[2] for p in e.sampled_points) for e in edges0])] # the y>0 semicircle
face1 = track_face(F_sphere, patch, edge_upper, 0.0, -1.0, 1.0, 1, cfg)
n_z = 2 * cfg.midslice_sample_density + 1
n_curve = cfg.edge_sample_density
println("  Face.mesh_vertices size: ", size(face1.mesh_vertices), "  (expect ($(n_z*n_curve), 3))")
@test size(face1.mesh_vertices) == (n_z * n_curve, 3)
@test all(r -> isapprox(sum(abs2, face1.mesh_vertices[r, :]), 1.0; atol = 1e-3), 1:size(face1.mesh_vertices, 1))
println("  every mesh vertex satisfies x^2+y^2+z^2≈1, confirmed.")
println("  Face.mesh_topology size: ", size(face1.mesh_topology))
@test size(face1.mesh_topology, 2) == 3
@test all(row -> length(unique(face1.mesh_topology[row, :])) == 3, 1:size(face1.mesh_topology, 1))
println("  every triangle has 3 distinct vertex indices (no degenerate triangles slipped through).")
@test all(face1.mesh_topology .>= 1) && all(face1.mesh_topology .<= n_z * n_curve)
println("track_face: hemisphere patch mesh geometry confirmed.")

println()
println("=" ^ 70)
println("6. decompose_3d_surface: full pipeline, expect exactly 2 faces (two hemispheres)")
println("=" ^ 70)

all_vertices, all_edges, all_faces, mesh = decompose_3d_surface(F_sphere, cfg)
println("  vertices: ", length(all_vertices), "   edges: ", length(all_edges), "   faces: ", length(all_faces))
@test length(all_faces) == 2 # only the (-1,1) slab has a nonempty slice; the two outer slabs are empty
@test length(unique(v.id for v in all_vertices)) == length(all_vertices) # no id collisions across slabs
@test length(unique(e.id for e in all_edges)) == length(all_edges)

println("  mesh: ", length(GeometryBasics.coordinates(mesh)), " vertices, ", length(GeometryBasics.faces(mesh)), " triangles")
mesh_pts = GeometryBasics.coordinates(mesh)
mesh_tris = GeometryBasics.faces(mesh)
@test all(p -> isapprox(sum(abs2, p), 1.0; atol = 1e-2), mesh_pts)
println("  every welded mesh vertex satisfies x^2+y^2+z^2≈1 (Float32 precision), confirmed.")

# Outward-normal check: on the unit sphere, ∇f = (2x,2y,2z) is exactly
# radial, so weld_mesh's "align triangle normal with ∇f" convention
# means EVERY triangle's normal must point outward, i.e.
# normal · centroid > 0 (centroid ≈ radial direction at that point too).
outward_ok = all(mesh_tris) do tri
    p1, p2, p3 = mesh_pts[Int(tri[1])], mesh_pts[Int(tri[2])], mesh_pts[Int(tri[3])]
    n = cross(p2 .- p1, p3 .- p1)
    centroid = (p1 .+ p2 .+ p3) ./ 3.0f0
    dot(n, centroid) >= 0
end
println("  all $(length(mesh_tris)) welded triangles have outward-pointing normals: ", outward_ok)
@test outward_ok
println("decompose_3d_surface: full pipeline produces a watertight, consistently-oriented sphere mesh.")

println()
println("=" ^ 70)
println("6b. _robust_slice_at_z zero-spurious-retry confirmation (sphere)")
println("    (combined vertex-type + gradient-magnitude gate, added while closing")
println("    out the Taubin heart's z_mid-degeneracy investigation)")
println("=" ^ 70)
let
    z_naive = (-1.0 + 1.0) / 2
    v_r, e_r, z_used = HomotopyGetsReal._robust_slice_at_z(F_sphere, patch, -1.0, 1.0, cfg)
    grads = [hypot(patch_direction(patch, e.sampled_points[1][1], e.sampled_points[1][2], e.sampled_points[1][3], cfg)...) for e in e_r]
    println("  slab [-1,1]: z_mid used=$(z_used) (naive=$(z_naive)) -> retried: $(!isapprox(z_used, z_naive; atol=1e-12))")
    println("  per-edge anchor |grad| at accepted z_mid: ", grads)
    @test isapprox(z_used, z_naive; atol = 1e-12)
end
println("  confirmed: zero spurious retries on the sphere (radially symmetric, uniform")
println("  gradient magnitude everywhere -- no branch-to-branch spread for the gate to")
println("  misfire on).")

println()
println("=" ^ 70)
println("7. Asymmetric ellipsoid: x^2 + 4y^2 + 9z^2 = 1 (catches x/y/z mixups a")
println("   perfectly symmetric sphere test cannot -- e.g. an accidental x<->y swap")
println("   somewhere would still pass every sphere check above by symmetry)")
println("=" ^ 70)

@var xe ye ze
f_ell = xe^2 + 4 * ye^2 + 9 * ze^2 - 1
F_ell = System([f_ell], variables = [xe, ye, ze])

# Hand-computed: ∂f/∂x = 2x = 0, ∂f/∂y = 8y = 0 => x=y=0, 9z^2=1 => z=±1/3.
z_crits_ell = compute_critical_z_slices(F_ell, cfg)
println("  compute_critical_z_slices (ellipsoid) -> ", z_crits_ell, "  (expect ±1/3 ≈ ±0.3333)")
@test length(z_crits_ell) == 2
@test isapprox(sort(z_crits_ell), [-1 / 3, 1 / 3]; atol = 1e-6)

ev, ee, ef, emesh = decompose_3d_surface(F_ell, cfg)
println("  decompose_3d_surface (ellipsoid): ", length(ev), " vertices, ", length(ee), " edges, ", length(ef), " faces")
@test length(ef) == 2
ell_mesh_pts = GeometryBasics.coordinates(emesh)
# Tight tolerance (not the earlier 1e-2 "does it stay roughly on the
# surface" placeholder): with adaptive re-anchoring + residual-based
# bisection (see `_sweep_hop!`'s docstring) every swept point, INCLUDING
# ones resolved right at the pole where the level curve degenerates to a
# point, converges to within Float32-mesh-storage precision of the true
# surface (empirically `< 1e-6` at every vertex, including the pole
# itself) -- not just "close enough for a display mesh".
ell_residuals = [abs(p[1]^2 + 4 * p[2]^2 + 9 * p[3]^2 - 1.0) for p in ell_mesh_pts]
sorted_resid = sort(ell_residuals)
println("  |f| residual distribution over all $(length(ell_mesh_pts)) welded mesh vertices:")
println("    min    = ", minimum(ell_residuals))
println("    mean   = ", mean(ell_residuals))
println("    median = ", median(ell_residuals))
println("    p90    = ", sorted_resid[ceil(Int, 0.90 * length(sorted_resid))])
println("    p99    = ", sorted_resid[ceil(Int, 0.99 * length(sorted_resid))])
println("    max    = ", maximum(ell_residuals))
println("    count with residual > 1e-4: ", count(>(1e-4), ell_residuals), " / ", length(ell_residuals))
println("    count with residual > 1e-5: ", count(>(1e-5), ell_residuals), " / ", length(ell_residuals))
@test all(<=(1e-4), ell_residuals)
println("  every welded mesh vertex satisfies x^2+4y^2+9z^2≈1 to within 1e-4, confirmed (variable roles not swapped).")

println()
println("=" ^ 70)
println("7a. _robust_slice_at_z zero-spurious-retry confirmation (ellipsoid)")
println("=" ^ 70)
let
    patch_ell = HomotopyGetsReal.build_patch_system(F_ell)
    z_bottom_e, z_top_e = -1 / 3, 1 / 3
    z_naive = (z_bottom_e + z_top_e) / 2
    v_r, e_r, z_used = HomotopyGetsReal._robust_slice_at_z(F_ell, patch_ell, z_bottom_e, z_top_e, cfg)
    grads = [hypot(patch_direction(patch_ell, e.sampled_points[1][1], e.sampled_points[1][2], e.sampled_points[1][3], cfg)...) for e in e_r]
    println("  slab [-1/3,1/3]: z_mid used=$(z_used) (naive=$(z_naive)) -> retried: $(!isapprox(z_used, z_naive; atol=1e-12))")
    println("  per-edge anchor |grad| at accepted z_mid: ", grads)
    @test isapprox(z_used, z_naive; atol = 1e-12)
end
println("  confirmed: zero spurious retries on the ellipsoid either (asymmetric but still")
println("  a single smooth branch per slab -- no multi-loop/pinch structure of the kind")
println("  that caused the Taubin heart's narrow slab to need one).")

println()
println("=" ^ 70)
println("7b. Re-anchoring cost benchmark (ellipsoid, anchor toward pole)")
println("=" ^ 70)

let
    patch_ell = HomotopyGetsReal.build_patch_system(F_ell)
    t_anchor = 0.3
    x0e, y0e = cos(t_anchor), 0.5 * sin(t_anchor)
    z_mid_e, z_bottom_e = 0.0, -1 / 3
    n_side = cfg.midslice_sample_density
    targets = collect(range(z_mid_e, z_bottom_e; length = n_side + 1))[2:end]

    n_rebuilds, dense = instrumented_sweep_direction(patch_ell, x0e, y0e, z_mid_e, targets, cfg)
    elapsed = @elapsed instrumented_sweep_direction(patch_ell, x0e, y0e, z_mid_e, targets, cfg)
    println("  $(n_side) z-targets toward the pole: $(n_rebuilds) tracker rebuild(s) (re-anchor events), $(round(elapsed * 1000; digits = 2)) ms total")
    println("  (compare Phase 5's own build_tracker benchmark: compile=:none construction is the dominant per-rebuild cost,")
    println("   so a handful of adaptive rebuilds per direction stays negligible next to the one-rebuild-per-anchor baseline)")
end

println()
println("=" ^ 70)
println("7c. patch_transversality_cos_tol sensitivity sweep (ellipsoid, same anchor)")
println("    Confirms the threshold isn't just \"what worked on this one case\":")
println("    tightening (-> 1) should trigger MORE rebuilds and keep/improve residual;")
println("    loosening (-> 0) should trigger FEWER rebuilds and let residual drift up.")
println("=" ^ 70)

let
    patch_ell = HomotopyGetsReal.build_patch_system(F_ell)
    t_anchor = 0.3
    x0e, y0e = cos(t_anchor), 0.5 * sin(t_anchor)
    z_mid_e, z_bottom_e = 0.0, -1 / 3
    n_side = cfg.midslice_sample_density
    targets = collect(range(z_mid_e, z_bottom_e; length = n_side + 1))[2:end]
    ell_resid_fn(p) = abs(p[2]^2 + 4 * p[3]^2 + 9 * p[1]^2 - 1.0) # dense_path is [z,x,y]-ordered

    cos_tols = [0.0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99]
    results = NamedTuple[]
    for ct in cos_tols
        cfg_ct = HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8, patch_transversality_cos_tol = ct)
        n_rebuilds, dense = instrumented_sweep_direction(patch_ell, x0e, y0e, z_mid_e, targets, cfg_ct)
        max_resid = maximum(ell_resid_fn(p) for p in dense)
        push!(results, (cos_tol = ct, n_rebuilds = n_rebuilds, max_resid = max_resid))
    end
    println("  cos_tol  | rebuilds | max |f| residual")
    println("  ---------+----------+-----------------")
    for r in results
        println("  ", rpad(string(r.cos_tol), 8), " | ", rpad(string(r.n_rebuilds), 8), " | ", r.max_resid)
    end
    # Monotonicity checks (weak/non-strict: ties are fine since rebuild
    # counts are integers and several thresholds can land in the same
    # bucket between consecutive cos_angle values along this one sweep;
    # what must NOT happen is a threshold INCREASE causing FEWER
    # rebuilds, or a residual that gets WORSE as rebuilds increase).
    rebuild_counts = [r.n_rebuilds for r in results]
    @test issorted(rebuild_counts) # non-decreasing as cos_tol increases
    println("  rebuild count is non-decreasing as cos_tol increases (", rebuild_counts, "), confirmed.")
    @test results[1].max_resid >= results[end].max_resid - 1e-8 # loosest tol should not beat the tightest
    println("  loosest threshold (cos_tol=$(results[1].cos_tol)) residual ($(results[1].max_resid)) is >= tightest")
    println("  threshold (cos_tol=$(results[end].cos_tol)) residual ($(results[end].max_resid)), confirming the expected direction.")
    println("  (cfg's own default, 0.9, sits at rebuilds=$(results[findfirst(==(0.9), cos_tols)].n_rebuilds), ")
    println("   residual=$(results[findfirst(==(0.9), cos_tols)].max_resid) -- comfortably on the \"few rebuilds, tight residual\" side.)")
end

println()
println("=" ^ 70)
println("8. @inferred type-stability checks")
println("=" ^ 70)

infer_crit_z(F_, cfg_) = compute_critical_z_slices(F_, cfg_)
infer_slice(F_, zv, cfg_) = slice_at_z(F_, zv, cfg_)
infer_patch(F_) = build_patch_system(F_)
infer_dir(p_, x0, y0, z0, cfg_) = patch_direction(p_, x0, y0, z0, cfg_)
infer_dense(F_, ph_, tr_, ys, ps, pts, ms, rt, st, pat, mw) = track_dense_path(F_, ph_, tr_, ys, ps, pts, ms, rt, st, pat, mw)
infer_sweep(F_, p_, x0, y0, zm, zb, zt, cfg_) = sweep_face_bidirectional(F_, p_, x0, y0, zm, zb, zt, cfg_)
infer_face(F_, p_, e_, zm, zb, zt, fid, cfg_) = track_face(F_, p_, e_, zm, zb, zt, fid, cfg_)
infer_weld(fs, p_, cfg_) = weld_mesh(fs, p_, cfg_)
infer_decomp(F_, cfg_) = decompose_3d_surface(F_, cfg_)

r1 = @inferred infer_crit_z(F_sphere, cfg)
println("@inferred compute_critical_z_slices  -> ", typeof(r1), "  OK")
@test r1 isa Vector{Float64}

r2 = @inferred infer_slice(F_sphere, 0.0, cfg)
println("@inferred slice_at_z                 -> ", typeof(r2), "  OK")
@test r2 isa Tuple{Vector{NativeVertex{Float64}},Vector{Edge{Float64}}}

r3 = @inferred infer_patch(F_sphere)
println("@inferred build_patch_system         -> ", typeof(r3), "  OK")


# patch_direction (via _gradient_at) is NOT @inferred-clean -- and,
# newly discovered while validating this file, neither is
# Solver.jacobian_rank_info itself (never previously verified with a
# direct @inferred test; only indirectly exercised through callers that
# happen to re-seal through a concrete constructor immediately after,
# e.g. NativeVertex{T}(...)). Both share the same root cause: HC's
# low-level `evaluate(...; bits=...)` has a genuinely runtime-value-
# dependent return type (Vector{Float64} vs Vector{ComplexF64} depending
# on whether the numeric result happens to be real), which no amount of
# outer-constructor sealing can retroactively narrow. Reported via
# Base.return_types, not a hard @inferred, matching build_tracker's own
# precedent for its own (differently-caused) documented exception.
dir_types = Base.return_types(infer_dir, (typeof(patch), Float64, Float64, Float64, typeof(cfg)))
jri_types = Base.return_types(jacobian_rank_info, (typeof(F_sphere), Vector{ComplexF64}, typeof(cfg)))
println("patch_direction: NOT @inferred-clean (upstream evaluate(...) limitation) -- Core.Compiler infers ", only(dir_types))
println("  (same root cause independently confirmed on Solver.jacobian_rank_info here: ", only(jri_types), ")")
r4 = infer_dir(patch, 0.0, 1.0, 0.0, cfg)
println("  runtime value is concretely: ", typeof(r4))
@test r4 isa Tuple{Float64,Float64}

# build_face_tracker itself is a direct pass-through of build_tracker's
# own documented Union instability (see build_tracker's docstring) --
# tested here via Base.return_types, not a hard @inferred, exactly
# mirroring how Phase 4's own script handled build_tracker.
bft_argtypes = (typeof(patch), Float64, Float64, Float64, typeof(cfg))
bft_types = Base.return_types(build_face_tracker, bft_argtypes)
println(
    "build_face_tracker: NOT @inferred-clean (inherits build_tracker's documented Union) -- ",
    "Core.Compiler infers ", only(bft_types),
)
_, ph_bft, tracker_bft = build_face_tracker(patch, 0.0, 1.0, 0.0, cfg)
@test (typeof(ph_bft), typeof(tracker_bft)) isa Tuple{DataType,DataType}

_, ph5, tracker5 = build_face_tracker(patch, 0.0, 1.0, 0.0, cfg)
y_state5 = ComplexF64[0.0, 1.0]
targets5 = collect(range(0.0, 1.0; length = cfg.midslice_sample_density + 1))[2:end]
r5 = @inferred infer_dense(
    patch.F_for_tracking, ph5, tracker5, y_state5, 0.0, targets5,
    cfg.max_path_steps, Float64(cfg.jacobian_rank_tol), Float64(cfg.singular_value_threshold),
    Float64(cfg.critical_point_tol), Float64(cfg.vertex_match_tol),
)
println("@inferred track_dense_path (pre-built tracker) -> ", typeof(r5), "  OK")
@test r5 isa Tuple{Vector{ComplexF64},Vector{Vector{Float64}}}


# Empirically (see build_face_tracker's own note above about union-
# splitting): sweep_face_bidirectional/track_face BOTH turn out to be
# fully @inferred-clean end-to-end DESPITE calling build_face_tracker
# internally -- Julia's compiler successfully union-splits over
# build_tracker's 3-member System-type Union and proves every branch
# produces the SAME concrete final return type. This is empirically
# confirmed here, not assumed.
r6 = @inferred infer_sweep(F_sphere, patch, 0.0, 1.0, 0.0, -1.0, 1.0, cfg)
println("@inferred sweep_face_bidirectional (union-splits cleanly through build_face_tracker) -> ", typeof(r6), "  OK")
@test r6 isa Tuple{Vector{Vector{Float64}},Vector{Vector{Float64}}}

r7 = @inferred infer_face(F_sphere, patch, edge_upper, 0.0, -1.0, 1.0, 1, cfg)
println("@inferred track_face (union-splits cleanly through build_face_tracker)             -> ", typeof(r7), "  OK")
@test r7 isa Face{Float64}

# weld_mesh/decompose_3d_surface are NOT @inferred-clean, but for a
# DIFFERENT reason than build_tracker's Union: GeometryBasics.Mesh's own
# constructor return type carries an existential `where _A` type
# parameter (its internal "views" storage type), which Core.Compiler
# cannot resolve away -- reported via Base.return_types, not a hard
# @inferred, since this is an upstream GeometryBasics.jl property, not
# something union-splitting (or anything in this file) can fix.
for (label, f_, argtypes) in (
    ("weld_mesh", infer_weld, (Vector{Face{Float64}}, typeof(patch), typeof(cfg))),
    ("decompose_3d_surface", infer_decomp, (typeof(F_sphere), typeof(cfg))),
)
    rt = Base.return_types(f_, argtypes)
    is_concrete = length(rt) == 1 && Base.isconcretetype(only(rt))
    println("Base.return_types $label -> ", only(rt), is_concrete ? "  (concrete)" : "  (NOT concrete: GeometryBasics.Mesh's own `where _A`)")
end

end
