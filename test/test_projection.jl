@testset "Projection (Phase 8)" begin
    using Random
    using LinearAlgebra

println("=" ^ 70)
println("1. random_orthogonal_matrix: SO(3) membership, Haar construction, reproducibility")
println("=" ^ 70)

for seed in 1:5
    Q = random_orthogonal_matrix(Float64, 3; rng = Xoshiro(seed))
    @test size(Q) == (3, 3)
    @test norm(Q' * Q - I) < 1e-12
    @test isapprox(det(Q), 1.0; atol = 1e-12)
end
println("5 seeds: orthonormal to <1e-12, det = +1 (SO(3), never a reflection).")

Q_rep_a = random_orthogonal_matrix(Float64, 3; rng = Xoshiro(42))
Q_rep_b = random_orthogonal_matrix(Float64, 3; rng = Xoshiro(42))
@test Q_rep_a == Q_rep_b
Q_big = random_orthogonal_matrix(BigFloat, 3; rng = Xoshiro(1))
@test Q_big isa Matrix{BigFloat}
@test Float64.(Q_big) == random_orthogonal_matrix(Float64, 3; rng = Xoshiro(1))
println("Seeded reproducibility and T-genericity (BigFloat entries exactly represent the Float64 draw) confirmed.")

println()
println("=" ^ 70)
println("2. projection-argument validation (loud ArgumentErrors, house style)")
println("=" ^ 70)

@var x y z
f_sph = x^2 + y^2 + z^2 - 1
F_sph = System([f_sph], variables = [x, y, z])
cfg = HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8)

@test_throws ArgumentError decompose_3d_surface(F_sph, cfg; projection = :bogus)
@test_throws ArgumentError decompose_3d_surface(F_sph, cfg; projection = ones(2, 2))
@test_throws ArgumentError decompose_3d_surface(F_sph, cfg; projection = 2.0 .* Matrix(1.0I, 3, 3))
Q_reflect = Matrix(1.0I, 3, 3)
Q_reflect[1, 1] = -1.0
@test_throws ArgumentError decompose_3d_surface(F_sph, cfg; projection = Q_reflect)
println("Rejected: unknown symbol, wrong size, non-orthonormal, det = -1 reflection.")

@test HomotopyConfig{Float64}().min_slab_width == 1e-3
@test HomotopyConfig{BigFloat}().min_slab_width isa BigFloat
println("min_slab_width config field present and T-typed.")

println()
println("=" ^ 70)
println("3. _verify_projection_ok: catches the z - x^2 crash class, does NOT over-scope")
println("=" ^ 70)

f_par = z - x^2
F_par = System([f_par], variables = [x, y, z])

# Degenerate identity projection: ∂f'/∂y ≡ 0 must become a loud ArgumentError
# naming the problem, not HC's OverflowError from deep inside solve.
err = try
    decompose_3d_surface(F_par, cfg; projection = Matrix(1.0I, 3, 3))
    nothing
catch e
    e
end
@test err isa ArgumentError
@test occursin("degenerate", sprint(showerror, err))
println("Identity projection on z - x^2: ArgumentError naming the vanishing partial. ✓")

# Deliberate scoping check: a chart-z-independent surface (cylinder, fz ≡ 0)
# is a WORKING input class (empty critical system) and must NOT be flagged --
# only the two AUGMENTING partials are checked.
F_cyl = System([x^2 + y^2 - 1], variables = [x, y, z])
@test HomotopyGetsReal._verify_projection_ok(F_cyl, cfg) === nothing
println("Cylinder (fz ≡ 0) passes the check: only augmenting partials are gated. ✓")

println()
println("=" ^ 70)
println("4. z - x^2 regression: crashes bare, works under projection = :random")
println("=" ^ 70)

# The Phase 8 motivating bug: bare fixed-axis decomposition dies inside HC's
# start-system construction (OverflowError as of HC 2.17).
@test_throws Exception decompose_3d_surface(F_par, cfg)
println("Bare decompose_3d_surface(z - x^2) still throws (the pre-Phase-8 crash), confirmed.")

vp, ep, fp, mp = decompose_3d_surface(F_par, cfg; projection = :random, rng = Xoshiro(42))
pts_p = GeometryBasics.coordinates(mp)
@test !isempty(fp)
@test !isempty(pts_p)
res_p = maximum(abs(Float64(p[3]) - Float64(p[1])^2) for p in pts_p)
println("projection = :random (seed 42): $(length(fp)) faces, $(length(pts_p)) mesh points, max world |z - x^2| = $(res_p)")
@test res_p < 1e-3
println("Every mapped-back mesh point satisfies the ORIGINAL equation: crash class closed. ✓")

println()
println("=" ^ 70)
println("5. sphere under projection = :random: well-formed world-frame mesh")
println("=" ^ 70)

vs, es, fs, ms = decompose_3d_surface(F_sph, cfg; projection = :random, rng = Xoshiro(7))
pts_s = GeometryBasics.coordinates(ms)
tris_s = GeometryBasics.faces(ms)
@test !isempty(fs) && !isempty(pts_s) && !isempty(tris_s)
# Rotation-invariant world check: the mapped-back mesh must still be the unit sphere.
@test all(p -> isapprox(sum(abs2, p), 1.0; atol = 1e-2), pts_s)
@test all(t -> length(unique((Int(t[1]), Int(t[2]), Int(t[3])))) == 3, tris_s)
outward_ok = all(tris_s) do t
    p1, p2, p3 = pts_s[Int(t[1])], pts_s[Int(t[2])], pts_s[Int(t[3])]
    n = cross(p2 .- p1, p3 .- p1)
    dot(n, (p1 .+ p2 .+ p3) ./ 3.0f0) >= 0
end
@test outward_ok
@test length(unique(v.id for v in vs)) == length(vs)
println("Rotated sphere: $(length(fs)) faces, $(length(pts_s)) mesh points, on-sphere to 1e-2, no degenerate")
println("triangles, all normals outward (det = +1 preserves the winding convention through map-back).")
println("(Face count deliberately not hard-asserted: chart topology may differ from the fixed-axis run.)")

if get(ENV, "HOMOTOPYGETSREAL_RUN_SLOW_TESTS", "0") == "1"
    println()
    println("=" ^ 70)
    println("6. SLOW: rotated Taubin heart, 5 fixed seeds -- the Phase 8 success criterion")
    println("   (the z=0 repeated-factor degeneracy signature must not fire structurally;")
    println("    report actual per-slab retry/throw counts, not just \"it worked\")")
    println("=" ^ 70)

    @var xh yh zh
    f_heart = (xh^2 + (1.2 * yh)^2 + zh^2 - 1)^3 - xh^2 * zh^3 - 0.1 * (1.2 * yh)^2 * zh^3
    F_heart = System([f_heart], variables = [xh, yh, zh])
    cfg_heart = HomotopyConfig{Float64}(
        bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3),
        edge_sample_density = 8, midslice_sample_density = 8,
    )

    for seed in 1:5
        Q = random_orthogonal_matrix(Float64, 3; rng = Xoshiro(seed))
        F_rot = HomotopyGetsReal._rotate_system(F_heart, Q)
        cfg_rot = HomotopyGetsReal._chart_config(cfg_heart, Q)
        patch_rot = HomotopyGetsReal.build_patch_system(F_rot)
        z_bounds = HomotopyGetsReal._slab_bounds(F_rot, cfg_rot)

        n_retried = 0
        n_threw = 0
        slab_report = String[]
        for i in 1:(length(z_bounds) - 1)
            zb, zt = z_bounds[i], z_bounds[i + 1]
            try
                _, e_r, z_used = HomotopyGetsReal._robust_slice_at_z(F_rot, patch_rot, zb, zt, cfg_rot)
                retried = !isapprox(z_used, (zb + zt) / 2; atol = 1e-12)
                retried && (n_retried += 1)
                push!(slab_report, "[$(round(zb; digits = 3)),$(round(zt; digits = 3))]:$(length(e_r))e$(retried ? "*RETRY" : "")")
            catch
                n_threw += 1
                push!(slab_report, "[$(round(zb; digits = 3)),$(round(zt; digits = 3))]:THREW")
            end
        end
        println("  seed $seed: $(length(z_bounds) - 1) slabs | retried = $n_retried | threw = $n_threw")
        println("    ", join(slab_report, "  "))
        @test n_threw == 0
    end
    println("  Zero throws across all seeds: the transversal-singular-curve gate refinements hold")
    println("  beyond the seed they were diagnosed on.")

    println()
    println("  --- seed 1: full decompose + world-residual localization ---")
    Q1 = random_orthogonal_matrix(Float64, 3; rng = Xoshiro(1))
    t_rot = @elapsed (vh, eh, fh, mh) = decompose_3d_surface(F_heart, cfg_heart; projection = Q1)
    pts_h = GeometryBasics.coordinates(mh)
    f_ev(px, py, pz) = (px^2 + (1.2 * py)^2 + pz^2 - 1)^3 - px^2 * pz^3 - 0.1 * (1.2 * py)^2 * pz^3
    res_h = [abs(f_ev(Float64(p[1]), Float64(p[2]), Float64(p[3]))) for p in pts_h]
    sorted_h = sort(res_h)
    n_bad = count(>(1e-4), res_h)
    println("  [$(round(t_rot; digits = 1)) s]  $(length(vh)) verts, $(length(eh)) edges, $(length(fh)) faces, $(length(pts_h)) mesh pts")
    println("  world |f|: median = $(sorted_h[cld(end, 2)])  p99 = $(sorted_h[ceil(Int, 0.99 * length(sorted_h))])  max = $(maximum(res_h))")
    println("  points > 1e-4: $n_bad / $(length(res_h))")
    @test sorted_h[cld(end, 2)] < 1e-6           # bulk of the mesh is excellent
    @test n_bad / length(res_h) < 0.05           # defects are rare ...
    bad_pts = [p for (p, r) in zip(pts_h, res_h) if r > 1e-4]
    if !isempty(bad_pts)
        band = maximum(abs(Float64(p[3])) for p in bad_pts)
        println("  localization: all >1e-4 points have |z_world| <= $(round(band; digits = 3)) (surface spans |z| <= 1.24)")
        @test band <= 0.2                        # ... and confined to the singular-curve band
    end
    println("  Documented known limitation confirmed quantitatively: residual defects exist only in")
    println("  the band around the transversal singular ellipse (no singular-curve decomposition yet).")
end

end
