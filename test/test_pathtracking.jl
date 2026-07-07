@testset "PathTracking (Phase 4)" begin

println("=" ^ 70)
println("Setup: nodal cubic f(x,y) = y^2 - x^3 - x^2 (same as Phase 2/3)")
println("=" ^ 70)

@var x y
f = y^2 - x^3 - x^2
F_curve = System([f], variables = [x, y])
cfg64 = HomotopyConfig{Float64}()

# Derive the known vertex x-coordinates from Phase 2's own functions
# (NOT connect_the_dots!/decompose_1d_curve) so nothing here is a magic
# literal duplicated from the Phase 3 script's printed output.
F_aug = System([f, differentiate(f, y)], [x, y])
crit_pts = compute_critical_points(F_aug, cfg64)
bnd_pts = intersect_bounding_object(F_curve, cfg64)
x_boundary = real(bnd_pts[1].coordinates[1])
println("Derived boundary x-coordinate: ", x_boundary)

rank_tol64 = Float64(cfg64.jacobian_rank_tol)
sv_thresh64 = Float64(cfg64.singular_value_threshold)
poor_acc_tol64 = Float64(cfg64.critical_point_tol)
min_width64 = Float64(cfg64.vertex_match_tol)

H_sys = System([f], variables = [y], parameters = [x])

println()
println("=" ^ 70)
println("1. track_path / track_bidirectional directly, bypassing connect_the_dots!")
println("=" ^ 70)

println("\n-- Interval (-1, 0): Critical@(-1,0) <-> Singular@(0,0) --")
x_left1, x_right1 = -1.0, 0.0
x_mid1 = (x_left1 + x_right1) / 2
y_mids1 = compute_midslice(F_curve, x_left1, x_right1, cfg64)
println("compute_midslice($x_left1, $x_right1) -> ", y_mids1)
@test length(y_mids1) == 2

ph1, tracker1 = build_tracker(H_sys, x_mid1, cfg64)
y_mid1_vec = ComplexF64[y_mids1[1]]
full_path1, y_land_left1, y_land_right1 = track_bidirectional(
    F_curve, ph1, tracker1, y_mid1_vec, x_mid1, x_left1, x_right1,
    cfg64.max_path_steps, rank_tol64, sv_thresh64, poor_acc_tol64, min_width64,
)
println("  landing at x_left=$x_left1  -> y=", y_land_left1, "  (expect ~0, the Critical vertex)")
println("  landing at x_right=$x_right1 -> y=", y_land_right1, "  (expect ~0, the Singular node)")
println("  full_path length: ", length(full_path1))
@test isapprox(real(y_land_left1[1]), 0.0; atol = 1e-3)
@test isapprox(real(y_land_right1[1]), 0.0; atol = 1e-3)
@test length(full_path1) >= 2

println("\n-- Interval (0, $x_boundary): Singular@(0,0) <-> Boundary@($x_boundary, ±4) --")
x_left2, x_right2 = 0.0, x_boundary
x_mid2 = (x_left2 + x_right2) / 2
y_mids2 = compute_midslice(F_curve, x_left2, x_right2, cfg64)
println("compute_midslice($x_left2, $x_right2) -> ", y_mids2)
@test length(y_mids2) == 2

for y_mid in y_mids2
    ph2, tracker2 = build_tracker(H_sys, x_mid2, cfg64)
    full_path2, y_land_left2, y_land_right2 = track_bidirectional(
        F_curve, ph2, tracker2, ComplexF64[y_mid], x_mid2, x_left2, x_right2,
        cfg64.max_path_steps, rank_tol64, sv_thresh64, poor_acc_tol64, min_width64,
    )
    println("  y_mid=$y_mid -> landing left(x=0): ", y_land_left2, "   landing right(x=$x_right2): ", y_land_right2)
    @test isapprox(real(y_land_left2[1]), 0.0; atol = 1e-3)
    @test isapprox(abs(real(y_land_right2[1])), 4.0; atol = 1e-3)
end
println("track_path/track_bidirectional landing values confirmed against Phase 2's known vertices.")

println()
println("=" ^ 70)
println("2. is_near_singular checks (both rank_tol AND sv_thresh, per _classify_vertex_type)")
println("=" ^ 70)

expected_rank = length(F_curve.expressions) # = 1, a plane curve
at_node = is_near_singular(F_curve, ComplexF64[0.0, 0.0], expected_rank, rank_tol64, sv_thresh64)
at_critical = is_near_singular(F_curve, ComplexF64[-1.0, 0.0], expected_rank, rank_tol64, sv_thresh64)
at_smooth = is_near_singular(F_curve, ComplexF64[x_mid1, y_mids1[1]], expected_rank, rank_tol64, sv_thresh64)

println("  is_near_singular at (0,0)              [genuine node]        -> ", at_node, "   (expect true)")
println("  is_near_singular at (-1,0)              [Critical, full rank] -> ", at_critical, "   (expect false)")
println("  is_near_singular at ($x_mid1, $(y_mids1[1])) [smooth interior]   -> ", at_smooth, "   (expect false)")

@test at_node == true
@test at_critical == false
@test at_smooth == false
println("is_near_singular correctly distinguishes the genuine node from Critical/smooth points.")

println()
println("=" ^ 70)
println("3. build_tracker threads cfg.path_tracker_precision into TrackerOptions.min_step_size")
println("=" ^ 70)

ph_default, tracker_default = build_tracker(H_sys, x_mid1, cfg64)
println("  default cfg.path_tracker_precision = ", cfg64.path_tracker_precision)
println("  tracker.options.min_step_size      = ", tracker_default.options.min_step_size)
@test tracker_default.options.min_step_size == Float64(cfg64.path_tracker_precision)

cfg_tight_prec = HomotopyConfig{Float64}(path_tracker_precision = 1e-20)
ph_tight, tracker_tight = build_tracker(H_sys, x_mid1, cfg_tight_prec)
println("  tightened cfg.path_tracker_precision = ", cfg_tight_prec.path_tracker_precision)
println("  tracker.options.min_step_size        = ", tracker_tight.options.min_step_size)
@test tracker_tight.options.min_step_size == 1e-20
println("Wiring confirmed: min_step_size exactly reflects cfg.path_tracker_precision.")
println("Note (per build_tracker's own docstring caveat): this confirms the WIRING, not a")
println("guaranteed change in *tracking outcome* on this simple curve -- HomotopyContinuation.jl")
println("exposes no literal requested-accuracy knob, min_step_size is the closest available lever.")

println()
println("=" ^ 70)
println("4. @inferred type-stability checks")
println("=" ^ 70)

infer_near_sing(F_, pt, er, rt, st) = is_near_singular(F_, pt, er, rt, st)
infer_build(hs, xs, cfg_) = build_tracker(hs, xs, cfg_)
infer_track(F_, ph_, tr_, ys, xs, xt, ms, rt, st, pat, mw) = track_path(F_, ph_, tr_, ys, xs, xt, ms, rt, st, pat, mw)
infer_bidir(F_, ph_, tr_, ym, xm, xl, xr, ms, rt, st, pat, mw) =
    track_bidirectional(F_, ph_, tr_, ym, xm, xl, xr, ms, rt, st, pat, mw)

r1 = @inferred infer_near_sing(F_curve, ComplexF64[0.0, 0.0], expected_rank, rank_tol64, sv_thresh64)
println("@inferred is_near_singular      -> ", typeof(r1), "  OK")

r2 = infer_build(H_sys, x_mid1, cfg64)
r2_ok = Base.return_types(infer_build, (System,Float64,typeof(cfg64)))
println(
    "build_tracker: NOT @inferred-clean (documented, upstream) -- ",
    "Core.Compiler infers ", only(r2_ok), " (a Union over CompiledSystem/InterpretedSystem/MixedSystem)",
)
println("  runtime value is concretely: ", typeof(r2))
@test r2 isa Tuple{ParameterHomotopy,Tracker}

ph3, tracker3 = build_tracker(H_sys, x_mid1, cfg64)
r3 = @inferred infer_track(
    F_curve, ph3, tracker3, y_mid1_vec, x_mid1, x_left1,
    cfg64.max_path_steps, rank_tol64, sv_thresh64, poor_acc_tol64, min_width64,
)
println("@inferred track_path            -> ", typeof(r3), "  OK")

ph4, tracker4 = build_tracker(H_sys, x_mid1, cfg64)
r4 = @inferred infer_bidir(
    F_curve, ph4, tracker4, y_mid1_vec, x_mid1, x_left1, x_right1,
    cfg64.max_path_steps, rank_tol64, sv_thresh64, poor_acc_tol64, min_width64,
)
println("@inferred track_bidirectional   -> ", typeof(r4), "  OK")

@test r1 isa Bool
@test r3 isa Tuple{Vector{ComplexF64},Vector{Vector{Float64}}}
@test r4 isa Tuple{Vector{Vector{Float64}},Vector{ComplexF64},Vector{ComplexF64}}

end
