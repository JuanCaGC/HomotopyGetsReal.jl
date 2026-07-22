@testset "Taubin heart (Phase 5 slow)" begin

"""
    _instrumented_sweep_hop!(hop_path, patch, cfg, state, y0, p0, p1, budget, tols, inner_bisections, outer_bisections)

Read-only "shadow" reimplementation of `FaceTracking._sweep_hop!`,
identical in every branch/condition to the production function, except:
  - it calls `PathTracking._track_path_segment!` directly (instead of
    going through `FaceTracking.track_dense_path`, which discards
    intermediate bisection points via its own `hop_path[end]`
    convention) so the NUMBER of points it pushes into a fresh, local
    `inner_hop_path` can be inspected -- that count minus 1 is exactly
    the number of times `_track_path_segment!`'s OWN
    poor-accuracy/`is_near_singular` bisection fired for this one hop
    (mirrors `track_dense_path`'s own internal accounting, just not
    discarded here);
  - it increments `outer_bisections` every time `_sweep_hop!`'s own
    residual-based bisection gate fires.
This is NOT a modification of `src/FaceTracking.jl` -- it is a parallel
diagnostic copy, used only for this integration test's introspection
into otherwise-private call counts, exactly like
`scratch_phase5_check.jl`'s own `instrumented_sweep_direction` harness
(which this one subsumes with two additional counters).
"""
function _instrumented_sweep_hop!(
    hop_path::Vector{Vector{Float64}},
    patch::NamedTuple,
    cfg::HomotopyConfig{Float64},
    state::NamedTuple,
    y0::Vector{ComplexF64},
    p0::Float64,
    p1::Float64,
    budget::Base.RefValue{Int},
    tols::NamedTuple,
    inner_bisections::Base.RefValue{Int},
    outer_bisections::Base.RefValue{Int},
)
    expected_rank = length(patch.F_for_tracking.expressions)
    inner_hop_path = Vector{Float64}[]
    inner_budget = Ref(max(cfg.max_path_steps, 1))
    y1 = HomotopyGetsReal._track_path_segment!(
        inner_hop_path, patch.F_for_tracking, state.ph, state.tracker, p0, y0, p1, inner_budget,
        expected_rank, tols.rank_tol64, tols.sv_thresh64, tols.poor_acc_tol64, tols.min_width64,
    )
    inner_bisections[] += length(inner_hop_path) - 1
    seg1 = inner_hop_path[end]
    budget[] -= 1
    x_land, y_land = seg1[2], seg1[3]
    resid = HomotopyGetsReal._residual_at(patch, x_land, y_land, p1, cfg)
    poor = !isfinite(resid) || abs(resid) > cfg.critical_point_tol

    if poor && budget[] > 0 && abs(p1 - p0) > tols.min_width64
        outer_bisections[] += 1
        pmid = (p0 + p1) / 2
        state_mid, ymid = _instrumented_sweep_hop!(hop_path, patch, cfg, state, y0, p0, pmid, budget, tols, inner_bisections, outer_bisections)
        return _instrumented_sweep_hop!(hop_path, patch, cfg, state_mid, ymid, pmid, p1, budget, tols, inner_bisections, outer_bisections)
    end

    push!(hop_path, seg1)
    fx_val, fy_val, _ = HomotopyGetsReal._gradient_at(patch, x_land, y_land, p1, cfg)
    denom = hypot(fx_val, fy_val) * hypot(state.fx0_64, state.fy0_64)
    cos_angle = denom == 0.0 ? 0.0 : (fx_val * state.fx0_64 + fy_val * state.fy0_64) / denom
    if cos_angle < tols.cos_tol64
        anchor_x, anchor_y = HomotopyGetsReal._project_to_slice(patch, x_land, y_land, p1, cfg)
        _, ph_new, tracker_new = HomotopyGetsReal.build_face_tracker(patch, anchor_x, anchor_y, p1, cfg)
        fx0_new, fy0_new, _ = HomotopyGetsReal._gradient_at(patch, anchor_x, anchor_y, p1, cfg)
        new_state = (ph = ph_new, tracker = tracker_new, fx0_64 = Float64(fx0_new), fy0_64 = Float64(fy0_new))
        return new_state, ComplexF64[ComplexF64(anchor_x), ComplexF64(anchor_y)]
    end
    return state, y1
end

"""
    instrumented_sweep_direction(patch, x0, y0, z_start, z_targets, cfg)
        -> (n_rebuilds, n_outer_bisections, n_inner_bisections, dense_path)

One-direction instrumented sweep (reproduces `_sweep_direction`'s own
per-target loop), returning all three fire-counts plus the actual swept
path, so callers can both aggregate statistics AND check the
surface-membership residual on the real output.
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
    inner_bisections = Ref(0)
    outer_bisections = Ref(0)
    dense_path = Vector{Vector{Float64}}(undef, length(z_targets))
    for (i, p_target) in enumerate(z_targets)
        hop_path = Vector{Float64}[]
        budget = Ref(max(cfg.max_path_steps, 1))
        ph_before = state.ph
        state, y_state = _instrumented_sweep_hop!(hop_path, patch, cfg, state, y_state, p_cur, p_target, budget, tols, inner_bisections, outer_bisections)
        state.ph !== ph_before && (n_rebuilds += 1)
        dense_path[i] = hop_path[end]
        p_cur = p_target
    end
    return (n_rebuilds = n_rebuilds, n_outer_bisections = outer_bisections[], n_inner_bisections = inner_bisections[], dense_path = dense_path)
end

"""
    instrumented_track_face(F, patch, edge, z_mid, z_bottom, z_top, cfg) -> NamedTuple

Reproduces `FaceTracking.track_face`'s own loop over
`edge.sampled_points` (one `sweep_face_bidirectional`-equivalent call
per curve sample, in both z-directions), using
`instrumented_sweep_direction` instead, and aggregates fire-counts
across the WHOLE face plus the max |f| residual over every swept point
(the same ground-truth check used throughout Phase 5).
"""
function instrumented_track_face(F, patch, edge, z_mid, z_bottom, z_top, cfg)
    total_rebuilds = 0
    total_outer = 0
    total_inner = 0
    max_resid = 0.0
    for pt in edge.sampled_points
        x0, y0 = HomotopyGetsReal._project_to_slice(patch, pt[1], pt[2], z_mid, cfg)
        n_side = cfg.midslice_sample_density
        targets_down = collect(range(z_mid, Float64(z_bottom); length = n_side + 1))[2:end]
        targets_up = collect(range(z_mid, Float64(z_top); length = n_side + 1))[2:end]
        rd = instrumented_sweep_direction(patch, x0, y0, z_mid, targets_down, cfg)
        ru = instrumented_sweep_direction(patch, x0, y0, z_mid, targets_up, cfg)
        total_rebuilds += rd.n_rebuilds + ru.n_rebuilds
        total_outer += rd.n_outer_bisections + ru.n_outer_bisections
        total_inner += rd.n_inner_bisections + ru.n_inner_bisections
        for p in vcat(rd.dense_path, ru.dense_path)
            r = abs(HomotopyGetsReal._residual_at(patch, p[2], p[3], p[1], cfg))
            max_resid = max(max_resid, r)
        end
    end
    return (n_rebuilds = total_rebuilds, n_outer_bisections = total_outer, n_inner_bisections = total_inner, max_resid = max_resid)
end

println("=" ^ 70)
println("Setup: Taubin heart (x^2+(1.2y)^2+z^2-1)^3 - x^2 z^3 - 0.1(1.2y)^2 z^3 = 0")
println("=" ^ 70)

@var x y z
f = (x^2 + (1.2 * y)^2 + z^2 - 1)^3 - x^2 * z^3 - 0.1 * (1.2 * y)^2 * z^3
F_heart = System([f], variables = [x, y, z]) # z LAST, required convention

# Bounding box tightened from the default (-4,4)^3 based on a direct
# numeric scan of f's real zero set (not HomotopyContinuation -- a plain
# sign-change scan over a grid): the surface's genuine real extent is
# x in [-1.15,1.15], y in [-0.86,0.86], z in [-1.0,1.24]. Loose bbox_z
# margins (-1.3,1.3) keep the two outermost slabs (below the bottom cusp,
# above the top lobe tips) deliberately non-empty-but-small rather than
# exactly tangent, matching the sphere/ellipsoid tests' own precedent of
# including a genuinely-empty outer slab rather than clipping exactly at
# the surface's own extent.
cfg = HomotopyConfig{Float64}(
    bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3),
    edge_sample_density = 8, midslice_sample_density = 8,
)

println()
println("=" ^ 70)
println("1. compute_critical_z_slices: multiple, closely-spaced critical z-values")
println("=" ^ 70)

t_crit = @elapsed z_crits = compute_critical_z_slices(F_heart, cfg)
println("  elapsed: $(round(t_crit; digits = 2))s")
println("  clustered critical z-values -> ", sort(z_crits))

# Also report the RAW (pre-cluster_scalars) critical-point z-coordinates,
# so the clustering decision itself is visible, not just its outcome.
crit_vertices = HomotopyGetsReal.compute_critical_points(F_heart, cfg)
println()
println("  raw critical points (post per-call cluster_vertices dedup, pre cluster_scalars):")
for v in crit_vertices
    println("    id=$(v.id)  coords=", round.(real.(v.coordinates); digits = 6), "  type=$(v.v_type)")
end
z_raw = sort(Float64[real(v.coordinates[3]) for v in crit_vertices])
println("  raw z-values: ", z_raw)
println("  pairwise gaps: ", diff(z_raw))
println("  cfg.vertex_match_tol (the merge threshold): ", cfg.vertex_match_tol)
println()
println("  Interpretation: gaps of ~8.9e-16 and ~0.0 (the two lobe-symmetric pairs, at")
println("  (0,±0.215,1.0648) and (±0.514,0,1.2367)) are correctly MERGED (both far below")
println("  vertex_match_tol=1e-4); gaps of ~0.065 and ~0.172 (genuinely distinct z-slices)")
println("  are correctly KEPT SEPARATE (both far above vertex_match_tol) -- no incorrect")
println("  merging or incorrect separation at this tolerance.")

@test length(z_crits) == 4
@test isapprox(sort(z_crits), [-1.0, 1.0, 1.0647678179140714, 1.2366591700121616]; atol = 1e-6)
println("compute_critical_z_slices: 4 distinct critical z-values confirmed (bottom cusp,")
println("top singular notch, inter-lobe saddle pair, outer lobe-tip pair).")

println()
println("=" ^ 70)
println("2. decompose_3d_surface: full pipeline")
println("=" ^ 70)

t_decomp = @elapsed (all_vertices, all_edges, all_faces, mesh) = decompose_3d_surface(F_heart, cfg)
println("  elapsed: $(round(t_decomp; digits = 2))s")
println("  vertices: ", length(all_vertices), "   edges: ", length(all_edges), "   faces: ", length(all_faces))
println("  z_bounds used: ", sort(unique(vcat([cfg.bbox_z[1]], z_crits, [cfg.bbox_z[2]]))))
@test length(unique(v.id for v in all_vertices)) == length(all_vertices)
@test length(unique(e.id for e in all_edges)) == length(all_edges)

println()
println("=" ^ 70)
println("3. Re-anchoring / bisection stress test near the tips")
println("   (is_near_singular's inner bisection AND patch_transversality_cos_tol's")
println("    outer re-anchoring, per slab)")
println("=" ^ 70)

slab_stats = NamedTuple[]
let
    patch = HomotopyGetsReal.build_patch_system(F_heart)
    z_bounds = sort(unique(vcat([cfg.bbox_z[1]], z_crits, [cfg.bbox_z[2]])))
    for i in 1:(length(z_bounds)-1)
        z_bottom, z_top = z_bounds[i], z_bounds[i+1]
        vertices_2d, edges_2d, z_mid = HomotopyGetsReal._robust_slice_at_z(F_heart, patch, z_bottom, z_top, cfg)
        if isempty(edges_2d)
            println("  slab [$( round(z_bottom;digits=4)), $(round(z_top;digits=4))] (z_mid=$(round(z_mid;digits=4))): empty (no real curve), skipped")
            continue
        end
        n_artificial = count(v -> v.v_type == Artificial, vertices_2d)
        slab_rebuilds = 0
        slab_outer = 0
        slab_inner = 0
        slab_max_resid = 0.0
        for edge in edges_2d
            stats = instrumented_track_face(F_heart, patch, edge, z_mid, z_bottom, z_top, cfg)
            slab_rebuilds += stats.n_rebuilds
            slab_outer += stats.n_outer_bisections
            slab_inner += stats.n_inner_bisections
            slab_max_resid = max(slab_max_resid, stats.max_resid)
        end
        n_hops_total = length(edges_2d) * 2 * cfg.edge_sample_density * cfg.midslice_sample_density
        push!(slab_stats, (z_bottom = z_bottom, z_top = z_top, z_mid = z_mid, n_artificial = n_artificial, max_resid = slab_max_resid))
        z_mid_naive = (z_bottom + z_top) / 2
        retried = !isapprox(z_mid, z_mid_naive; atol = 1e-12)
        flag = retried ? "  <-- _robust_slice_at_z RETRIED (naive midpoint was degenerate, see section 6)" : ""
        println(
            "  slab [$(round(z_bottom;digits=4)), $(round(z_top;digits=4))] (width=$(round(z_top-z_bottom;digits=4)), z_mid used=$(round(z_mid;digits=6)), $(length(edges_2d)) edge(s), $(n_artificial) Artificial vertex/vertices): ",
            "cos_tol re-anchors=$(slab_rebuilds), outer(residual) bisections=$(slab_outer), ",
            "inner(is_near_singular) bisections=$(slab_inner), out of $(n_hops_total) total hops, ",
            "max|f| along swept points=$(slab_max_resid)$(flag)",
        )
    end
end
println()
println("  Interpretation: the two slabs genuinely adjacent to the closely-spaced upper")
println("  tips ([1.0,1.0648] and [1.0648,1.2367], both <0.18 wide) are the intended stress")
println("  test -- both mechanisms fire a moderate, non-zero number of times (dozens of")
println("  re-anchors, hundreds of outer bisections, up to 24 inner is_near_singular")
println("  bisections) and STILL keep max|f| at ~2.4e-6, consistent with the ellipsoid's")
println("  own validated behavior. The [-1.0,1.0] slab, whose naive exact-midpoint z_mid=0")
println("  used to be catastrophically wrong (max|f| up to 1.64, see git history / prior")
println("  report), is now automatically retried by _robust_slice_at_z at a nearby")
println("  perturbed z_mid -- see section 6 for the before/after comparison.")

@test all(s -> s.max_resid < 1e-4, slab_stats)
println("  every non-empty slab (including the previously-degenerate one) now keeps")
println("  max|f| along every swept point below 1e-4.")

println()
println("=" ^ 70)
println("4. |f| residual distribution across the final welded mesh")
println("=" ^ 70)

mesh_pts = GeometryBasics.coordinates(mesh)
mesh_tris = GeometryBasics.faces(mesh)
f_eval(px, py, pz) = (px^2 + (1.2 * py)^2 + pz^2 - 1)^3 - px^2 * pz^3 - 0.1 * (1.2 * py)^2 * pz^3
heart_residuals = [abs(f_eval(p[1], p[2], p[3])) for p in mesh_pts]
sorted_resid = sort(heart_residuals)
println("  |f| residual distribution over all $(length(mesh_pts)) welded mesh vertices:")
println("    min    = ", minimum(heart_residuals))
println("    mean   = ", mean(heart_residuals))
println("    median = ", median(heart_residuals))
println("    p90    = ", sorted_resid[ceil(Int, 0.90 * length(sorted_resid))])
println("    p99    = ", sorted_resid[ceil(Int, 0.99 * length(sorted_resid))])
println("    max    = ", maximum(heart_residuals))
println("    count with residual > 1e-2: ", count(>(1e-2), heart_residuals), " / ", length(heart_residuals))
println("    count with residual > 1e-3: ", count(>(1e-3), heart_residuals), " / ", length(heart_residuals))

# z-partitioned breakdown: isolate the [-1,1] slab (the one _robust_slice_at_z
# had to retry, per section 3) from every other slab, to confirm the fix
# closed the gap everywhere, not just on average.
main_slab_mask = [-1.0 - 1e-3 <= p[3] <= 1.0 + 1e-3 for p in mesh_pts]
resid_main = heart_residuals[main_slab_mask]
resid_rest = heart_residuals[.!main_slab_mask]
println()
println("  z-partitioned breakdown:")
println("    [-1,1] slab (the one _robust_slice_at_z retried, $(count(main_slab_mask)) vertices): ",
        "max=$(maximum(resid_main)), mean=$(mean(resid_main)), median=$(median(resid_main)), ",
        "count>1e-2=$(count(>(1e-2), resid_main))")
println("    every other slab (near the genuine tips, $(count(.!main_slab_mask)) vertices): ",
        "max=$(maximum(resid_rest)), mean=$(mean(resid_rest)), median=$(median(resid_rest)), ",
        "count>1e-2=$(count(>(1e-2), resid_rest))")

@test all(<=(1e-4), heart_residuals)
println("  every welded mesh vertex, INCLUDING the previously-degenerate [-1,1] slab, now")
println("  satisfies f≈0 within 1e-4 -- consistent with every prior Phase 5 test")
println("  (sphere/ellipsoid). The fix closed the gap, not just masked it on average.")

println()
println("=" ^ 70)
println("5. Degenerate-triangle filter + winding-order spot check")
println("=" ^ 70)

println("  mesh: ", length(mesh_pts), " vertices, ", length(mesh_tris), " triangles")
n_degenerate = count(tri -> length(unique((tri[1], tri[2], tri[3]))) < 3, mesh_tris)
println("  degenerate triangles (repeated vertex index) found in final mesh: ", n_degenerate)
@test n_degenerate == 0

patch_heart = HomotopyGetsReal.build_patch_system(F_heart)
n_spot = min(30, length(mesh_tris))
rng_idx = round.(Int, range(1, length(mesh_tris); length = n_spot))
# Alignment predicate is scale-RELATIVE, not a bare `>= 0`: weld_mesh corrects
# winding at T (Float64) precision, but this test recomputes the dot product
# from Float32-rounded mesh coordinates, so an exactly-degenerate triangle
# (normal ~ orthogonal to ∇f) can land epsilon-negative (observed: -1.4e-22,
# sixteen orders below the ~1e-6 magnitude of a genuine flip).
n_aligned = count(idx -> begin
    tri = mesh_tris[idx]
    p1, p2, p3 = mesh_pts[tri[1]], mesh_pts[tri[2]], mesh_pts[tri[3]]
    n = cross(p2 .- p1, p3 .- p1)
    gx, gy, gz = HomotopyGetsReal._gradient_at(patch_heart, Float64(p1[1]), Float64(p1[2]), Float64(p1[3]), cfg)
    g = [gx, gy, gz]
    dot(n, g) >= -1e-10 * norm(n) * norm(g)
end, rng_idx)
worst_dot_sign = minimum(idx -> begin
    tri = mesh_tris[idx]
    p1, p2, p3 = mesh_pts[tri[1]], mesh_pts[tri[2]], mesh_pts[tri[3]]
    n = cross(p2 .- p1, p3 .- p1)
    gx, gy, gz = HomotopyGetsReal._gradient_at(patch_heart, Float64(p1[1]), Float64(p1[2]), Float64(p1[3]), cfg)
    dot(n, [gx, gy, gz])
end, rng_idx)
println("  spot-checked $(n_spot) triangles: $(n_aligned)/$(n_spot) have normal . ∇f >= 0 (expect $(n_spot)/$(n_spot))")
println("  worst (most negative) normal . ∇f among spot-checked triangles: ", worst_dot_sign)
@test n_aligned == n_spot
println("  winding-order convention (normal aligned with +∇f) confirmed on this multi-lobe,")
println("  non-convex, pinch-containing mesh -- not just the sphere/ellipsoid's simpler cases.")

println()
println("=" ^ 70)
println("6. z_mid-degeneracy gap: root cause + _robust_slice_at_z fix confirmation")
println("=" ^ 70)

let
    @var xt yt
    f_at_z0 = HomotopyContinuation.subs(f, [x, y, z] => [xt, yt, 0])
    println("  Root cause (unchanged from the original diagnosis): f(x,y,0) = ", f_at_z0)
    println("  (x^2+1.44y^2-1)^3 expanded:                          ", (xt^2 + 1.44 * yt^2 - 1)^3)
    println("  These are IDENTICAL: f(x,y,0) is an exact PERFECT CUBE of the quadric")
    println("  g(x,y) = x^2+1.44y^2-1 -- a non-reduced curve. The [-1,1] slab's naive")
    println("  z_mid=(z_bottom+z_top)/2 lands EXACTLY there (both z-critical bounds, -1 and")
    println("  +1, are symmetric about 0 by construction -- see the original report for the")
    println("  full argument on why this is structural to the Taubin-heart family, not a")
    println("  coefficient-specific fluke).")

    println()
    println("  BEFORE the fix (naive z_mid=0, calling slice_at_z directly, bypassing the")
    println("  new retry wrapper):")
    vertices_z0, edges_z0 = slice_at_z(F_heart, 0.0, cfg)
    for v in vertices_z0
        origin = get(v.metadata, :origin, nothing)
        println("    vertex id=$(v.id) coords=$(round.(real.(v.coordinates); digits=4)) type=$(v.v_type) metadata[:origin]=$(origin)")
    end
    n_fallback = count(v -> v.v_type == Artificial && get(v.metadata, :origin, nothing) == :endpoint_fallback, vertices_z0)
    println("  -> $(n_fallback) vertex/vertices tagged :endpoint_fallback (the precise signal")
    println("  _robust_slice_at_z keys off, per Topology._resolve_endpoint's new tag --")
    println("  distinct from any ordinary cluster_vertices-merge-driven Artificial vertex,")
    println("  which would NOT carry this tag).")

    println()
    println("  The vertex-type gate alone is NOT sufficient here: retrying only on")
    println("  :endpoint_fallback+Singular co-occurrence lands on z_mid=0.02 (topologically")
    println("  clean -- 0 Artificial/Singular vertices) but that candidate is STILL deep in a")
    println("  numerically near-degenerate neighborhood (gradient ~ z^2 near the true")
    println("  repeated-factor plane), which wrecks the downstream sweep. Confirmed directly:")
    patch_heart = HomotopyGetsReal.build_patch_system(F_heart)
    vertices_002, edges_002 = slice_at_z(F_heart, 0.02, cfg)
    face_002 = HomotopyGetsReal.track_face(F_heart, patch_heart, edges_002[1], 0.02, -1.0, 1.0, 1, cfg)
    resid_002 = maximum(abs(HomotopyGetsReal._residual_at(patch_heart, face_002.mesh_vertices[i,1], face_002.mesh_vertices[i,2], face_002.mesh_vertices[i,3], cfg)) for i in 1:size(face_002.mesh_vertices,1))
    println("    z_mid=0.02 (vertex-type-clean, gate REJECTED by gradient check): downstream")
    println("    sweep max|f|=$(resid_002) -- catastrophic despite clean topology.")
    @test resid_002 > 1e-2

    println()
    println("  AFTER the full two-part fix (_robust_slice_at_z(F, patch, -1.0, 1.0, cfg), the")
    println("  actual call decompose_3d_surface now makes for this slab):")
    vertices_fixed, edges_fixed, z_mid_fixed = HomotopyGetsReal._robust_slice_at_z(F_heart, patch_heart, -1.0, 1.0, cfg)
    println("    z_mid used: ", z_mid_fixed, "  (perturbed away from both the naive 0.0 AND the")
    println("    vertex-type-clean-but-gradient-suspect 0.02)")
    for v in vertices_fixed
        println("    vertex id=$(v.id) coords=$(round.(real.(v.coordinates); digits=4)) type=$(v.v_type)")
    end
    n_fallback_fixed = count(v -> v.v_type == Artificial && get(v.metadata, :origin, nothing) == :endpoint_fallback, vertices_fixed)
    face_fixed = HomotopyGetsReal.track_face(F_heart, patch_heart, edges_fixed[1], z_mid_fixed, -1.0, 1.0, 1, cfg)
    resid_fixed = maximum(abs(HomotopyGetsReal._residual_at(patch_heart, face_fixed.mesh_vertices[i,1], face_fixed.mesh_vertices[i,2], face_fixed.mesh_vertices[i,3], cfg)) for i in 1:size(face_fixed.mesh_vertices,1))
    println("  -> $(n_fallback_fixed) :endpoint_fallback vertices, downstream sweep max|f|=$(resid_fixed)")
    println("  -- both topologically clean AND numerically well-conditioned.")
    @test n_fallback_fixed == 0
    @test length(vertices_fixed) == 2
    @test length(edges_fixed) == 2
    @test resid_fixed < 1e-4
end

println()
println("  Retry cost actually incurred across the WHOLE surface (from section 3's per-slab")
println("  report above): see the \"z_mid used\" column -- any slab whose reported z_mid")
println("  differs from its naive (z_bottom+z_top)/2 needed at least one retry of either")
println("  gate. No slab came close to exhausting cfg.max_z_mid_retries=$(cfg.max_z_mid_retries).")

println()
println("=" ^ 70)
println("7. Phase 9a: naked-edge characterization baseline (the defect 9b must close)")
println("=" ^ 70)

# Verified deterministic across repeated decomposes and separate Julia
# sessions (2026-07). These cracks are the CHARACTERIZED CURRENT DEFECT of
# tolerance-only welding across slab boundaries: adjacent slabs' faces land on
# the same crit-slice curve at different sample points, so weld_mesh never
# merges them. Phase 9b's incidence-based stitching must drive these to zero;
# if Phase 9a's own changes shift ANY of these numbers, that is an immediate
# regression signal.
naked = HomotopyGetsReal._naked_mesh_edges(mesh)
count_naked_near(zc) = count(((i, j),) -> abs((mesh_pts[i][3] + mesh_pts[j][3]) / 2 - zc) < 0.01, naked)
println("  naked edges: $(length(naked)) total | z=-1: $(count_naked_near(-1.0))  z=1: $(count_naked_near(1.0))  z=1.0648: $(count_naked_near(1.0648))  z=1.2367: $(count_naked_near(1.2367))")
@test length(naked) == 188
@test count_naked_near(-1.0) == 10
@test count_naked_near(1.0) == 70
@test count_naked_near(1.0648) == 84
@test count_naked_near(1.2367) == 24

println()
println("=" ^ 70)
println("8. Phase 9a: fixed-axis incidence run -- crit-slice census, coverage,")
println("   and the assignment-distance distributions that will drive 9b's thresholds")
println("=" ^ 70)

t_inc = @elapsed (vi, ei, fi, mi, inc) = decompose_3d_surface(F_heart, cfg; incidence = true)
println("  [$(round(t_inc; digits = 1)) s]  $(length(fi)) faces, $(length(inc.crit_slices)) crit-slices")

# z_bounds = [-1.3, -1, 1, 1.0648, 1.2367, 1.3] -> interior boundaries j = 2..5.
@test length(inc.crit_slices) == 4
cs_by_z = Dict(round(cs.z; digits = 4) => cs for cs in inc.crit_slices)
@test length(cs_by_z[-1.0].edges) == 0          # cusp: point-slice
@test cs_by_z[-1.0].is_degenerate
@test length(cs_by_z[1.0].edges) == 4           # singular notch: decomposes
@test length(cs_by_z[1.0648].edges) == 12       # saddle pair: nodal curve
@test length(cs_by_z[1.2367].edges) == 0        # lobe tips: point-slice
@test cs_by_z[1.2367].is_degenerate

# id namespace (Phase 10, Stage 2): edges stay globally unique (unchanged
# positional scheme); vertex ids may now legitimately collide across
# containers via the shared VertexRegistry -- this fixture is the one that
# actually EXERCISES that (z=1.0 and z=1.0648 crit-slices carry real
# vertices, unlike the sphere/ellipsoid controls in test_incidence.jl), so
# any collision here must be a genuine geometric match, not an accident.
all_registered_vertices = vcat(vi, [v for cs in inc.crit_slices for v in cs.vertices])
by_id = Dict{Int,Vector{Vector{ComplexF64}}}()
for v in all_registered_vertices
    push!(get!(by_id, v.id, Vector{ComplexF64}[]), v.coordinates)
end
collision_ok = all(coords -> all(c -> norm(c .- coords[1]) <= cfg.vertex_match_tol, coords), values(by_id))
n_colliding_ids = count(coords -> length(coords) > 1, values(by_id))
@test collision_ok
known_vertex_ids = Set(keys(by_id))

all_eids = vcat([e.id for e in ei], [e.id for cs in inc.crit_slices for e in cs.edges])
@test length(unique(all_eids)) == length(all_eids)

all_endpoint_ids = vcat(
    [e.left_vertex_id for e in ei], [e.right_vertex_id for e in ei],
    [e.left_vertex_id for cs in inc.crit_slices for e in cs.edges],
    [e.right_vertex_id for cs in inc.crit_slices for e in cs.edges],
)
@test all(id -> id in known_vertex_ids, all_endpoint_ids)
println("  vertex ids: $(n_colliding_ids) legitimate cross-container collision(s) (registry ",
        "merges), all geometrically verified; edge ids unique; every endpoint resolves.")

# coverage: every non-bbox landing is assigned to SOMETHING (pure measurement,
# no distance thresholds in 9a -- see ColumnLanding's docstring)
all_non_bbox = [l for d in (inc.column_landings_bottom, inc.column_landings_top) for ls in values(d) for l in ls if l.kind != :bbox]
@test !isempty(all_non_bbox)
@test all(l -> l.kind != :none, all_non_bbox)
println("  coverage: $(length(all_non_bbox)) non-bbox landings, 100% assigned (0 :none)")

# THE 9b DESIGN DATA: per-boundary assignment-distance distributions.
z_bounds_h = sort(unique(vcat([cfg.bbox_z[1]], filter(zc -> cfg.bbox_z[1] <= zc <= cfg.bbox_z[2], z_crits), [cfg.bbox_z[2]])))
face_slab(f) = findfirst(i -> z_bounds_h[i] < f.mid_slice_z < z_bounds_h[i+1], 1:(length(z_bounds_h) - 1))
println("  assignment-distance distributions (median / p90 / max) per interior boundary:")
for j in 2:(length(z_bounds_h) - 1)
    dists_bottom = Float64[l.dist for f in fi if face_slab(f) == j for l in inc.column_landings_bottom[f.id]]
    dists_top = Float64[l.dist for f in fi if face_slab(f) == j - 1 for l in inc.column_landings_top[f.id]]
    dd = sort(vcat(dists_bottom, dists_top))
    isempty(dd) && continue
    kinds = Dict{Symbol,Int}()
    for f in fi, (slabidx, d) in ((j, inc.column_landings_bottom), (j - 1, inc.column_landings_top))
        face_slab(f) == slabidx || continue
        for l in d[f.id]
            kinds[l.kind] = get(kinds, l.kind, 0) + 1
        end
    end
    println("    boundary z=$(round(z_bounds_h[j]; digits = 4)): n=$(length(dd))  median=$(round(dd[cld(end, 2)]; sigdigits = 3))  p90=$(round(dd[ceil(Int, 0.9length(dd))]; sigdigits = 3))  max=$(round(dd[end]; sigdigits = 3))  kinds=$(kinds)")
end

println()
println("=" ^ 70)
println("8b. Phase 9c: stitched-mesh naked-edge closure (the actual, honest outcome)")
println("=" ^ 70)

# weld_mesh's incidence path (exercised automatically since `mi` above came
# from incidence = true) layers THREE mechanisms: Phase 9b's snap-unification
# (_snap_boundary_points!, closes fold/point-type boundaries COMPLETELY --
# cusp z=-1, lobe tips z=1.2367 -- matching the sphere/ellipsoid pole-
# convergence controls exactly, 0 naked edges, deterministic), then Phase
# 9c's coordinated loft (_identify_edge_runs + _append_loft_triangles!:
# full-coverage ribbon triangulation against a crit-slice edge's COMPLETE
# sample sequence, not a snap-derived subset, plus seam-capping) and
# _chained_edge_polylines + _split_t_junctions (T-junction fixing extended to
# cover junction points SHARED between adjacent crit-slice edges, not just
# within one edge's own polyline).
#
# Investigated directly (2026-07) whether a higher-fidelity, higher-risk
# alternative -- BertiniReal-style targeted homotopy tracking toward known
# crit-slice vertices ("Option A"), replacing free-sweep track_face hops --
# would do meaningfully better than this loft approach ("Option B"):
# measured precision showed no mechanism for Option A to beat DIRECT reuse
# of already-Newton-polished crit-slice coordinates (parity at best, since
# it necessarily adds tracker predictor-corrector error on top of converging
# to the SAME target); not pursued.
#
# Measured on the multi-edge, multi-face EDGE-type boundaries (singular
# notch z=1.0, saddle pair z=1.0648) -- the ONLY boundaries this reduces,
# since fold/point-type boundaries were already fully closed by Phase 9b
# alone: total naked edges 58 (Phase 9b only) -> 31-35 across 4 repeated
# fixed-axis decomposes (2026-07) -- a genuine, substantial reduction, NOT a
# full closure. The run-to-run spread (31/33/34/35, z=1.0 in [16,19], z=1.0648
# in [15,17]) is the SAME already-documented cross-process HC.jl solver
# jitter this file's rotated-seed-1 section exists to characterize (affects
# crit-slice sample coordinates by tiny amounts, which this file's finer
# per-column confidence/bridge-DP decisions are visibly more sensitive to
# than Phase 9b's coarser snap ever was) -- hence loose, not exact, bounds
# below for the multi-edge boundaries specifically. Diagnosed but NOT fully
# resolved: roughly a third of the residual is genuine cross-edge-junction
# cases the chaining did not close, a third seam artifacts at run
# boundaries the capping did not fully reach, and a third UNDIAGNOSED after
# three separate investigation attempts -- confirmed not a resolution
# artifact (matching Phase 9b's own finding for edge_sample_density). Full
# closure is DEFERRED as a future, separately-scoped phase, same treatment
# as full singular-curve decomposition.
naked_stitched = HomotopyGetsReal._naked_mesh_edges(mi)
mi_pts = GeometryBasics.coordinates(mi)
count_stitched_near(zc) = count(((i, j),) -> abs((mi_pts[i][3] + mi_pts[j][3]) / 2 - zc) < 0.01, naked_stitched)
println("  stitched-mesh naked edges: $(length(naked_stitched)) total | z=-1: $(count_stitched_near(-1.0))  z=1: $(count_stitched_near(1.0))  z=1.0648: $(count_stitched_near(1.0648))  z=1.2367: $(count_stitched_near(1.2367))")
println("  (unstitched baseline was 188 = 10/70/84/24, Phase 9b-only baseline was 58 = 0/26/32/0 -- see section 7)")
@test count_stitched_near(-1.0) == 0        # fold-type: deterministic, exact (Phase 9b alone already closes this)
@test count_stitched_near(1.2367) == 0      # fold-type: deterministic, exact (Phase 9b alone already closes this)
@test length(naked_stitched) < 58           # genuine improvement over Phase 9b -- not exact, see jitter note above
@test length(naked_stitched) <= 40          # regression guard around the measured 31-35 range
println("  fold/point-type boundaries (cusp, tips) FULLY closed, matching sphere/ellipsoid")
println("  pole-convergence controls exactly. Multi-face edge-type boundaries substantially")
println("  reduced (58 -> 31-35) but not closed -- see comment above for the diagnosis and")
println("  the deferred fix; loose bounds here reflect measured HC.jl cross-process jitter,")
println("  not a weaker guarantee about the mechanism itself.")

println()
println("=" ^ 70)
println("9. Phase 9a: rotated seed-1 incidence probe (evidence-gathering, loose asserts)")
println("=" ^ 70)

Q1 = random_orthogonal_matrix(Float64, 3; rng = Xoshiro(1))
t_rot_inc = @elapsed (vr, er, fr, mr, incr) = decompose_3d_surface(F_heart, cfg; projection = Q1, incidence = true)
@test incr isa SurfaceIncidence{Float64}
@test !isempty(incr.crit_slices)
n_none = count(l -> l.kind == :none, [l for d in (incr.column_landings_bottom, incr.column_landings_top) for ls in values(d) for l in ls])
n_all = length([l for d in (incr.column_landings_bottom, incr.column_landings_top) for ls in values(d) for l in ls])
println("  [$(round(t_rot_inc; digits = 1)) s]  $(length(fr)) faces, $(length(incr.crit_slices)) crit-slices, $(n_none)/$(n_all) landings :none")
println("  per-boundary census + distance stats (rotated chart -- the 9b uncertainty probe):")
for cs in incr.crit_slices
    dd = Float64[]
    for d in (incr.column_landings_bottom, incr.column_landings_top), ls in values(d), l in ls
        # attribute landings to this boundary by matching the assigned cell's namespace
        if (l.kind == :edge && any(e -> e.id == l.id, cs.edges)) ||
           (l.kind == :crit_slice_vertex && any(v -> v.id == l.id, cs.vertices))
            push!(dd, l.dist)
        end
    end
    sort!(dd)
    stats = isempty(dd) ? "no cell-assigned landings" :
        "n=$(length(dd)) median=$(round(dd[cld(end, 2)]; sigdigits = 3)) p90=$(round(dd[ceil(Int, 0.9length(dd))]; sigdigits = 3)) max=$(round(dd[end]; sigdigits = 3))"
    nfall = count(v -> v.v_type == HomotopyGetsReal.Artificial && get(v.metadata, :origin, nothing) == :endpoint_fallback, cs.vertices)
    println("    boundary j=$(cs.boundary_index) z=$(round(cs.z; digits = 4)): $(length(cs.edges)) edges, $(length(cs.vertices)) verts ($(nfall) fallback), degenerate=$(cs.is_degenerate) | $stats")
end

end
