# dev/scratch/scratch_historical_curves_check.jl
#
# Measurement-only investigation (2026-07): run three historical plane
# curves (from an earlier thesis draft, with BertiniReal reference numbers
# already on record) through the CURRENT decompose_1d_curve pipeline and
# report vertex/edge/critical-point counts, timing, error/degenerate-case
# observations, and a visual sanity check -- for manual comparison against
# the BertiniReal numbers, not an automated pass/fail judgment.

using HomotopyContinuation
using LinearAlgebra
using GeometryBasics
using GLMakie
using HomotopyGetsReal

"""
    _scan_extent(fexpr, xv, yv; lo=-6.0, hi=6.0, n=400) -> (xrange, yrange)

Cheap, non-HC numeric sign-change scan (same style as the Taubin heart's own
bbox-tightening comment: "a plain sign-change scan over a grid, not
HomotopyContinuation") over a coarse grid, to check whether the curve's real
extent comfortably fits the default HomotopyConfig{Float64}() bbox of
(-4,4)x(-4,4) before running anything.
"""
function _scan_extent(fexpr, xv, yv; lo = -6.0, hi = 6.0, n = 400)
    xs = range(lo, hi; length = n)
    ys = range(lo, hi; length = n)
    xmin, xmax, ymin, ymax = Inf, -Inf, Inf, -Inf
    found = false
    prev_row = fill(NaN, n)
    for (j, yv0) in enumerate(ys)
        for (i, xv0) in enumerate(xs)
            val = HomotopyContinuation.evaluate(fexpr, [xv, yv] => [xv0, yv0])
            v = real(val)
            if i > 1 && sign(v) != sign(prev_row[i-1]) && !isnan(prev_row[i-1])
                found = true
                xmin, xmax = min(xmin, xv0), max(xmax, xv0)
                ymin, ymax = min(ymin, yv0), max(ymax, yv0)
            end
            prev_row[i] = v
        end
    end
    # Also scan columns (catch horizontal sign changes the row scan alone can miss).
    prev_col = fill(NaN, n)
    for (i, xv0) in enumerate(xs)
        for (j, yv0) in enumerate(ys)
            val = HomotopyContinuation.evaluate(fexpr, [xv, yv] => [xv0, yv0])
            v = real(val)
            if j > 1 && sign(v) != sign(prev_col[j-1]) && !isnan(prev_col[j-1])
                found = true
                xmin, xmax = min(xmin, xv0), max(xmax, xv0)
                ymin, ymax = min(ymin, yv0), max(ymax, yv0)
            end
            prev_col[j] = v
        end
    end
    return found, (xmin, xmax), (ymin, ymax)
end

println("=" ^ 70)
println("Pre-scan: real extent of each curve (plain grid sign-scan, no HC)")
println("=" ^ 70)

@var xv yv
f1 = (xv^2 + yv^2 - 1)^3 - 27 * xv^2 * yv^2
f2 = (xv^3 - xv * yv^2 + yv + 1)^2 * (xv^2 + yv^2 - 1) + yv^2 - 5
f3 = xv^8 - 28 * xv^6 * yv^2 + 70 * xv^4 * yv^4 - 28 * xv^2 * yv^6 + yv^8 + 15 * xv^4 * yv^2 - 15 * xv^2 * yv^4

for (label, fe) in (("Curve 1", f1), ("Curve 2", f2), ("Curve 3", f3))
    found, (xmin, xmax), (ymin, ymax) = _scan_extent(fe, xv, yv)
    println("  $label: found=$found  x in [$xmin, $xmax]  y in [$ymin, $ymax]")
end
println()
println("  NOTE on Curve 3: separately re-scanned at +-6, +-10, +-20 -- the detected")
println("  extent kept growing with the scan range every time (6.0 -> 10.0 -> 20.0),")
println("  unlike curves 1/2 (stable at ~2.4 even scanned out to +-10). Curve 3's")
println("  degree-8 leading term x^8-28x^6y^2+70x^4y^4-28x^2y^6+y^8 is EXACTLY")
println("  Re((x+iy)^8) = r^8*cos(8*theta) -- its own zero set is 16 rays through the")
println("  origin at angles solving cos(8*theta)=0, valid at EVERY radius. The lower-order")
println("  perturbation (15x^4y^2-15x^2y^4, degree 6) dominates near the origin but is")
println("  negligible at large r, so Curve 3's true real zero set is UNBOUNDED --")
println("  asymptotic to those 16 rays at infinity, not a compact curve. Any bbox choice")
println("  (ours or BertiniReal's, unknown) only captures a finite window of an intrinsically")
println("  unbounded curve, so vertex/edge counts for Curve 3 specifically are NOT just a")
println("  density/precision comparison -- they also depend on which finite window was used.")

println()
println("=" ^ 70)
println("Running each curve through decompose_1d_curve")
println("=" ^ 70)

function _report(label, F, cfg; used_default::Bool)
    println()
    println("-" ^ 70)
    println(label, used_default ? " (default HomotopyConfig{Float64}())" : " (CUSTOM cfg -- see note)")
    println("-" ^ 70)

    t_cold = @elapsed (vertices, edges) = decompose_1d_curve(F, cfg)
    t_warm = @elapsed decompose_1d_curve(F, cfg)

    n_critical = count(v -> v.v_type == HomotopyGetsReal.Critical, vertices)
    n_singular = count(v -> v.v_type == HomotopyGetsReal.Singular, vertices)
    n_boundary = count(v -> v.v_type == HomotopyGetsReal.Boundary, vertices)
    n_artificial = count(v -> v.v_type == HomotopyGetsReal.Artificial, vertices)

    println("  wall-clock: first-call=$(round(t_cold; digits=3))s  warm=$(round(t_warm; digits=3))s")
    println("  critical points (Critical+Singular, post-clustering): $(n_critical + n_singular)  (Critical=$(n_critical), Singular=$(n_singular))")
    println("  vertices total: $(length(vertices))  (Critical=$(n_critical), Boundary=$(n_boundary), Singular=$(n_singular), Artificial=$(n_artificial))")
    println("  edges: $(length(edges))")

    n_degenerate_edges = count(e -> length(unique(e.sampled_points)) < 2, edges)
    n_singular_edges = count(e -> e.is_singular, edges)
    println("  degenerate edges (< 2 distinct sample points): $(n_degenerate_edges)")
    println("  edges flagged is_singular: $(n_singular_edges)")

    fallback_vertices = count(v -> get(v.metadata, :origin, nothing) == :endpoint_fallback, vertices)
    println("  :endpoint_fallback artificial vertices: $(fallback_vertices)")

    return vertices, edges, (t_cold = t_cold, t_warm = t_warm)
end

cfg_default = HomotopyConfig{Float64}()
println("Default HomotopyConfig{Float64}() bbox: x=$(cfg_default.bbox_x) y=$(cfg_default.bbox_y)")

F1 = System([f1], variables = [xv, yv])
v1, e1, t1 = _report("Curve 1: (x^2+y^2-1)^3 - 27x^2y^2 = 0", F1, cfg_default; used_default = true)

F2 = System([f2], variables = [xv, yv])
v2, e2, t2 = _report("Curve 2: (x^3-xy^2+y+1)^2 (x^2+y^2-1) + y^2 - 5 = 0", F2, cfg_default; used_default = true)

F3 = System([f3], variables = [xv, yv])
v3, e3, t3 = _report("Curve 3: x^8 -28x^6y^2 +70x^4y^4 -28x^2y^6 +y^8 +15x^4y^2 -15x^2y^4 = 0", F3, cfg_default; used_default = true)

println()
println("=" ^ 70)
println("Curve 3 bbox-sensitivity check (given the confirmed unboundedness above):")
println("re-running with a LARGER bbox to show how much the counts move.")
println("=" ^ 70)
cfg_wide = HomotopyConfig{Float64}(bbox_x = (-8.0, 8.0), bbox_y = (-8.0, 8.0))
v3w, e3w, t3w = _report("Curve 3 (bbox_x=bbox_y=(-8,8), CUSTOM, for sensitivity comparison only)", F3, cfg_wide; used_default = false)

println()
println("=" ^ 70)
println("SUMMARY TABLE (current package vs. BertiniReal reference, for manual comparison)")
println("=" ^ 70)
println("  Curve 1: ours crit=$(count(v->v.v_type in (HomotopyGetsReal.Critical,HomotopyGetsReal.Singular), v1)) vertices=$(length(v1)) edges=$(length(e1))   |   BertiniReal crit=6 vertices=32 edges=20")
println("  Curve 2: ours crit=$(count(v->v.v_type in (HomotopyGetsReal.Critical,HomotopyGetsReal.Singular), v2)) vertices=$(length(v2)) edges=$(length(e2))   |   BertiniReal crit=8 vertices=44 edges=28")
println("  Curve 3: ours crit=$(count(v->v.v_type in (HomotopyGetsReal.Critical,HomotopyGetsReal.Singular), v3)) vertices=$(length(v3)) edges=$(length(e3)))  [default bbox]   |   BertiniReal crit=25 vertices=223 edges=125")
println("  Curve 3: ours crit=$(count(v->v.v_type in (HomotopyGetsReal.Critical,HomotopyGetsReal.Singular), v3w)) vertices=$(length(v3w)) edges=$(length(e3w))  [bbox=(-8,8), unbounded curve -- see note above]")

println()
println("=" ^ 70)
println("Visual sanity check: rendering each curve")
println("=" ^ 70)
mkpath(joinpath(@__DIR__, "renders_historical"))
for (label, F, v, e) in (("curve1", F1, v1, e1), ("curve2", F2, v2, e2), ("curve3", F3, v3, e3), ("curve3_wide_bbox", F3, v3w, e3w))
    fig = plot_curve_decomposition(v, e)
    path = joinpath(@__DIR__, "renders_historical", "$(label).png")
    GLMakie.save(path, fig)
    println("  saved $path")
end
