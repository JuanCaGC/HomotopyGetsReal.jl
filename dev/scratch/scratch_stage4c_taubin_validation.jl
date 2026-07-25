# dev/scratch/scratch_stage4c_taubin_validation.jl
#
# Stage 4c required validation (2026-07): run decompose_3d_surface on the
# Taubin heart fixture (test/test_taubin.jl's own equation/cfg, reused
# verbatim) WITHOUT deflate first to enumerate today's real Singular
# vertices, then WITH deflate=true (both incidence=true, so both the
# slice-level and surface-level diagnostics run) to report per-vertex
# isosingular data at both levels, explicitly labeled per the approved
# reporting format. Not committed as a permanent test -- investigation
# script only, matching this project's established pattern.

using HomotopyContinuation
using LinearAlgebra
using GeometryBasics
using GLMakie
using HomotopyGetsReal

@var x y z
f = (x^2 + (1.2 * y)^2 + z^2 - 1)^3 - x^2 * z^3 - 0.1 * (1.2 * y)^2 * z^3
F_heart = System([f], variables = [x, y, z])

cfg = HomotopyConfig{Float64}(
    bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3),
    edge_sample_density = 8, midslice_sample_density = 8,
)

println("="^70)
println("Step 1: baseline (no deflate) -- enumerate today's real Singular vertices")
println("="^70)
t0 = @elapsed (v0, e0, f0, m0, inc0) = decompose_3d_surface(F_heart, cfg; incidence = true)
println("elapsed: $(round(t0; digits=2))s")
sing0 = filter(v -> v.v_type == Singular, v0)
println("total vertices: ", length(v0), "   Singular: ", length(sing0))
for v in sing0
    println("  id=$(v.id)  coords=", round.(real.(v.coordinates); digits=6))
end

# Gate instrumentation (2026-07, this pass): the user asked an explicit, precise
# question -- does the new _deflation_applicable gate exclude EVERY Faug-Singular
# candidate encountered in this fixture (meaning resolve_isosingular_dimension never
# actually fires anywhere in this run), or does it let some through? The final
# top-level Singular COUNT alone cannot answer this (candidates seen mid-pipeline,
# e.g. inside a single crit-slice's own compute_critical_points call, are not
# necessarily what survives to the final merged/deduped surface vertex list). To get
# a real answer, monkey-patch _deflation_applicable in this session only (not the
# source file) to count every call and every true/false outcome, covering every
# call site that runs under deflate=true (crit-slices via _decompose_crit_slice ->
# slice_at_z -> decompose_1d_curve -> compute_critical_points, the robust-slice
# retry path, and the direct surface-level _surface_critical_vertices call below).
const _gate_calls = Ref(0)
const _gate_true = Ref(0)
function HomotopyGetsReal._deflation_applicable(
    F_original::HomotopyGetsReal.System,
    x0::AbstractVector,
    cfg::HomotopyGetsReal.HomotopyConfig{T};
    expected_rank::Int = length(F_original.variables),
) where {T<:AbstractFloat}
    _gate_calls[] += 1
    d = HomotopyGetsReal.estimate_corank(F_original, x0, cfg; expected_rank = expected_rank)
    result = if d == 0
        true
    else
        minor_size = expected_rank - d + 1
        minor_size <= length(F_original.expressions) && minor_size <= length(F_original.variables)
    end
    result && (_gate_true[] += 1)
    return result
end

println()
println("="^70)
println("Step 2: deflate=true, incidence=true -- full diagnostic run")
println("="^70)
t1 = @elapsed (v1, e1, f1, m1, inc1) = decompose_3d_surface(F_heart, cfg; incidence = true, deflate = true)
println("elapsed: $(round(t1; digits=2))s   (baseline was $(round(t0; digits=2))s)")
sing1 = filter(v -> v.v_type == Singular, v1)
println("total vertices: ", length(v1), "   Singular: ", length(sing1))

println()
println("Slice-level vertex set unchanged by deflate=true (sanity check): ",
    length(sing0) == length(sing1) && sort(collect(v.id for v in sing0)) == sort(collect(v.id for v in sing1)))

println()
println("Gate instrumentation so far (Step 2's decompose_3d_surface call only):")
println("  _deflation_applicable called: $(_gate_calls[]) time(s)  (== every Faug-Singular candidate encountered anywhere in the pipeline)")
println("  of those, applicable=true (resolve_isosingular_dimension actually invoked): $(_gate_true[])")

println()
println("="^70)
println("Per-vertex table: SLICE-LEVEL diagnostic (deflating against F_2d, the z-slice curve)")
println("="^70)
red_flags = String[]
for v in sing1
    coords = round.(real.(v.coordinates); digits=6)
    if haskey(v.metadata, :isosingular_verdict)
        verdict = v.metadata[:isosingular_verdict]
        dim = v.metadata[:isosingular_dimension]
        seq = v.metadata[:isosingular_deflation_sequence]
        rounds = length(seq) - 1
        println("  id=$(v.id)  coords=$coords")
        println("    slice-level:   verdict=$verdict  dim=$dim  rounds=$rounds  corank_seq=$seq")
        if verdict != Resolved
            push!(red_flags, "id=$(v.id) slice-level verdict=$verdict (expected Resolved)")
        elseif !(rounds in 1:2)
            push!(red_flags, "id=$(v.id) slice-level rounds=$rounds outside established [1,2] (soft flag)")
        end
    else
        println("  id=$(v.id)  coords=$coords")
        println("    slice-level:   NO isosingular metadata (did not go through a deflate=true decompose_1d_curve call this run)")
    end
end

println()
println("="^70)
println("Per-vertex table: SURFACE-LEVEL diagnostic (deflating against full 3-var F)")
println("="^70)
println("(from fold_anchors / _surface_critical_vertices, cross-referenced to the same physical points by coordinate proximity)")
gate_calls_before_surface = _gate_calls[]
gate_true_before_surface = _gate_true[]
fold_anchors_deflated = HomotopyGetsReal._surface_critical_vertices(F_heart, cfg; deflate = true)
surface_singular = filter(v -> v.v_type == Singular, fold_anchors_deflated)
println("surface-level total vertices: $(length(fold_anchors_deflated))   Singular: $(length(surface_singular))")
println("gate calls during this surface-level call: $(_gate_calls[] - gate_calls_before_surface)   applicable=true: $(_gate_true[] - gate_true_before_surface)")
for fv in surface_singular
    coords = round.(real.(fv.coordinates); digits=6)
    # match to a slice-level vertex by coordinate proximity, for side-by-side reporting
    match = findfirst(v -> norm(real.(v.coordinates) .- real.(fv.coordinates)) < cfg.vertex_match_tol * 10, sing1)
    match_label = match === nothing ? "(no slice-level match within tolerance)" : "matches slice-level id=$(sing1[match].id)"
    if haskey(fv.metadata, :isosingular_verdict)
        verdict = fv.metadata[:isosingular_verdict]
        dim = fv.metadata[:isosingular_dimension]
        seq = fv.metadata[:isosingular_deflation_sequence]
        rounds = length(seq) - 1
        println("  coords=$coords  $match_label")
        println("    surface-level: verdict=$verdict  dim=$dim  rounds=$rounds  corank_seq=$seq")
        if verdict != Resolved
            push!(red_flags, "surface-level at $coords verdict=$verdict (expected Resolved)")
        elseif !(rounds in 1:2)
            push!(red_flags, "surface-level at $coords rounds=$rounds outside established [1,2] (soft flag)")
        end
    else
        println("  coords=$coords  $match_label")
        println("    surface-level: gated out (not applicable), no isosingular metadata")
    end
end

println()
println("="^70)
println("EXPLICIT FINDING: did resolve_isosingular_dimension actually fire anywhere in this run?")
println("="^70)
if _gate_calls[] > 0 && _gate_true[] == 0
    println("  NO. _deflation_applicable was called $(_gate_calls[]) time(s) across this entire run")
    println("  (crit-slices + robust-slice retries + the direct surface-level call), and EVERY")
    println("  single call returned false. Taubin's real singular structure, under this cfg, does")
    println("  not include any point where deflation actually fires -- the crash fix is confirmed")
    println("  (no ArgumentError, clean completion), but this run provides NO cost/behavior")
    println("  evidence for genuine surface singularities going through resolve_isosingular_dimension.")
    println("  This is NOT the same finding as 'validated real deflation behavior on Taubin and it")
    println("  was clean' -- it is 'the gate correctly filtered out every candidate this fixture has.'")
elseif _gate_calls[] == 0
    println("  N/A. _deflation_applicable was never even called -- this run encountered zero")
    println("  Faug-Singular candidates anywhere (slice-level or surface-level), so there was")
    println("  nothing for the gate to filter in the first place. Same conclusion as above: no")
    println("  cost/behavior evidence for genuine singularities from this run.")
else
    println("  YES. $(_gate_true[]) of $(_gate_calls[]) gate call(s) returned true, meaning")
    println("  resolve_isosingular_dimension genuinely ran on at least one real candidate in this")
    println("  fixture -- see the per-vertex tables above for the actual verdicts/rounds.")
end

println()
println("="^70)
println("Red-flag summary (thresholds defined BEFORE this run, per the approved proposal)")
println("="^70)
if isempty(red_flags)
    println("  none -- every verdict Resolved, every round/attempt count within established range")
else
    for r in red_flags
        println("  FLAG: ", r)
    end
end
