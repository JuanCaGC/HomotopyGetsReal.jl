# dev/scratch/scratch_stage4c_historical_curves.jl
#
# Stage 4c required validation (2026-07), historical-curve leg: re-run the
# three historical plane curves (same equations/cfg as
# scratch_historical_curves_check.jl) through decompose_1d_curve with
# deflate=true, now that the _deflation_applicable gate is in place, and
# report per-vertex isosingular data plus the same red-flag summary format
# used for the Taubin validation. Investigation script only, not committed
# as a permanent test, matching this project's established pattern.

using HomotopyContinuation
using LinearAlgebra
using HomotopyGetsReal

@var xv yv
f1 = (xv^2 + yv^2 - 1)^3 + 27 * xv^2 * yv^2
f2 = (xv^3 - xv * yv^2 + yv + 1)^2 * (xv^2 + yv^2 - 1) + yv^2 - 5
f3 = xv^8 - 28 * xv^6 * yv^2 + 70 * xv^4 * yv^4 - 28 * xv^2 * yv^6 + yv^8 + 15 * xv^4 * yv^2 - 15 * xv^2 * yv^4

F1 = System([f1], variables = [xv, yv])
F2 = System([f2], variables = [xv, yv])
F3 = System([f3], variables = [xv, yv])

cfg_default = HomotopyConfig{Float64}()
cfg_wide = HomotopyConfig{Float64}(bbox_x = (-8.0, 8.0), bbox_y = (-8.0, 8.0))

function _run(label, F, cfg)
    println()
    println("=" ^ 70)
    println(label)
    println("=" ^ 70)

    t0 = @elapsed (vertices, edges) = decompose_1d_curve(F, cfg; deflate = true)
    sing = filter(v -> v.v_type == Singular, vertices)
    println("elapsed: $(round(t0; digits=3))s   total vertices: $(length(vertices))   Singular: $(length(sing))")

    red_flags = String[]
    for v in sing
        coords = round.(real.(v.coordinates); digits=6)
        if haskey(v.metadata, :isosingular_verdict)
            verdict = v.metadata[:isosingular_verdict]
            dim = v.metadata[:isosingular_dimension]
            seq = v.metadata[:isosingular_deflation_sequence]
            rounds = length(seq) - 1
            println("  id=$(v.id)  coords=$coords  verdict=$verdict  dim=$dim  rounds=$rounds  corank_seq=$seq")
            if verdict != Resolved
                push!(red_flags, "$label id=$(v.id) verdict=$verdict (expected Resolved)")
            elseif !(rounds in 1:2)
                push!(red_flags, "$label id=$(v.id) rounds=$rounds outside established [1,2] (soft flag)")
            end
        else
            println("  id=$(v.id)  coords=$coords  -- gated out (not applicable), no isosingular metadata")
        end
    end
    return vertices, edges, red_flags
end

all_red_flags = String[]

_, _, rf1 = _run("Curve 1: (x^2+y^2-1)^3 + 27x^2y^2 = 0  (default bbox)", F1, cfg_default)
append!(all_red_flags, rf1)

_, _, rf2 = _run("Curve 2: (x^3-xy^2+y+1)^2 (x^2+y^2-1) + y^2 - 5 = 0  (default bbox)", F2, cfg_default)
append!(all_red_flags, rf2)

_, _, rf3 = _run("Curve 3: x^8-28x^6y^2+70x^4y^4-28x^2y^6+y^8+15x^4y^2-15x^2y^4 = 0  (default bbox)", F3, cfg_default)
append!(all_red_flags, rf3)

_, _, rf3w = _run("Curve 3 (bbox=(-8,8), sensitivity comparison)", F3, cfg_wide)
append!(all_red_flags, rf3w)

println()
println("=" ^ 70)
println("Red-flag summary across all historical-curve runs")
println("=" ^ 70)
if isempty(all_red_flags)
    println("  none -- every Singular vertex either gated out cleanly or Resolved within [1,2] rounds")
else
    for r in all_red_flags
        println("  FLAG: ", r)
    end
end
