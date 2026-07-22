# dev/scratch/scratch_continuity_ok_diagnosis.jl
#
# Diagnostic investigation (2026-07, report-only per user instruction): does
# SurfaceIncidence.continuity_ok's flagged discontinuities (measured
# previously: 10/14 faces fixed-axis Taubin, 15/22 rotated seed 1) represent
# REAL branch-jumping (a sweep landing on the wrong topological branch -- a
# genuine correctness bug) or resolution-limited ambiguity (columns close to
# a genuine branch junction at current sampling densities, which
# _check_continuity!'s own docstring already says can't be distinguished
# below ~2x local chord spacing)?
#
# NOT a fix -- purely measurement. Reuses production functions unchanged.

using HomotopyContinuation
using LinearAlgebra
using GeometryBasics
using HomotopyGetsReal
using Random
using HomotopyGetsReal: _residual_at, _median_spacing

@var x y z
f = (x^2 + (1.2 * y)^2 + z^2 - 1)^3 - x^2 * z^3 - 0.1 * (1.2 * y)^2 * z^3
F_heart = System([f], variables = [x, y, z])
patch = HomotopyGetsReal.build_patch_system(F_heart)

"""
    _diagnose(cfg, projection, label) -> (inc, faces, records)

Runs decompose_3d_surface(...; incidence = true) and, for every flagged
face, collects one diagnostic record PER violating consecutive-column pair:
the two landing points (world coordinates), their direct 3D gap, their
surface residuals, the local chord spacing elsewhere on the SAME boundary
row (the scale the 2x-chord-spacing resolution limit is stated relative
to), and the distance from each point to the nearest genuine crit-slice
vertex / fold-anchor junction.
"""
function _diagnose(cfg::HomotopyConfig{Float64}, projection, label::String)
    kwargs = projection === nothing ? NamedTuple() : (; projection = projection)
    vi, ei, faces, mi, inc = decompose_3d_surface(F_heart, cfg; incidence = true, kwargs...)
    n_z = 2 * cfg.midslice_sample_density + 1
    id_to_face = Dict(fc.id => fc for fc in faces)

    edge_id_to_cs = Dict(e.id => cs for cs in inc.crit_slices for e in cs.edges)
    vertex_id_to_cs = Dict(v.id => cs for cs in inc.crit_slices for v in cs.vertices)
    junction_points = Vector{Float64}[]
    for cs in inc.crit_slices, v in cs.vertices
        push!(junction_points, Float64[real(v.coordinates[1]), real(v.coordinates[2]), cs.z])
    end
    for v in inc.critical_vertices
        push!(junction_points, Float64[real(v.coordinates[1]), real(v.coordinates[2]), real(v.coordinates[3])])
    end

    flagged = sort([fid for (fid, ok) in inc.continuity_ok if !ok])
    println("[$label] $(length(flagged)) / $(length(inc.continuity_ok)) faces flagged: $(flagged)")

    records = NamedTuple[]
    for fid in flagged
        face = id_to_face[fid]
        n_curve = size(face.mesh_vertices, 1) ÷ n_z
        for (side, viol_dict, row_off, landings_dict) in (
            (:bottom, inc.continuity_violations_bottom, 0, inc.column_landings_bottom),
            (:top, inc.continuity_violations_top, n_z - 1, inc.column_landings_top),
        )
            haskey(landings_dict, fid) || continue
            landings = landings_dict[fid]
            viol_cols = get(viol_dict, fid, Int[])
            isempty(viol_cols) && continue

            row_pts = Vector{Float64}[Vector{Float64}(face.mesh_vertices[row_off*n_curve+cc, :]) for cc in 1:n_curve]
            local_scale = _median_spacing(row_pts)

            for c in viol_cols
                p_c, p_c1 = row_pts[c], row_pts[c+1]
                d3 = norm(p_c .- p_c1)
                r_c = abs(_residual_at(patch, p_c[1], p_c[2], p_c[3], cfg))
                r_c1 = abs(_residual_at(patch, p_c1[1], p_c1[2], p_c1[3], cfg))

                l_c, l_c1 = landings[c], landings[c+1]
                cs_c = l_c.kind === :edge ? get(edge_id_to_cs, l_c.id, nothing) :
                       l_c.kind === :crit_slice_vertex ? get(vertex_id_to_cs, l_c.id, nothing) : nothing
                junc_dist = isempty(junction_points) ? NaN : minimum(norm(p_c .- jp) for jp in junction_points)

                push!(records, (
                    face = fid, side = side, c = c,
                    p_c = p_c, p_c1 = p_c1, d3 = d3,
                    r_c = r_c, r_c1 = r_c1,
                    local_scale = local_scale, ratio = local_scale > 0 ? d3 / local_scale : Inf,
                    kind_c = l_c.kind, id_c = l_c.id, dist_c = l_c.dist,
                    kind_c1 = l_c1.kind, id_c1 = l_c1.id, dist_c1 = l_c1.dist,
                    junc_dist = junc_dist,
                ))
            end
        end
    end
    return inc, faces, records
end

println("=" ^ 70)
println("PART 1: baseline measurement (matches test_taubin.jl's own cfg/seed)")
println("=" ^ 70)

cfg = HomotopyConfig{Float64}(
    bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3),
    edge_sample_density = 8, midslice_sample_density = 8,
)

inc_fx, faces_fx, rec_fx = _diagnose(cfg, nothing, "fixed-axis")
println()
Q1 = random_orthogonal_matrix(Float64, 3; rng = Xoshiro(1))
inc_rot, faces_rot, rec_rot = _diagnose(cfg, Q1, "rotated seed 1")

println()
println("=" ^ 70)
println("PART 2: per-violation evidence table (fixed-axis)")
println("=" ^ 70)
for r in rec_fx
    println("  face=$(r.face) side=$(r.side) c=$(r.c): d3=$(round(r.d3;sigdigits=4)) local_scale=$(round(r.local_scale;sigdigits=4)) ratio(d3/scale)=$(round(r.ratio;sigdigits=3))")
    println("    p_c =$(round.(r.p_c;digits=4))  resid=$(r.r_c)  landing=($(r.kind_c),$(r.id_c),dist=$(round(r.dist_c;sigdigits=3)))")
    println("    p_c1=$(round.(r.p_c1;digits=4))  resid=$(r.r_c1)  landing=($(r.kind_c1),$(r.id_c1),dist=$(round(r.dist_c1;sigdigits=3)))")
    println("    dist to nearest junction: p_c=$(round(r.junc_dist;sigdigits=4))")
end

println()
println("=" ^ 70)
println("PART 2b: per-violation evidence table (rotated seed 1)")
println("=" ^ 70)
for r in rec_rot
    println("  face=$(r.face) side=$(r.side) c=$(r.c): d3=$(round(r.d3;sigdigits=4)) local_scale=$(round(r.local_scale;sigdigits=4)) ratio(d3/scale)=$(round(r.ratio;sigdigits=3))")
    println("    p_c =$(round.(r.p_c;digits=4))  resid=$(r.r_c)  landing=($(r.kind_c),$(r.id_c),dist=$(round(r.dist_c;sigdigits=3)))")
    println("    p_c1=$(round.(r.p_c1;digits=4))  resid=$(r.r_c1)  landing=($(r.kind_c1),$(r.id_c1),dist=$(round(r.dist_c1;sigdigits=3)))")
    println("    dist to nearest junction: p_c=$(round(r.junc_dist;sigdigits=4))")
end

println()
println("=" ^ 70)
println("PART 3: Taubin-heart-specific structural ground truth (y -> -y mirror")
println("symmetry: f depends on y only via (1.2y)^2). A genuine cross-lobe jump")
println("should show p_c/p_c1 straddling y=0 with BOTH |y| well above the local")
println("noise floor -- not just a small perturbation near y=0 itself.")
println("=" ^ 70)
for (label, recs) in (("fixed-axis", rec_fx), ("rotated seed 1", rec_rot))
    println("  [$label]")
    for r in recs
        y_c, y_c1 = r.p_c[2], r.p_c1[2]
        straddles = sign(y_c) != sign(y_c1) && min(abs(y_c), abs(y_c1)) > 0.05
        println("    face=$(r.face) $(r.side) c=$(r.c): y_c=$(round(y_c;digits=4)) y_c1=$(round(y_c1;digits=4)) ",
                straddles ? "-> STRADDLES y=0 with both |y|>0.05 (candidate real cross-lobe jump)" :
                            "-> does not straddle meaningfully (consistent with same-lobe ambiguity)")
    end
end

println()
println("=" ^ 70)
println("PART 3a: violation category breakdown by (kind_c, kind_c1) pair --")
println("does the SAME critical_point/crit_slice_vertex kind-mismatch pattern")
println("found in the fixed-axis case account for a comparable share of the")
println("rotated case's violations too?")
println("=" ^ 70)
function _kind_breakdown(recs)
    counts = Dict{Tuple{Symbol,Symbol},Int}()
    for r in recs
        key = r.kind_c <= r.kind_c1 ? (r.kind_c, r.kind_c1) : (r.kind_c1, r.kind_c)
        counts[key] = get(counts, key, 0) + 1
    end
    return counts
end
println("  [fixed-axis] $(length(rec_fx)) total violations: $(_kind_breakdown(rec_fx))")
println("  [rotated seed 1] $(length(rec_rot)) total violations: $(_kind_breakdown(rec_rot))")

println()
println("=" ^ 70)
println("PART 3b: direct check of the fold-tip dual-representation hypothesis --")
println("for every fixed-axis violation pairing a :critical_point against a")
println(":crit_slice_vertex, do the TWO CELLS THEMSELVES (not just the swept")
println("columns landing near them) refer to the same physical location? IDs are")
println("re-derived from THIS run's own data (HC.jl assigns ids in a run-dependent")
println("order -- a hardcoded id from a previous run's printout would be wrong).")
println("=" ^ 70)
cp_by_id = Dict(v.id => v for v in inc_fx.critical_vertices)
csv_by_id = Dict(v.id => (v, cs) for cs in inc_fx.crit_slices for v in cs.vertices)
seen_pairs = Set{Tuple{Int,Int}}()
for r in rec_fx
    kc, kc1 = (r.kind_c, r.id_c), (r.kind_c1, r.id_c1)
    pair = if kc[1] === :critical_point && kc1[1] === :crit_slice_vertex
        (kc[2], kc1[2])
    elseif kc[1] === :crit_slice_vertex && kc1[1] === :critical_point
        (kc1[2], kc[2])
    else
        nothing
    end
    pair === nothing && continue
    pair in seen_pairs && continue
    push!(seen_pairs, pair)
    cp_id, csv_id = pair
    cp = cp_by_id[cp_id]
    csv, cs_owner = csv_by_id[csv_id]
    cp_coords = real.(cp.coordinates)
    csv_coords = Float64[real(csv.coordinates[1]), real(csv.coordinates[2]), cs_owner.z]
    d = norm(cp_coords .- csv_coords)
    println("  critical_point($cp_id) = $(round.(cp_coords;digits=6))  vs  crit_slice_vertex($csv_id) = $(round.(csv_coords;digits=6))  distance=$(d)")
end
println("  (for comparison, local chord spacing at z=1.2367 was ~0.00066 -- if a distance")
println("   above is comparably tiny, that :critical_point/:crit_slice_vertex pair is the")
println("   SAME physical point under two different cell ids, and _cells_adjacent's")
println("   conservative kind-mismatch rule is flagging a data-model artifact there, not")
println("   a real spatial jump.)")

println()
println("=" ^ 70)
println("PART 4: resolution sensitivity -- re-run at higher edge_sample_density")
println("and midslice_sample_density, same face ids (topology-preserving knobs),")
println("check whether flagged-face count / violation-pair count shrinks.")
println("=" ^ 70)
for (esd, msd) in ((8, 8), (16, 16), (24, 24))
    cfg_d = HomotopyConfig{Float64}(
        bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3),
        edge_sample_density = esd, midslice_sample_density = msd,
    )
    inc_d, faces_d, rec_d = _diagnose(cfg_d, nothing, "fixed-axis esd=$esd msd=$msd")
    n_flagged = count(!, values(inc_d.continuity_ok))
    n_viol = length(rec_d)
    by_z = Dict{Float64,Int}()
    for r in rec_d
        zc = round(r.p_c[3]; digits = 2)
        by_z[zc] = get(by_z, zc, 0) + 1
    end
    dup_ratio_stats = [round(r.ratio; sigdigits = 3) for r in rec_d]
    println("  esd=$esd msd=$msd: $(n_flagged) faces flagged, $(n_viol) violation pairs ($(faces_d |> length) total faces)")
    println("    by z: $(by_z)")
    println("    ratios (d3/local_scale): $(sort(dup_ratio_stats))")
end

println()
println("=" ^ 70)
println("PART 5: resolution sensitivity, ROTATED seed 1 -- same protocol, tracking")
println("specifically the 3 candidate real cross-lobe jumps flagged in Part 3")
println("(face 10 top c=3, face 15 bottom c=2, face 22 bottom c=6 at esd=msd=8).")
println("=" ^ 70)
for (esd, msd) in ((8, 8), (16, 16), (24, 24))
    cfg_d = HomotopyConfig{Float64}(
        bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3),
        edge_sample_density = esd, midslice_sample_density = msd,
    )
    inc_d, faces_d, rec_d = _diagnose(cfg_d, Q1, "rotated esd=$esd msd=$msd")
    n_flagged = count(!, values(inc_d.continuity_ok))
    n_viol = length(rec_d)
    n_straddle = count(rec_d) do r
        y_c, y_c1 = r.p_c[2], r.p_c1[2]
        sign(y_c) != sign(y_c1) && min(abs(y_c), abs(y_c1)) > 0.05
    end
    println("  esd=$esd msd=$msd: $(n_flagged) faces flagged, $(n_viol) violation pairs, ",
            "$(n_straddle) straddle-y=0-with-both-|y|>0.05 candidates ($(length(faces_d)) total faces)")
    for r in rec_d
        y_c, y_c1 = r.p_c[2], r.p_c1[2]
        straddles = sign(y_c) != sign(y_c1) && min(abs(y_c), abs(y_c1)) > 0.05
        straddles && println("    face=$(r.face) $(r.side) c=$(r.c): y_c=$(round(y_c;digits=4)) y_c1=$(round(y_c1;digits=4)) ratio=$(round(r.ratio;sigdigits=3)) d3=$(round(r.d3;sigdigits=4))")
    end
end
