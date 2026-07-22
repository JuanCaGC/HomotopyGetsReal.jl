# dev/scratch/scratch_phase9c_optionB_loft_check.jl
#
# Direction-change investigation (2026-07): before committing to Option A's
# larger, riskier targeted-tracking redesign (Phase 10, Stage 3), test
# whether Option B's FULL loft/zippering approach -- described but
# deliberately deferred in SurfaceDecomposition.jl's weld_mesh docstring
# ("a genuine conforming-triangulation/boundary-zippering step... a
# materially different and larger algorithm than snap-and-split") -- already
# closes the 58 residual naked edges at the Taubin heart's multi-face
# edge-type boundaries (z=1.0 notch, z=1.0648 saddle pair), which Phase 9b's
# snap-based weld_mesh could not (it can only relabel EXISTING swept-column
# landings onto a SUBSET of a crit-slice edge's own sample points -- it can
# never fabricate coverage of the points no column happened to land near).
#
# This script is a STANDALONE MEASUREMENT, not wired into
# decompose_3d_surface/weld_mesh. It reuses production functions
# (_snap_boundary_points!, _landing_confidence, _inward_row_points,
# _median_spacing, cluster_points_indexed, _naked_mesh_edges, _residual_at)
# unchanged, and adds ONE new piece of geometry: a classic two-polyline
# "bridge" (ruled-surface / zipper) triangulation via dynamic programming,
# which touches EVERY point of BOTH input polylines by construction --
# unlike Phase 9b's snap, which only touches the subset of crit-slice
# samples some column's monotone assignment happened to land on.
#
# SCOPE LIMITATION, stated up front (not discovered after the fact): this
# prototype only replaces the boundary-row triangle strip for a face+side
# whose ENTIRE confident :edge run is handled; the single quad at each END
# of a run that borders a non-loft column (a different landing kind, or a
# different crit-slice edge) is left using the OLD (already-snapped)
# boundary vertex, not the new loft vertex, since reconciling that seam
# would require the same fan-splitting machinery weld_mesh's
# _split_t_junctions already has for the analogous problem -- reusing it
# here is straightforward future work but out of scope for a same-day
# ablation measurement. Any residual naked edges concentrated at RUN
# ENDPOINTS (as opposed to spread across a run's interior) are exactly this
# known artifact, not a flaw in the core hypothesis being tested.

using HomotopyContinuation
using LinearAlgebra
using GeometryBasics
using HomotopyGetsReal
using HomotopyGetsReal: _snap_boundary_points!, _landing_confidence, _inward_row_points,
    _median_spacing, _naked_mesh_edges, _residual_at, _project_to_slice, _edge_spacing,
    cluster_points_indexed

println("=" ^ 70)
println("Setup: identical Taubin heart fixture to test_taubin.jl")
println("=" ^ 70)

@var x y z
f = (x^2 + (1.2 * y)^2 + z^2 - 1)^3 - x^2 * z^3 - 0.1 * (1.2 * y)^2 * z^3
F_heart = System([f], variables = [x, y, z])
cfg = HomotopyConfig{Float64}(
    bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3),
    edge_sample_density = 8, midslice_sample_density = 8,
)
patch = HomotopyGetsReal.build_patch_system(F_heart)

t_inc = @elapsed (vi, ei, fi, mi, inc) = decompose_3d_surface(F_heart, cfg; incidence = true)
println("  [$(round(t_inc; digits=1)) s]  $(length(fi)) faces, $(length(inc.crit_slices)) crit-slices")

naked_baseline = _naked_mesh_edges(mi)
mesh_pts = [Vector{Float64}(GeometryBasics.coordinates(mi)[i]) for i in eachindex(GeometryBasics.coordinates(mi))]
count_near(zc) = count(((i, j),) -> abs((mesh_pts[i][3] + mesh_pts[j][3]) / 2 - zc) < 0.01, naked_baseline)
println("  Phase 9b baseline (production weld_mesh, incidence-aware): $(length(naked_baseline)) naked edges")
println("    z=-1: $(count_near(-1.0))  z=1.0: $(count_near(1.0))  z=1.0648: $(count_near(1.0648))  z=1.2367: $(count_near(1.2367))")

cs_by_z = Dict(round(cs.z; digits = 4) => cs for cs in inc.crit_slices)
target_boundaries = [cs_by_z[1.0], cs_by_z[1.0648]]
target_js = Set(cs.boundary_index for cs in target_boundaries)
println("  Targeting boundary_index(es) $(sort(collect(target_js))) for the loft experiment (z=1.0, z=1.0648).")

println()
println("=" ^ 70)
println("Building each target crit-slice edge's OWN full point sequence (Q),")
println("re-projected exactly like _snap_boundary_points!'s edge_targets")
println("=" ^ 70)

edge_targets = Dict{Int,Vector{Vector{Float64}}}()
edge_to_cs = Dict{Int,typeof(target_boundaries[1])}()
for cs in target_boundaries, e in cs.edges
    tg = Vector{Vector{Float64}}(undef, length(e.sampled_points))
    for (k, p) in enumerate(e.sampled_points)
        xr, yr = _project_to_slice(patch, p[1], p[2], cs.z, cfg)
        tg[k] = Float64[xr, yr, cs.z]
    end
    edge_targets[e.id] = tg
    edge_to_cs[e.id] = cs
end
println("  $(length(edge_targets)) crit-slice edges across the two target boundaries, ",
        "$(sum(length(v) for v in values(edge_targets))) total Q points.")

"""
    _bridge_triangulate(P, Q) -> Vector{NTuple{2,Symbol}}... (moves)

Classic minimum-total-edge-length triangulation of the ribbon between two
ordered polylines `P` (length n) and `Q` (length m), via O(n*m) DP -- the
same family of algorithm as `_monotone_snap_targets`'s DP, but for a FULL
zipper (every point of both P and Q becomes a mesh vertex) rather than a
many-to-one snap assignment. Returns the sequence of moves backtracked from
(n,m) to (1,1): `:P` means "advance the P index" (triangle uses P[i-1],P[i],Q[j]),
`:Q` means "advance the Q index" (triangle uses P[i],Q[j-1],Q[j]).
"""
function _bridge_moves(P::Vector{Vector{Float64}}, Q::Vector{Vector{Float64}})
    n, m = length(P), length(Q)
    cost = fill(Inf, n, m)
    move = fill(:none, n, m)
    cost[1, 1] = norm(P[1] .- Q[1])
    for i in 2:n
        cost[i, 1] = cost[i-1, 1] + norm(P[i] .- Q[1])
        move[i, 1] = :P
    end
    for j in 2:m
        cost[1, j] = cost[1, j-1] + norm(P[1] .- Q[j])
        move[1, j] = :Q
    end
    for i in 2:n, j in 2:m
        optP = cost[i-1, j] + norm(P[i] .- Q[j])
        optQ = cost[i, j-1] + norm(P[i] .- Q[j])
        if optP <= optQ
            cost[i, j] = optP
            move[i, j] = :P
        else
            cost[i, j] = optQ
            move[i, j] = :Q
        end
    end
    moves = Symbol[]
    i, j = n, m
    while (i, j) != (1, 1)
        mv = move[i, j]
        push!(moves, mv)
        mv === :P ? (i -= 1) : (j -= 1)
    end
    return reverse(moves), cost[n, m]
end

"""
    _loft_ribbon(P, Q) -> Vector{NTuple{3,Vector{Float64}}}

Triangles (as literal point triples, not indices -- caller assigns global
vertex ids after clustering) for the ribbon between P and Q, built by
walking `_bridge_moves`' backtrace. Every P[i] and every Q[j] appears in at
least one returned triangle (full coverage by construction).
"""
function _loft_ribbon(P::Vector{Vector{Float64}}, Q::Vector{Vector{Float64}})
    moves, cost = _bridge_moves(P, Q)
    tris = NTuple{3,Vector{Float64}}[]
    i, j = 1, 1
    for mv in moves
        if mv === :P
            push!(tris, (P[i], P[i+1], Q[j]))
            i += 1
        else
            push!(tris, (P[i], Q[j], Q[j+1]))
            j += 1
        end
    end
    return tris, cost
end

println()
println("=" ^ 70)
println("Sanity check on _bridge_moves/_loft_ribbon: equal-length straight lines")
println("=" ^ 70)
P_test = [Float64[i, 0, 0] for i in 0:4]
Q_test = [Float64[i, 1, 0] for i in 0:4]
tris_test, cost_test = _loft_ribbon(P_test, Q_test)
verts_in_tris = Set{Vector{Float64}}()
for t in tris_test, v in t
    push!(verts_in_tris, v)
end
@assert length(tris_test) == length(P_test) + length(Q_test) - 2
@assert all(p -> p in verts_in_tris, P_test)
@assert all(q -> q in verts_in_tris, Q_test)
println("  $(length(tris_test)) triangles (expected $(length(P_test)+length(Q_test)-2)), ",
        "every P and Q point used at least once. OK.")

println()
println("=" ^ 70)
println("Identifying confident, target-boundary :edge column runs per face/side")
println("(mirrors _snap_boundary_points!'s own pass-1 grouping logic)")
println("=" ^ 70)

id_to_fi = Dict(f.id => i for (i, f) in enumerate(fi))
edge_spacing = Dict{Int,Float64}(e.id => _edge_spacing(e) for cs in inc.crit_slices for e in cs.edges)

struct LoftRun
    fi::Int
    side::Symbol       # :bottom or :top
    cols::UnitRange{Int}   # inclusive column range in the face's curve order
    edge_id::Int
end

runs = LoftRun[]
n_skipped_low_confidence = 0
for (fidx, face) in enumerate(fi)
    n_pts = size(face.mesh_vertices, 1)
    n_pts == 0 && continue
    n_z = 2 * cfg.midslice_sample_density + 1
    n_curve = n_pts ÷ n_z
    for (landings_dict, side) in ((inc.column_landings_bottom, :bottom), (inc.column_landings_top, :top))
        haskey(landings_dict, face.id) || continue
        landings = landings_dict[face.id]
        length(landings) == n_curve || continue
        row_pts = _inward_row_points(face, side, n_z, n_curve)
        scale_ref = _median_spacing(row_pts)
        boundary_row = side === :bottom ? [Vector{Float64}(face.mesh_vertices[c, :]) for c in 1:n_curve] :
                                           [Vector{Float64}(face.mesh_vertices[(n_z-1)*n_curve+c, :]) for c in 1:n_curve]

        c = 1
        while c <= n_curve
            l = landings[c]
            is_target_edge = l.kind === :edge && haskey(edge_to_cs, l.id)
            if is_target_edge && _landing_confidence(boundary_row[c], l, edge_spacing, scale_ref, patch, cfg)
                c0 = c
                eid = l.id
                while c <= n_curve && landings[c].kind === :edge && landings[c].id == eid &&
                      _landing_confidence(boundary_row[c], landings[c], edge_spacing, scale_ref, patch, cfg)
                    c += 1
                end
                push!(runs, LoftRun(fidx, side, c0:(c-1), eid))
            else
                is_target_edge && (global n_skipped_low_confidence += 1)
                c += 1
            end
        end
    end
end
println("  $(length(runs)) confident target-boundary :edge runs found across all faces/sides.")
println("  $(n_skipped_low_confidence) individual columns skipped (target edge, but NOT confident).")
for r in runs
    println("    face $(fi[r.fi].id) $(r.side) cols $(r.cols) -> crit-slice edge $(r.edge_id)")
end

println()
println("=" ^ 70)
println("Rebuilding the mesh: production triangulation everywhere, EXCEPT the")
println("interior of each identified run's boundary-row quads, replaced by the")
println("loft ribbon against that edge's FULL Q sequence")
println("=" ^ 70)

all_points = Vector{Float64}[]
provenance = Tuple{Int,Int}[]
index_of = Dict{Tuple{Int,Int},Int}()
for (fidx, face) in enumerate(fi)
    for r in 1:size(face.mesh_vertices, 1)
        push!(all_points, face.mesh_vertices[r, :])
        push!(provenance, (fidx, r))
        index_of[(fidx, r)] = length(all_points)
    end
end

# Same snap-unification pass production weld_mesh runs (still closes the
# fold boundaries and gives partial credit at z=1.0/1.0648 exactly as before;
# irrelevant for columns whose triangles we're about to replace entirely).
_snap_boundary_points!(all_points, index_of, fi, inc, patch, cfg)

# Dropped-quad lookup: (fidx, side, c) -> true if quad column c (i.e. the
# quad spanning curve-columns c,c+1 at the outer row) is FULLY inside some
# run's interior (both c and c+1 in the SAME run) and must be dropped from
# the ordinary triangulation.
dropped = Dict{Tuple{Int,Symbol,Int},Bool}()
for r in runs, c in first(r.cols):(last(r.cols) - 1)
    dropped[(r.fi, r.side, c)] = true
end

global_triangles = NTuple{3,Int}[]
loft_new_points = Vector{Float64}[]
loft_provenance_tag = Tuple{Symbol,Int,Int}[]   # (:loft, run_index, local_point_index) for bookkeeping only

for (fidx, face) in enumerate(fi)
    n_pts = size(face.mesh_vertices, 1)
    n_pts == 0 && continue
    n_z = 2 * cfg.midslice_sample_density + 1
    n_curve = n_pts ÷ n_z
    for r_local in 1:(n_z - 1), c in 1:(n_curve - 1)
        is_outer_bottom = r_local == 1
        is_outer_top = r_local == n_z - 1
        skip = (is_outer_bottom && get(dropped, (fidx, :bottom, c), false)) ||
               (is_outer_top && get(dropped, (fidx, :top, c), false))
        skip && continue
        v1 = (r_local - 1) * n_curve + c
        v2 = r_local * n_curve + c
        v3 = r_local * n_curve + (c + 1)
        v4 = (r_local - 1) * n_curve + (c + 1)
        g1, g2, g3, g4 = index_of[(fidx, v1)], index_of[(fidx, v2)], index_of[(fidx, v3)], index_of[(fidx, v4)]
        g1 != g2 && g2 != g3 && g3 != g1 && push!(global_triangles, (g1, g2, g3))
        g1 != g3 && g3 != g4 && g4 != g1 && push!(global_triangles, (g1, g3, g4))
    end
end

total_loft_cost = 0.0
q_raw_idxs = Int[]   # diagnostic: raw all_points indices of every pushed Q-derived vertex
# (b, part 1 of 2) per-edge-id, ordered flat index of EACH Q sample k (1..m),
# FIRST occurrence wins -- exactly _snap_boundary_points!'s own
# target_flat_idx bookkeeping pattern, just built from the FULL Q sequence
# (every k, since the loft touches all of them) instead of only "used" k's.
q_target_flat_idx = Dict{Tuple{Int,Int},Int}()
# (a) seam-cap bookkeeping: for each run, the flat index of its P[1]/P[end]
# and the Q vertex adjacent to them in the CHOSEN orientation, so the cap
# triangles (below) can reference the ribbon's own real endpoints.
run_end_info = Vector{NamedTuple}(undef, length(runs))

for (ridx, r) in enumerate(runs)
    face = fi[r.fi]
    n_pts = size(face.mesh_vertices, 1)
    n_z = 2 * cfg.midslice_sample_density + 1
    n_curve = n_pts ÷ n_z
    # P = the row ONE STEP IN from the boundary (row 2 for :bottom, row n_z-1
    # for :top) -- NOT the boundary row itself, which is exactly the row being
    # discarded/replaced. This is the same row _snap_boundary_points!'s own
    # scale_ref uses (_inward_row_points), so reusing it directly rather than
    # re-deriving the row-index arithmetic by hand.
    P_full = _inward_row_points(face, r.side, n_z, n_curve)
    P = P_full[r.cols]
    Q = edge_targets[r.edge_id]
    Q_rev = reverse(Q)
    tris_fwd, cost_fwd = _loft_ribbon(P, Q)
    tris_rev, cost_rev = _loft_ribbon(P, Q_rev)
    forward = cost_fwd <= cost_rev
    tris, cost = forward ? (tris_fwd, cost_fwd) : (tris_rev, cost_rev)
    Q_used = forward ? Q : Q_rev
    m = length(Q)
    orig_k(k_used) = forward ? k_used : (m + 1 - k_used)
    global total_loft_cost += cost

    p1_flat = 0
    pend_flat = 0
    q_first_flat = 0
    q_last_flat = 0
    for (pa, pb, pc) in tris
        idxs = Int[]
        for p in (pa, pb, pc)
            push!(all_points, p)
            push!(idxs, length(all_points))
            if p in Q
                push!(q_raw_idxs, length(all_points))
                k_used = findfirst(qp -> qp == p, Q_used)
                k = orig_k(k_used)
                haskey(q_target_flat_idx, (r.edge_id, k)) || (q_target_flat_idx[(r.edge_id, k)] = length(all_points))
                k_used == 1 && q_first_flat == 0 && (q_first_flat = length(all_points))
                k_used == m && q_last_flat == 0 && (q_last_flat = length(all_points))
            elseif p === P[1] || p == P[1]
                p1_flat == 0 && (p1_flat = length(all_points))
            elseif p === P[end] || p == P[end]
                pend_flat == 0 && (pend_flat = length(all_points))
            end
        end
        idxs[1] != idxs[2] && idxs[2] != idxs[3] && idxs[3] != idxs[1] && push!(global_triangles, (idxs[1], idxs[2], idxs[3]))
    end
    run_end_info[ridx] = (p1_flat = p1_flat, pend_flat = pend_flat, q_first_flat = q_first_flat, q_last_flat = q_last_flat)
end
println("  $(length(runs)) runs lofted, total DP bridge cost (sum of chosen-orientation edge lengths) = $(round(total_loft_cost; digits=4))")

# (a) part 2: seam-capping triangles. A run's LEFT end needs a cap iff its
# first column isn't the face's own first column (c0>1, so a KEPT transition
# quad still borders it using the OLD boundary-row vertex); symmetrically for
# the RIGHT end. The cap triangle (old boundary vertex, ribbon's own P[1] or
# P[end], ribbon's own adjacent Q vertex) reproduces the SPECIFIC edge
# (old-boundary-vertex -- P[1]) that the dropped quad used to contribute a
# second triangle to -- exactly the edge the "seam" naked-edge diagnosis
# identified as under-covered (1 triangle instead of 2) once its neighboring
# quad was dropped.
n_caps = 0
for (ridx, r) in enumerate(runs)
    face = fi[r.fi]
    n_pts = size(face.mesh_vertices, 1)
    n_z = 2 * cfg.midslice_sample_density + 1
    n_curve = n_pts ÷ n_z
    row_offset = r.side === :bottom ? 0 : (n_z - 1)
    info = run_end_info[ridx]
    c0, c1 = first(r.cols), last(r.cols)
    if c0 > 1 && info.p1_flat != 0 && info.q_first_flat != 0
        old_v = index_of[(r.fi, row_offset * n_curve + c0)]
        push!(global_triangles, (old_v, info.p1_flat, info.q_first_flat))
        global n_caps += 1
    end
    if c1 < n_curve && info.pend_flat != 0 && info.q_last_flat != 0
        old_v = index_of[(r.fi, row_offset * n_curve + c1)]
        push!(global_triangles, (old_v, info.pend_flat, info.q_last_flat))
        global n_caps += 1
    end
end
println("  $(n_caps) seam-cap triangles added.")
println("  $(length(global_triangles)) triangles total (production + loft + caps, pre-clustering).")

reps, membership = cluster_points_indexed(all_points, cfg.vertex_match_tol)
fixed_triangles = NTuple{3,Int}[]
for (g1, g2, g3) in global_triangles
    a, b, c = membership[g1], membership[g2], membership[g3]
    a != b && b != c && c != a && push!(fixed_triangles, (a, b, c))
end

# (b, part 2 of 2): reuse _split_t_junctions, extended to chain ADJACENT
# crit-slice edges (sharing a vertex) into ONE combined polyline, so a
# junction point shared between edge_id A and edge_id B is no longer
# outside _split_t_junctions' single-polyline scope (confirmed by direct
# check: 10/16 residual Q-Q naked edges had endpoints on DIFFERENT
# edge_ids -- exactly this junction case). Chains are built once per
# target CritSlice via a simple greedy walk over the edge-adjacency graph
# (vertex id -> touching edges); each chain's combined sequence is the
# concatenation of its edges' own Q sequences in walk order (reversing an
# edge's own sequence when the walk traverses it right-to-left),
# deduplicating the shared vertex at each join.
function _edge_seq(eid, edge_targets, q_target_flat_idx, membership)
    m = length(edge_targets[eid])
    return Int[membership[q_target_flat_idx[(eid, k)]] for k in 1:m if haskey(q_target_flat_idx, (eid, k))]
end

function _chain_polylines(cs, edge_targets, q_target_flat_idx, membership)
    left_of = Dict(e.id => e.left_vertex_id for e in cs.edges)
    right_of = Dict(e.id => e.right_vertex_id for e in cs.edges)
    touching = Dict{Int,Vector{Int}}()
    for e in cs.edges
        push!(get!(touching, e.left_vertex_id, Int[]), e.id)
        push!(get!(touching, e.right_vertex_id, Int[]), e.id)
    end
    visited = Set{Int}()
    out = Dict{Int,Vector{Int}}()
    for e0 in cs.edges
        e0.id in visited && continue
        chain = Tuple{Int,Bool}[(e0.id, true)]
        push!(visited, e0.id)
        cur_v = right_of[e0.id]
        while true
            cands = [eid for eid in get(touching, cur_v, Int[]) if !(eid in visited)]
            isempty(cands) && break
            nxt = first(cands)
            push!(visited, nxt)
            fwd = left_of[nxt] == cur_v
            push!(chain, (nxt, fwd))
            cur_v = fwd ? right_of[nxt] : left_of[nxt]
        end
        cur_v = left_of[e0.id]
        while true
            cands = [eid for eid in get(touching, cur_v, Int[]) if !(eid in visited)]
            isempty(cands) && break
            prv = first(cands)
            push!(visited, prv)
            fwd = right_of[prv] == cur_v
            pushfirst!(chain, (prv, fwd))
            cur_v = fwd ? left_of[prv] : right_of[prv]
        end
        seq = Int[]
        for (eid, fwd) in chain
            s = _edge_seq(eid, edge_targets, q_target_flat_idx, membership)
            s2 = fwd ? s : reverse(s)
            for r in s2
                (isempty(seq) || seq[end] != r) && push!(seq, r)
            end
        end
        isempty(seq) || (out[e0.id] = seq)
    end
    return out
end

edge_polylines = Dict{Int,Vector{Int}}()
for cs in target_boundaries
    merge!(edge_polylines, _chain_polylines(cs, edge_targets, q_target_flat_idx, membership))
end
n_before_split = length(fixed_triangles)
total_polyline_coverage = sum(length(v) for v in values(edge_polylines))
println("  chained into $(length(edge_polylines)) connected polylines (from $(length(edge_targets)) individual edges).")
fixed_triangles = HomotopyGetsReal._split_t_junctions(fixed_triangles, edge_polylines)
println("  _split_t_junctions: $(length(edge_polylines)) edge polylines, ",
        "$(total_polyline_coverage) total covered Q-points (of $(sum(length(v) for v in values(edge_targets))) possible), ",
        "$(n_before_split) -> $(length(fixed_triangles)) triangles.")

points3 = [GeometryBasics.Point3f(Float32(p[1]), Float32(p[2]), Float32(p[3])) for p in reps]
tris3 = [GeometryBasics.TriangleFace{Int}(t[1], t[2], t[3]) for t in fixed_triangles]
mesh_loft = GeometryBasics.Mesh(points3, tris3)

naked_loft = _naked_mesh_edges(mesh_loft)
mesh_pts_loft = [Vector{Float64}(GeometryBasics.coordinates(mesh_loft)[i]) for i in eachindex(GeometryBasics.coordinates(mesh_loft))]
count_near_loft(zc) = count(((i, j),) -> abs((mesh_pts_loft[i][3] + mesh_pts_loft[j][3]) / 2 - zc) < 0.01, naked_loft)

println()
println("=" ^ 70)
println("RESULT: naked-edge count, production (9b snap) vs Option B (full loft)")
println("=" ^ 70)
println("  Phase 9b (production):  total $(length(naked_baseline))  |  z=-1: $(count_near(-1.0))  z=1.0: $(count_near(1.0))  z=1.0648: $(count_near(1.0648))  z=1.2367: $(count_near(1.2367))")
println("  Option B (full loft):   total $(length(naked_loft))  |  z=-1: $(count_near_loft(-1.0))  z=1.0: $(count_near_loft(1.0))  z=1.0648: $(count_near_loft(1.0648))  z=1.2367: $(count_near_loft(1.2367))")

println()
println("=" ^ 70)
println("Diagnostic: classify each residual loft naked edge by endpoint type")
println("(Q = a crit-slice edge's own point; other = face-interior/seam point)")
println("=" ^ 70)
q_global = Set(membership[i] for i in q_raw_idxs)
n_qq = count(((i, j),) -> (i in q_global) && (j in q_global), naked_loft)
n_seam = count(((i, j),) -> xor(i in q_global, j in q_global), naked_loft)
n_other = length(naked_loft) - n_qq - n_seam
println("  Q-Q edges (both endpoints are crit-slice sample points): $(n_qq)")
println("  seam edges (one endpoint Q, one endpoint P/other):        $(n_seam)")
println("  other edges (neither endpoint a loft Q point):            $(n_other)")

println()
println("=" ^ 70)
println("Cross-edge-junction hypothesis check: for the Q-Q naked edges, do the")
println("two endpoints belong to the SAME crit-slice edge_id (a genuine within-")
println("edge skip _split_t_junctions should have caught) or DIFFERENT edge_ids")
println("(a junction point shared between adjacent crit-slice edges -- outside")
println("_split_t_junctions' single-polyline scope, since it requires eid_p==eid_q)?")
println("=" ^ 70)
rep_to_edge_ids = Dict{Int,Set{Int}}()
for ((eid, _), flat_idx) in q_target_flat_idx
    g = membership[flat_idx]
    push!(get!(rep_to_edge_ids, g, Set{Int}()), eid)
end
qq_edges = [(i, j) for (i, j) in naked_loft if (i in q_global) && (j in q_global)]
n_same_edge = 0
n_diff_edge = 0
for (i, j) in qq_edges
    eids_i = get(rep_to_edge_ids, i, Set{Int}())
    eids_j = get(rep_to_edge_ids, j, Set{Int}())
    shares = !isempty(intersect(eids_i, eids_j))
    shares ? (global n_same_edge += 1) : (global n_diff_edge += 1)
    println("    ($i,$j): endpoint edge_ids $(sort(collect(eids_i))) vs $(sort(collect(eids_j))) -> ",
            shares ? "SAME edge (should have been fixed)" : "DIFFERENT edges (junction point, out of scope)")
end
println("  -> $(n_same_edge) same-edge (unexpected, _split_t_junctions should have caught), ",
        "$(n_diff_edge) different-edge (junction points, confirms/refutes hypothesis)")

println()
println("=" ^ 70)
println("Precision: residuals of the crit-slice Q points actually used by the loft")
println("(these are the SAME already-Newton-polished points Phase 9a produced --")
println(" zero additional numerical solving performed on them here)")
println("=" ^ 70)
q_residuals = Float64[]
for r in runs
    for p in edge_targets[r.edge_id]
        push!(q_residuals, abs(_residual_at(patch, p[1], p[2], p[3], cfg)))
    end
end
unique!(q_residuals)
sort!(q_residuals)
println("  n=$(length(q_residuals))  min=$(q_residuals[1])  median=$(q_residuals[cld(end,2)])  ",
        "p90=$(q_residuals[ceil(Int, 0.9*length(q_residuals))])  max=$(q_residuals[end])")
println("  (Option A's targeted-hop landings would ADD tracker predictor-corrector")
println("   error on top of converging TO these same coordinates -- see report for")
println("   the full argument; not run here per the user's own permission to reason")
println("   from existing evidence rather than re-implementing Option A for this check.)")
