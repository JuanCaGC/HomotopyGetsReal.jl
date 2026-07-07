# src/Topology.jl
#
# Phase 3: steps 3-6 of the six-step algorithmic framework, specialized
# to plane curves (2-variable varieties) -- the "1D Curve Decomposition
# Pipeline":
#   3. Interslice / MidSlice!      -> compute_midslice
#   4. ConnectTheDots!             -> connect_the_dots!
#   5. Merge / GetMergeCandidates  -> decompose_1d_curve's cross-source
#                                      cluster_vertices call, plus the new
#                                      Clustering.cluster_scalars primitive
#   6. sample_edge
#
# This is the "MidSlice-First Strategy": path tracking always starts
# from a smooth, provably-non-singular interior point (a
# `compute_midslice` witness), never from a critical/singular vertex
# itself. The 3D analogue (z-slicing between critical Z-values to build
# a midslice MESH, `HomotopyConfig.midslice_sample_density`) is a
# different, later thing (Phase 5) that will call `decompose_1d_curve`
# once per z-slice -- it is NOT used here.
#
# Precision boundary (same reasoning as Phase 2's Solver.jl):
#   - `solve` (inside compute_midslice) and the `PathTracking.jl` engine
#     (`build_tracker`/`track_bidirectional`, used inside
#     connect_the_dots!) are genuinely Float64/ComplexF64-only. All path
#     tracking below happens in Float64 regardless of `T`.
#   - When `T != Float64`, only the FINAL landing point at each end of a
#     tracked branch is Newton-polished (via `Solver._newton_polish`,
#     reused unchanged) to genuine `T`-precision, immediately before it
#     is matched against `vertices` or turned into a new `NativeVertex`.
#     Raw intermediate path samples collected during Float64 tracking
#     are stored as `real.(...)` cast to `T` without a full re-polish.
#   - `compute_midslice`'s own `solve` call and `cluster_scalars`'s
#     scalar clustering are otherwise the only other numerically
#     "interesting" operations here; the former is Float64-tracked (its
#     result is a bare `Vector{Complex{T}}` witness that is *not* found
#     via `_newton_polish`, since it is deliberately a smooth interior
#     point, not a solution being classified/matched) and the latter is
#     already fully `T`-generic arithmetic (see Clustering.jl).
#
# Phase 4 note: the adaptive bidirectional path-tracking engine that
# `connect_the_dots!` uses (`_is_near_singular_f64`, `_track_segment!`,
# `_track_outward!` in earlier Phase 3 code) has been extracted,
# generalized, and hardened into `src/PathTracking.jl`
# (`is_near_singular`, `build_tracker`, `track_path`,
# `track_bidirectional`) so it is public and reusable by Phase 5's
# future 3D face-tracking, instead of being private to this file. The
# algorithm itself (and this file's behavior) is unchanged by that
# refactor -- see `PathTracking.jl`'s module docstring for the exact
# diff of what changed vs. what stayed the same.

"""
    compute_midslice(F::System, x_left::T, x_right::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> Vector{Complex{T}}

Step 3 of the six-step framework (plane-curve specialization -- see the
module docstring for why this is *not* the 3D z-slicing operation of
the same conceptual name in a later phase).

`F` must be a raw plane curve: `length(F.variables) == 2 &&
length(F.expressions) == 1` (same convention as
[`intersect_bounding_object`](@ref)'s raw-`F` input), else
`ArgumentError`.

Fixes `x_var` (the first of `F.variables`, matching
[`compute_critical_points`](@ref)'s existing positional destructuring
convention) at the midpoint `x_mid = (x_left + x_right) / 2`, path-tracks
the resulting univariate system in `y_var` (the second variable) in
Float64, and returns every root within `cfg.critical_point_tol` of real
(the tolerance is *reused*, not reintroduced as a new `HomotopyConfig`
field -- see the Phase 3 architecture discussion) as a
`Vector{Complex{T}}`.

No internal smoothness check is performed: this is guaranteed by
construction as long as callers only ever invoke this strictly between
two adjacent *distinct* x-values (see
[`Clustering.cluster_scalars`](@ref)), never at a known vertex's own x.
"""
function compute_midslice(F::System, x_left::T, x_right::T, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    length(F.variables) == 2 && length(F.expressions) == 1 || throw(ArgumentError(
        "compute_midslice: expected F with exactly 1 equation in exactly 2 variables " *
        "(a raw plane curve); got $(length(F.expressions)) equation(s) in $(length(F.variables)) variable(s).",
    ))
    x_var, y_var = F.variables
    f = F.expressions[1]

    x_mid = (x_left + x_right) / T(2)
    # Substitute a Float64 copy of x_mid before solving -- see the
    # analogous fix (and its rationale) in Solver.jl's
    # intersect_bounding_object: HC's polyhedral start-system
    # construction has no method for `Complex{BigFloat}` coefficients,
    # so a system built from a literal BigFloat substitution cannot be
    # passed to `solve` directly.
    Fsub = System([subs(f, x_var => Float64(x_mid))], variables = [y_var])

    result = solve(Fsub; show_progress = false)
    raw_sols = solutions(result; only_nonsingular = false)

    tol64 = Float64(cfg.critical_point_tol)
    roots = Complex{T}[]
    for s in raw_sols
        maximum(abs, imag.(s)) <= tol64 || continue
        push!(roots, Complex{T}(s[1]))
    end
    return roots
end

"""
    _resolve_endpoint(F, x_var, y_var, x_val::T, y_guess::Complex{T}, vertices::Vector{NativeVertex{T}}, cfg::HomotopyConfig{T}) where {T}
        -> (vertex_id::Int, y_final::Complex{T})

Shared endpoint-resolution logic used by [`connect_the_dots!`](@ref) at
both ends of a tracked branch:

1. Newton-polishes `y_guess` to genuine `T`-precision (a no-op when
   `T === Float64`) against the reduced *square* system obtained by
   fixing `x_var => x_val` in `F` -- exactly the same
   fix-one-variable-then-polish-the-rest pattern
   [`intersect_bounding_object`](@ref) already uses, reusing
   `Solver._newton_polish` unchanged.
2. Searches the full `vertices` list (never by sorted-list position --
   see the Phase 3 naive-adjacency counterexample) for the closest
   match by actual `(x, y)` Euclidean distance, using
   `cfg.vertex_match_tol` as the match radius.
3. If no vertex matches, mutates `vertices` in place: appends a new
   `Artificial`-typed `NativeVertex{T}` (id = current max id + 1) via
   the config-aware constructor from Phase 2, and returns its id.

Tags this fallback vertex's `metadata[:origin] = :endpoint_fallback`
(purely additive -- `NativeVertex`'s `metadata::Dict{Symbol,Any}` field
was designed in Phase 1 exactly for this kind of per-vertex provenance
note, no struct/constructor change needed). This distinguishes THIS
specific "a tracked path failed to close back onto any known vertex"
provenance from the OTHER, unrelated way an `Artificial` vertex can
arise in [`decompose_1d_curve`](@ref) -- `Clustering.cluster_vertices`
merging a `Critical` and a `Boundary` point that coincidentally land
within `vertex_match_tol` of each other (a benign near-tangency, not a
sign of anything wrong with the curve). Consumers that care about the
difference (currently: [`SurfaceDecomposition._robust_slice_at_z`](@ref),
which treats only THIS tag as a signal that a candidate `z_mid` produced
an untrustworthy slice) check `metadata[:origin]`, not just `v_type ==
Artificial`.
"""
function _resolve_endpoint(
    F::System,
    x_var,
    y_var,
    x_val::T,
    y_guess::Complex{T},
    vertices::Vector{NativeVertex{T}},
    cfg::HomotopyConfig{T},
) where {T<:AbstractFloat}
    f = F.expressions[1]
    Fsub = System([subs(f, x_var => x_val)], variables = [y_var])
    y_polished = _newton_polish(Fsub, Complex{T}[y_guess], cfg)[1]
    coord = Complex{T}[Complex{T}(x_val), y_polished]

    best_idx = 0
    best_dist = typemax(T)
    for (i, v) in enumerate(vertices)
        d = norm(v.coordinates .- coord)
        if d < best_dist
            best_dist = d
            best_idx = i
        end
    end

    if best_idx != 0 && best_dist <= cfg.vertex_match_tol
        return vertices[best_idx].id, y_polished
    else
        new_id = maximum(v.id for v in vertices) + 1
        push!(vertices, NativeVertex(cfg, new_id, coord, Artificial; metadata = Dict{Symbol,Any}(:origin => :endpoint_fallback)))
        return new_id, y_polished
    end
end

"""
    connect_the_dots!(
        F::System,
        x_left::T, x_mid::T, x_right::T,
        y_mid::Complex{T},
        edge_id::Int,
        vertices::Vector{NativeVertex{T}},
        cfg::HomotopyConfig{T},
    ) where {T<:AbstractFloat}
        -> Edge{T}

Step 4 of the six-step framework. Called once per branch, i.e. once per
`y`-root returned by [`compute_midslice`](@ref) at `x_mid = (x_left +
x_right) / 2`.

Tracks outward from the smooth midslice witness `(x_mid, y_mid)` toward
`x_left` and `x_right` independently via [`build_tracker`](@ref)/
[`track_bidirectional`](@ref) (Phase 4's `PathTracking.jl` engine,
parameterized by `x_var`), using `cfg.max_path_steps` as the raw
tracking step budget for *each* direction (adaptive bisection near
detected singularities/failures -- see
[`_track_path_segment!`](@ref) -- still bounded overall by this budget;
`cfg.edge_sample_density` is never used here, it is reserved
exclusively for [`sample_edge`](@ref)'s later equidistant resampling).

At each end, the landing point is resolved via [`_resolve_endpoint`](@ref)
against the *full* `vertices` list by actual `(x, y)` distance
(`cfg.vertex_match_tol`) -- never by position in a sorted list (the
Phase 3 naive-adjacency counterexample). `vertices` is mutated in place
(hence the `!`) whenever a tracked branch lands on no known vertex: a
new `Artificial`-typed `NativeVertex{T}` is appended so later calls in
the same `decompose_1d_curve` run see it too.

`edge_id` is supplied by the caller; `connect_the_dots!` never
self-assigns edge ids. `Edge.is_singular` is `true` iff either resolved
endpoint vertex (matched or newly-created) has `v_type == Singular`. The
returned `Edge{T}` has raw (non-equidistant) `sampled_points` already
populated, running left-endpoint -> ... -> midpoint -> ... ->
right-endpoint.
"""
function connect_the_dots!(
    F::System,
    x_left::T,
    x_mid::T,
    x_right::T,
    y_mid::Complex{T},
    edge_id::Int,
    vertices::Vector{NativeVertex{T}},
    cfg::HomotopyConfig{T},
) where {T<:AbstractFloat}
    length(F.variables) == 2 && length(F.expressions) == 1 || throw(ArgumentError(
        "connect_the_dots!: expected F with exactly 1 equation in exactly 2 variables " *
        "(a raw plane curve); got $(length(F.expressions)) equation(s) in $(length(F.variables)) variable(s).",
    ))
    x_var, y_var = F.variables
    f = F.expressions[1]

    H_sys = System([f], variables = [y_var], parameters = [x_var])
    xm64 = Float64(x_mid)
    ph, tracker = build_tracker(H_sys, xm64, cfg)
    y_mid64 = ComplexF64[ComplexF64(y_mid)]

    rank_tol64 = Float64(cfg.jacobian_rank_tol)
    sv_thresh64 = Float64(cfg.singular_value_threshold)
    poor_acc_tol64 = Float64(cfg.critical_point_tol)
    min_width64 = Float64(cfg.vertex_match_tol)

    full_path64, y_land_left, y_land_right = track_bidirectional(
        F, ph, tracker, y_mid64, xm64, Float64(x_left), Float64(x_right),
        cfg.max_path_steps, rank_tol64, sv_thresh64, poor_acc_tol64, min_width64,
    )

    left_id, _ = _resolve_endpoint(F, x_var, y_var, x_left, Complex{T}(y_land_left[1]), vertices, cfg)
    right_id, _ = _resolve_endpoint(F, x_var, y_var, x_right, Complex{T}(y_land_right[1]), vertices, cfg)

    left_vertex = vertices[findfirst(v -> v.id == left_id, vertices)]
    right_vertex = vertices[findfirst(v -> v.id == right_id, vertices)]

    full_path = Vector{T}[]
    push!(full_path, T[real(left_vertex.coordinates[1]), real(left_vertex.coordinates[2])])
    for p in full_path64
        push!(full_path, T.(p))
    end
    push!(full_path, T[real(right_vertex.coordinates[1]), real(right_vertex.coordinates[2])])

    is_sing = left_vertex.v_type == Singular || right_vertex.v_type == Singular

    return Edge{T}(
        id = edge_id,
        left_vertex_id = left_id,
        right_vertex_id = right_id,
        sampled_points = full_path,
        is_singular = is_sing,
    )
end

"""
    sample_edge(edge::Edge{T}, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> Edge{T}

Step 6 of the six-step framework: pure geometric resampling, with zero
`HomotopyContinuation.jl` involvement (it takes no `System`/`vertices`
argument, unlike every other function in this file). Arc-length
parametrizes `edge.sampled_points` and linearly interpolates to exactly
`cfg.edge_sample_density` equidistant points, returning a new `Edge{T}`
with the same `id`/`left_vertex_id`/`right_vertex_id`/`is_singular` and
replaced `sampled_points`.

Degenerate inputs (fewer than 2 raw points, or all raw points
coincident, e.g. a zero-length edge) are handled by simply repeating the
single available point `cfg.edge_sample_density` times, rather than
dividing by a zero total arc length.
"""
function sample_edge(edge::Edge{T}, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    pts = edge.sampled_points
    n_out = cfg.edge_sample_density

    if length(pts) < 2
        p = isempty(pts) ? T[] : copy(pts[1])
        return Edge{T}(
            id = edge.id, left_vertex_id = edge.left_vertex_id, right_vertex_id = edge.right_vertex_id,
            sampled_points = [copy(p) for _ in 1:n_out], is_singular = edge.is_singular,
        )
    end

    seglens = [norm(pts[i+1] .- pts[i]) for i in 1:(length(pts)-1)]
    cumlen = cumsum(vcat(T(0), seglens))
    total = cumlen[end]

    if total == 0
        return Edge{T}(
            id = edge.id, left_vertex_id = edge.left_vertex_id, right_vertex_id = edge.right_vertex_id,
            sampled_points = [copy(pts[1]) for _ in 1:n_out], is_singular = edge.is_singular,
        )
    end

    targets = range(T(0), total, length = n_out)
    new_pts = Vector{Vector{T}}(undef, n_out)
    seg_idx = 1
    for (k, t) in enumerate(targets)
        while seg_idx < length(cumlen) - 1 && cumlen[seg_idx+1] < t
            seg_idx += 1
        end
        t0, t1 = cumlen[seg_idx], cumlen[seg_idx+1]
        frac = t1 > t0 ? (t - t0) / (t1 - t0) : T(0)
        new_pts[k] = pts[seg_idx] .+ frac .* (pts[seg_idx+1] .- pts[seg_idx])
    end

    return Edge{T}(
        id = edge.id, left_vertex_id = edge.left_vertex_id, right_vertex_id = edge.right_vertex_id,
        sampled_points = new_pts, is_singular = edge.is_singular,
    )
end

"""
    decompose_1d_curve(F::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
        -> (vertices::Vector{NativeVertex{T}}, edges::Vector{Edge{T}})

Orchestrator for the full "1D Curve Decomposition Pipeline" (steps 3-6,
following on from Phase 2's steps 1-2). Takes the *raw* plane curve
system (same validation as [`compute_midslice`](@ref)):
`length(F.variables) == 2 && length(F.expressions) == 1`.

1. Builds the augmented critical-point system
   `System([f, differentiate(f, y_var)], [x_var, y_var])` and calls
   [`compute_critical_points`](@ref) on it; calls
   [`intersect_bounding_object`](@ref) on `F` directly.
2. Renumbers the second vertex set's ids past the first set's maximum id
   (their internal numbering otherwise collides, since each function
   starts counting at 1 independently), concatenates both sets, and
   calls `Clustering.cluster_vertices(..., cfg.vertex_match_tol)` --
   the actual cross-source "GetMergeCandidates" step (step 5).
3. Calls `Clustering.cluster_scalars` on the merged vertices' real
   x-coordinates to obtain the distinct x-slots that need connecting
   (see the Phase 3 naive-adjacency counterexample for why this step,
   rather than naively sorting+adjacent-pairing raw vertices, is
   required).
4. For each adjacent pair of distinct x-slots, calls
   [`compute_midslice`](@ref) once, then [`connect_the_dots!`](@ref)
   once per returned `y`-root, using a locally-owned incrementing
   `edge_id` counter. `vertices` is threaded through by reference so
   in-place `Artificial`-vertex insertions from one interval are visible
   to later intervals processed in the same call.
5. Runs [`sample_edge`](@ref) over every produced edge.
"""
function decompose_1d_curve(F::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    length(F.variables) == 2 && length(F.expressions) == 1 || throw(ArgumentError(
        "decompose_1d_curve: expected F with exactly 1 equation in exactly 2 variables " *
        "(a raw plane curve); got $(length(F.expressions)) equation(s) in $(length(F.variables)) variable(s).",
    ))
    x_var, y_var = F.variables
    f = F.expressions[1]

    F_aug = System([f, differentiate(f, y_var)], F.variables)
    crit_vertices = compute_critical_points(F_aug, cfg)
    bnd_vertices = intersect_bounding_object(F, cfg)

    offset = isempty(crit_vertices) ? 0 : maximum(v.id for v in crit_vertices)
    bnd_renumbered = NativeVertex{T}[
        NativeVertex{T}(id = v.id + offset, coordinates = v.coordinates, v_type = v.v_type, metadata = v.metadata)
        for v in bnd_vertices
    ]

    vertices = cluster_vertices(vcat(crit_vertices, bnd_renumbered), cfg.vertex_match_tol)

    xs = T[real(v.coordinates[1]) for v in vertices]
    distinct_xs = cluster_scalars(xs, cfg.vertex_match_tol)

    edges = Edge{T}[]
    next_edge_id = 1
    for i in 1:(length(distinct_xs)-1)
        x_left = distinct_xs[i]
        x_right = distinct_xs[i+1]
        x_mid = (x_left + x_right) / T(2)

        y_mids = compute_midslice(F, x_left, x_right, cfg)
        for y_mid in y_mids
            e = connect_the_dots!(F, x_left, x_mid, x_right, y_mid, next_edge_id, vertices, cfg)
            push!(edges, e)
            next_edge_id += 1
        end
    end

    final_edges = Edge{T}[sample_edge(e, cfg) for e in edges]
    return vertices, final_edges
end
