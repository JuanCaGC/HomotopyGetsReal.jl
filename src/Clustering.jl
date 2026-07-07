# src/Clustering.jl
#
# Phase 2: pure geometric deduplication/merge utility operating on
# NativeVertex{T} (not raw coordinate vectors, unlike the old
# prototype's `cluster_points` in prototipo_viejo_julia/Solver.jl).
#
# This is deliberately a separate file with a narrow dependency surface
# (Types.jl + LinearAlgebra only — no HomotopyContinuation.jl) so it
# stays a standalone geometric primitive. `compute_critical_points` and
# `intersect_bounding_object` in src/Solver.jl call `cluster_vertices`
# for their Phase 2 dedup needs, and step 5 of the six-step framework
# (Merge / GetMergeCandidates) is expected to reuse or extend this same
# primitive in a later phase — hence the API is designed generically
# around "a tolerance-indexed clustering of vertices with metadata-
# preserving merge" rather than anything specific to Phase 2's two
# callers.

using LinearAlgebra

"""
    cluster_vertices(vertices::Vector{NativeVertex{T}}, tol::T) where {T<:AbstractFloat}
        -> Vector{NativeVertex{T}}

Partitions `vertices` into clusters by transitive closure of the
relation "within `tol` of each other" (i.e. connected components of the
graph where an edge connects any two vertices whose Euclidean
coordinate distance is `<= tol`), and returns one representative
`NativeVertex{T}` per cluster.

This is a general clustering *primitive*: it does not know anything
about `Critical`/`Boundary`/`Singular` classification semantics beyond
the conflict-resolution rule documented below, so it is safe for both
Phase 2's straightforward dedup use and for a future, more elaborate
`Merge`/`GetMergeCandidates` step (6-step framework, step 5) to build on
top of.

# Representative construction
For each cluster:
- `id`: the smallest `id` among cluster members (deterministic and
  stable regardless of input order).
- `coordinates`: the componentwise centroid (arithmetic mean) of all
  member coordinates.
- `v_type`: `Singular` if *any* member is `Singular` (a numerically
  detected singularity is never allowed to be "averaged away" by
  merging with a nearby regular point); else the common `v_type` if all
  members agree; else `Artificial` to flag that this vertex is a
  synthetic merge of members with genuinely different classifications.
- `metadata`: built by [`merge_metadata`](@ref) — see its docstring for
  the exact merge rule. Nothing from any member is discarded.

Distance is computed with `LinearAlgebra.norm` on the (complex)
coordinate difference, so this works uniformly for `T = Float64` and
`T = BigFloat`.
"""
function cluster_vertices(vertices::Vector{NativeVertex{T}}, tol::T) where {T<:AbstractFloat}
    n = length(vertices)
    n == 0 && return NativeVertex{T}[]

    parent = collect(1:n)
    find(i) = begin
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        i
    end
    function union!(i, j)
        ri, rj = find(i), find(j)
        ri != rj && (parent[ri] = rj)
        return nothing
    end

    @inbounds for i in 1:n, j in (i+1):n
        if norm(vertices[i].coordinates .- vertices[j].coordinates) <= tol
            union!(i, j)
        end
    end

    clusters = Dict{Int,Vector{Int}}()
    for i in 1:n
        root = find(i)
        push!(get!(clusters, root, Int[]), i)
    end

    reps = Vector{NativeVertex{T}}(undef, length(clusters))
    for (k, member_idxs) in enumerate(values(clusters))
        reps[k] = _merge_cluster(vertices, member_idxs)
    end
    sort!(reps; by = v -> v.id)
    return reps
end

function _merge_cluster(vertices::Vector{NativeVertex{T}}, member_idxs::Vector{Int}) where {T<:AbstractFloat}
    members = @view vertices[member_idxs]

    rep_id = minimum(v.id for v in members)

    centroid = sum(v.coordinates for v in members) ./ T(length(members))

    v_type = if any(v.v_type === Singular for v in members)
        Singular
    else
        types = unique(v.v_type for v in members)
        length(types) == 1 ? types[1] : Artificial
    end

    metadata = merge_metadata([v.metadata for v in members], [v.id for v in members])

    return NativeVertex{T}(
        id = rep_id,
        coordinates = centroid,
        v_type = v_type,
        metadata = metadata,
    )
end

"""
    merge_metadata(metadatas::Vector{Dict{Symbol,Any}}, ids::Vector{Int})
        -> Dict{Symbol,Any}

Combines the `metadata` dictionaries of a set of clustered vertices
into one, without discarding any of the original per-member
information:

- `:cluster_size => length(metadatas)`
- `:cluster_member_ids => sort(ids)`
- For every key present in *every* member's metadata, [`_combine_values`](@ref)
  decides how to combine it (see its docstring for the exact rule --
  in particular, **integer-valued fields such as `:jacobian_rank` are
  never averaged**, since e.g. `(2 + 2 + 1) / 3 = 1.667` is not a valid
  rank; only genuinely continuous `AbstractFloat` fields like
  `:singular_values` or `:tolerance_used` are).
- Keys not shared by every member are kept individually as
  `Symbol("<key>_member<i>") => value` so partial/diagnostic
  information from a subset of members still survives the merge.
"""
function merge_metadata(metadatas::Vector{Dict{Symbol,Any}}, ids::Vector{Int})
    merged = Dict{Symbol,Any}()
    merged[:cluster_size] = length(metadatas)
    merged[:cluster_member_ids] = sort(ids)

    shared_keys = length(metadatas) == 0 ? Set{Symbol}() :
        reduce(intersect, (Set(keys(m)) for m in metadatas))

    for key in shared_keys
        values_ = [m[key] for m in metadatas]
        merged[key] = _combine_values(key, values_)
    end

    for (i, m) in enumerate(metadatas)
        for (key, val) in m
            key in shared_keys && continue
            merged[Symbol(key, :_member, i)] = val
        end
    end

    return merged
end

"""
    _combine_values(key::Symbol, values_::Vector)

Merge rule used by [`merge_metadata`](@ref) for a single metadata key
shared by every clustered member, in priority order:

1. `key === :jacobian_rank` (and every value is an `Integer`): take the
   **minimum** rank across members, not an average -- ranks are
   discrete, and `(2 + 2 + 1) / 3 = 1.667` is not a valid rank. The
   minimum is the conservative choice: it's consistent with the
   `Singular`-wins policy already used for `v_type` resolution in
   [`cluster_vertices`](@ref) (a more rank-deficient member should not
   have its degeneracy averaged away by nearby regular members).
2. Every value is `<:AbstractFloat` (genuinely continuous, e.g.
   `:tolerance_used`): arithmetic mean.
3. Every value is an `AbstractVector{<:AbstractFloat}` of the same
   length (e.g. `:singular_values`): elementwise arithmetic mean.
4. All values are `==`-equal (covers exact-duplicate `Integer`/`Symbol`/
   etc. fields that aren't otherwise handled above, e.g. `:fixed_variable`):
   keep that shared value as-is.
5. Otherwise (including `Integer`-valued fields other than
   `:jacobian_rank` that genuinely disagree, where we have no domain
   knowledge about whether averaging is meaningful): keep every
   member's value, unreduced, as a `Vector` -- nothing is silently
   discarded even when it can't be sensibly combined.
"""
function _combine_values(key::Symbol, values_::Vector)
    if key === :jacobian_rank && all(v -> v isa Integer, values_)
        return minimum(values_)
    elseif all(v -> v isa AbstractFloat, values_)
        return sum(values_) / length(values_)
    elseif all(v -> v isa AbstractVector{<:AbstractFloat}, values_) && length(unique(length.(values_))) == 1
        return sum(values_) ./ length(values_)
    elseif all(==(values_[1]), values_)
        return values_[1]
    else
        return values_
    end
end

"""
    cluster_scalars(xs::Vector{T}, tol::T) where {T<:AbstractFloat}
        -> Vector{T}

Sibling primitive to [`cluster_vertices`](@ref), for the Phase 3
"which x-values are actually distinct" problem: same union-find-by-
tolerance idea as `cluster_vertices`, collapsed to 1D (an edge connects
`x_i` and `x_j` iff `abs(x_i - x_j) <= tol`), returning one
representative (arithmetic mean) per connected component, **sorted**
ascending.

# Why a simple sorted adjacent-chain scan is enough (not O(n^2) union-find)
For 1D data, if `a <= b <= c` and `abs(a - c) <= tol`, then necessarily
`abs(a - b) <= tol` and `abs(b - c) <= tol` too (since `b` lies between
`a` and `c`). So any edge between two non-adjacent elements of the
*sorted* sequence already implies edges between every intermediate
adjacent pair. Consequently the connected components of the full
"within `tol`" graph are exactly the maximal runs of the sorted sequence
whose consecutive gaps are all `<= tol` -- an `O(n log n)` sort plus a
single `O(n)` scan, with no loss of generality relative to a full
pairwise union-find (which is what [`cluster_vertices`](@ref) must use
instead, since that equivalence does not hold in 2+ dimensions).

This is exactly the fix for the naive-adjacency counterexample: given
the nodal cubic's boundary vertices at `x ≈ 2.2268` (once at `y = -4`,
once at `y = +4`), `cluster_scalars` collapses them into a single
x-slot *before* `decompose_1d_curve` decides which intervals need a
`compute_midslice` call, so the pipeline never mistakes "same x" for
"same vertex" (that finer distinction is left to
`connect_the_dots!`'s full `(x, y)` coordinate matching).
"""
function cluster_scalars(xs::Vector{T}, tol::T) where {T<:AbstractFloat}
    n = length(xs)
    n == 0 && return T[]

    sorted_xs = sort(xs)
    reps = T[]
    cluster_sum = sorted_xs[1]
    cluster_count = 1
    for i in 2:n
        if sorted_xs[i] - sorted_xs[i-1] <= tol
            cluster_sum += sorted_xs[i]
            cluster_count += 1
        else
            push!(reps, cluster_sum / T(cluster_count))
            cluster_sum = sorted_xs[i]
            cluster_count = 1
        end
    end
    push!(reps, cluster_sum / T(cluster_count))
    return reps
end

"""
    cluster_points_indexed(points::Vector{Vector{T}}, tol::T) where {T<:AbstractFloat}
        -> (representatives::Vector{Vector{T}}, membership::Vector{Int})

Phase 5 primitive, needed by `SurfaceDecomposition.weld_mesh` for a job
neither [`cluster_vertices`](@ref) nor [`cluster_scalars`](@ref) can do:
`cluster_vertices` operates on `NativeVertex{T}`'s `Complex{T}`
coordinates (mesh vertices are plain real `T`-valued rows of a
`Face.mesh_vertices` matrix, not `NativeVertex`), and -- more
fundamentally -- it returns only deduplicated *representatives*, never
the mapping from each original input to the cluster it merged into.
`weld_mesh` needs exactly that mapping to remap every face's
`mesh_topology` (local row indices) into a single global vertex
numbering (the mesh-index analogue of Phase 3's cross-source id
renumbering, see `SurfaceDecomposition.weld_mesh`'s own docstring).

Same tolerance-indexed union-find-by-Euclidean-distance idea as
`cluster_vertices` (an edge connects `points[i]` and `points[j]` iff
`norm(points[i] .- points[j]) <= tol`, `O(n^2)` pairwise, exactly like
`cluster_vertices` -- no `NativeVertex`-specific metadata/`v_type`
merge logic is needed here, so this is a strictly simpler sibling, not
a generalization of `cluster_vertices`), but reports full provenance
instead of discarding it:

- `representatives[k]` is the componentwise centroid (arithmetic mean)
  of every `points[i]` in cluster `k`, in the order clusters are first
  encountered while scanning `1:n` (deterministic given input order,
  unlike `cluster_vertices`'s id-based sort -- there is no id concept
  for raw points).
- `membership[i]` is the index into `representatives` that `points[i]`
  merged into, for every `i` in `1:length(points)`, so callers can
  translate any original per-point reference (e.g. a local mesh-grid
  row index) into the deduplicated global numbering.
"""
function cluster_points_indexed(points::Vector{Vector{T}}, tol::T) where {T<:AbstractFloat}
    n = length(points)
    n == 0 && return Vector{T}[], Int[]

    parent = collect(1:n)
    find(i) = begin
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        i
    end
    function union!(i, j)
        ri, rj = find(i), find(j)
        ri != rj && (parent[ri] = rj)
        return nothing
    end

    @inbounds for i in 1:n, j in (i+1):n
        if norm(points[i] .- points[j]) <= tol
            union!(i, j)
        end
    end

    root_to_cluster = Dict{Int,Int}()
    clusters = Vector{Int}[]
    for i in 1:n
        root = find(i)
        cluster_idx = get!(root_to_cluster, root) do
            push!(clusters, Int[])
            length(clusters)
        end
        push!(clusters[cluster_idx], i)
    end

    representatives = Vector{Vector{T}}(undef, length(clusters))
    membership = Vector{Int}(undef, n)
    for (k, member_idxs) in enumerate(clusters)
        representatives[k] = sum(points[i] for i in member_idxs) ./ T(length(member_idxs))
        for i in member_idxs
            membership[i] = k
        end
    end

    return representatives, membership
end
