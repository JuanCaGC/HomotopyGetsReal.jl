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

Cluster nearby vertices and return one representative per component.

Form connected components under Euclidean distance `≤ tol` on complex
coordinates, then merge each cluster into a single [`NativeVertex`](@ref).
Representatives are sorted by `id`.
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

    # Representative rules: smallest id; coordinate centroid; Singular wins over other types,
    # else common v_type, else Artificial; metadata via merge_metadata (nothing discarded).
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

Cluster nearby scalar values and return sorted cluster centroids.

Sort `xs`, merge consecutive values within `tol`, and return the
arithmetic mean of each run in ascending order.
"""
function cluster_scalars(xs::Vector{T}, tol::T) where {T<:AbstractFloat}
    # 1D transitive closure reduces to adjacent runs after sorting: if a ≤ b ≤ c and
    # |a-c| ≤ tol, then |a-b| and |b-c| ≤ tol too. O(n log n) sort + O(n) scan suffices
    # (unlike cluster_vertices in 2+ dimensions). Collapses duplicate x-slots before
    # decompose_1d_curve; finer (x,y) matching is left to connect_the_dots!.
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

Cluster nearby points and return representatives plus a membership map.

Apply tolerance-based union-find on real coordinate vectors (same distance
rule as [`cluster_vertices`](@ref)). Return centroid representatives and
`membership[i]`, the cluster index for each input point.
"""
function cluster_points_indexed(points::Vector{Vector{T}}, tol::T) where {T<:AbstractFloat}
    # Used by SurfaceDecomposition.weld_mesh: unlike cluster_vertices, operates on raw
    # Vector{Vector{T}} mesh rows and returns membership so mesh_topology indices can be
    # remapped to a global vertex numbering. Representatives ordered by first encounter.
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
