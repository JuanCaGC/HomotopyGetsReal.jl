# src/Types.jl
#
# Phase 1: core, precision-parametric data structures shared by every
# step of the HomotopyGetsReal pipeline (compute_critical_points,
# intersect_bounding_object, interslice/MidSlice!, ConnectTheDots!,
# Merge/GetMergeCandidates, sample_edge).
#
# These structs are intentionally minimal in Phase 1: no defaults are
# tied to `HomotopyConfig` yet (that wiring happens in Phase 2), and the
# outer constructors below only take care of building well-typed fields
# from loosely-typed keyword input (e.g. literal arrays), which is
# exactly the kind of thing that silently produced `Vector{Vector{Any}}`
# fields in the old prototype.

"""
    VertexType

Classification of a [`NativeVertex`](@ref) within the pipeline.

- `Critical`: a critical point of the projection map (from
  `compute_critical_points`).
- `Boundary`: a point where a homotopy path crosses the bounding object
  (from `intersect_bounding_object`).
- `Singular`: a point where the Jacobian is rank-deficient beyond
  `HomotopyConfig.jacobian_rank_tol`/`singular_value_threshold`.
- `Artificial`: a synthetic vertex introduced by the pipeline itself
  (e.g. a merge/averaging result from `Merge`/`GetMergeCandidates`,
  or a slicing artifact) rather than one computed directly from the
  variety.
"""
@enum VertexType begin
    Critical
    Boundary
    Singular
    Artificial
end

"""
    NativeVertex{T<:AbstractFloat}

A single vertex of the numerical decomposition, in complex space (as
homotopy paths live in `ℂⁿ` even when we ultimately only care about the
real points).

# Fields
- `id::Int`: unique identifier used to reference this vertex from
  [`Edge`](@ref)s and [`Face`](@ref)s without embedding it by value.
- `coordinates::Vector{Complex{T}}`: the vertex's coordinates.
- `v_type::VertexType`: classification, see [`VertexType`](@ref).
- `metadata::Dict{Symbol,Any}`: **intentionally untyped** grab-bag for
  per-vertex bookkeeping that varies by vertex type and pipeline stage
  (e.g. `:tolerance_used`, `:jacobian_rank`, `:source_slice`,
  `:singular_values`, ...). Making this `Dict{Symbol,Any}` rather than a
  concrete field lets every step attach whatever diagnostic/derived
  information is convenient without forcing a schema change on
  `NativeVertex` itself. This is the *only* field allowed to be
  untyped; every other field in every struct in this file is concrete
  and parametric in `T`.
"""
struct NativeVertex{T<:AbstractFloat}
    id::Int
    coordinates::Vector{Complex{T}}
    v_type::VertexType
    metadata::Dict{Symbol,Any}
end

"""
    NativeVertex{T}(; id, coordinates, v_type, metadata = Dict{Symbol,Any}())

Keyword constructor. `coordinates` may be any vector of numbers
(real or complex, any precision); it is always converted to a concrete
`Vector{Complex{T}}` so the resulting struct is fully type-stable
regardless of what was passed in.
"""
function NativeVertex{T}(;
    id::Int,
    coordinates::AbstractVector,
    v_type::VertexType,
    metadata::Dict{Symbol,Any} = Dict{Symbol,Any}(),
) where {T<:AbstractFloat}
    coords = Vector{Complex{T}}(undef, length(coordinates))
    @inbounds for i in eachindex(coordinates)
        coords[i] = Complex{T}(coordinates[i])
    end
    return NativeVertex{T}(id, coords, v_type, metadata)
end

"""
    NativeVertex(cfg::HomotopyConfig{T}, id, coordinates, v_type; metadata = Dict{Symbol,Any}())

Config-aware outer constructor (Phase 2, additive — does not replace or
alter the Phase 1 keyword constructor above). Builds a `NativeVertex{T}`
via the existing keyword constructor, then stamps
`metadata[:tolerance_used]` with whichever tolerance from `cfg` is
relevant to `v_type`:

- `Critical`  -> `cfg.critical_point_tol`
- `Boundary`  -> `cfg.boundary_point_tol`
- `Singular`  -> `cfg.singular_value_threshold`
- `Artificial` -> `cfg.vertex_match_tol` (merge/averaging tolerance)

so that every vertex produced through this constructor is
self-documenting about which numerical knob accepted it, without
callers having to thread that bookkeeping through by hand.
"""
function NativeVertex(
    cfg::HomotopyConfig{T},
    id::Int,
    coordinates::AbstractVector,
    v_type::VertexType;
    metadata::Dict{Symbol,Any} = Dict{Symbol,Any}(),
) where {T<:AbstractFloat}
    tol = if v_type === Critical
        cfg.critical_point_tol
    elseif v_type === Boundary
        cfg.boundary_point_tol
    elseif v_type === Singular
        cfg.singular_value_threshold
    else # Artificial
        cfg.vertex_match_tol
    end
    metadata = copy(metadata)
    metadata[:tolerance_used] = tol
    return NativeVertex{T}(; id = id, coordinates = coordinates, v_type = v_type, metadata = metadata)
end

"""
    Edge{T<:AbstractFloat}

A tracked path segment connecting two vertices, together with the
points sampled along it -- populated (raw) by `ConnectTheDots!`,
refined to equidistant spacing by `sample_edge`.

# Fields
- `id::Int`: unique identifier for this edge.
- `left_vertex_id::Int`, `right_vertex_id::Int`: `id`s of the
  [`NativeVertex`](@ref) endpoints this edge connects.
- `sampled_points::Vector{Vector{T}}`: real-space points sampled along
  the edge. Deliberately `Vector{Vector{T}}` (a concrete, homogeneous
  vector-of-vectors) rather than `Vector{Vector{Any}}`, which was a bug
  in the old prototype that silently boxed every coordinate and defeated
  type inference downstream.
- `is_singular::Bool`: whether this edge passes through/near a singular
  locus (as opposed to being a smooth path).
"""
struct Edge{T<:AbstractFloat}
    id::Int
    left_vertex_id::Int
    right_vertex_id::Int
    sampled_points::Vector{Vector{T}}
    is_singular::Bool
end

"""
    Edge{T}(; id, left_vertex_id, right_vertex_id, sampled_points = Vector{T}[], is_singular = false)

Keyword constructor. `sampled_points` may be any iterable of iterables
of numbers; each inner point is converted to a concrete `Vector{T}`,
and the outer container to a concrete `Vector{Vector{T}}`, so the field
never degrades to `Vector{Vector{Any}}` regardless of how the caller
built the input (e.g. from `push!`-ing onto an untyped `[]` literal).
"""
function Edge{T}(;
    id::Int,
    left_vertex_id::Int,
    right_vertex_id::Int,
    sampled_points = Vector{T}[],
    is_singular::Bool = false,
) where {T<:AbstractFloat}
    pts = Vector{Vector{T}}(undef, length(sampled_points))
    @inbounds for (i, p) in enumerate(sampled_points)
        pts[i] = Vector{T}(collect(p))
    end
    return Edge{T}(id, left_vertex_id, right_vertex_id, pts, is_singular)
end

"""
    Face{T<:AbstractFloat}

A 2D cell of the decomposition, built from a mid-slice mesh
(`interslice`/`MidSlice!`) and bounded by a set of [`Edge`](@ref)s.

# Fields
- `id::Int`: unique identifier for this face.
- `mid_slice_z::T`: the z-coordinate of the mid-slice this face was
  extracted from.
- `boundary_edges::Vector{Int}`: `id`s of the [`Edge`](@ref)s that
  bound this face.
- `mesh_vertices::Matrix{T}`: mesh vertex coordinates, one row per
  vertex.
- `mesh_topology::Matrix{Int}`: mesh connectivity (e.g. triangle/quad
  index lists), one row per mesh element.
"""
struct Face{T<:AbstractFloat}
    id::Int
    mid_slice_z::T
    boundary_edges::Vector{Int}
    mesh_vertices::Matrix{T}
    mesh_topology::Matrix{Int}
end

"""
    Face{T}(; id, mid_slice_z, boundary_edges = Int[], mesh_vertices = zeros(T, 0, 3), mesh_topology = zeros(Int, 0, 3))

Keyword constructor. `mid_slice_z` is converted via `T(...)`, and
`mesh_vertices`/`mesh_topology` are converted to concrete
`Matrix{T}`/`Matrix{Int}` so the struct is type-stable regardless of
the element type of whatever array-like object was passed in.
"""
function Face{T}(;
    id::Int,
    mid_slice_z,
    boundary_edges::AbstractVector{<:Integer} = Int[],
    mesh_vertices::AbstractMatrix = zeros(T, 0, 3),
    mesh_topology::AbstractMatrix{<:Integer} = zeros(Int, 0, 3),
) where {T<:AbstractFloat}
    return Face{T}(
        id,
        T(mid_slice_z),
        Vector{Int}(boundary_edges),
        Matrix{T}(mesh_vertices),
        Matrix{Int}(mesh_topology),
    )
end
