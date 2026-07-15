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

Classify a [`NativeVertex`](@ref) by how it was produced.

- `Critical`: critical point of the projection map.
- `Boundary`: intersection with the bounding object.
- `Singular`: rank-deficient Jacobian at the point.
- `Artificial`: synthetic vertex from merging or slicing.
"""
@enum VertexType begin
    Critical
    Boundary
    Singular
    Artificial
end

"""
    NativeVertex{T<:AbstractFloat}

A vertex of the numerical decomposition with complex coordinates.

Vertices are referenced by `id` from [`Edge`](@ref) and [`Face`](@ref).
Coordinates live in complex space because homotopy paths are tracked in `ℂⁿ`.

# Fields

- `id`: Unique identifier.
- `coordinates`: Vertex coordinates as `Vector{Complex{T}}`.
- `v_type`: Classification; see [`VertexType`](@ref).
- `metadata`: Per-vertex diagnostics and bookkeeping.
"""
struct NativeVertex{T<:AbstractFloat}
    id::Int
    coordinates::Vector{Complex{T}}
    v_type::VertexType
    # metadata is Dict{Symbol,Any} so pipeline stages can attach varying diagnostics
    # (e.g. :jacobian_rank, :singular_values) without schema changes on NativeVertex.
    metadata::Dict{Symbol,Any}
end

"""
    NativeVertex{T}(; id, coordinates, v_type, metadata = Dict{Symbol,Any}())

Construct a [`NativeVertex`](@ref) with type-stable complex coordinates.

Convert any numeric `coordinates` vector to `Vector{Complex{T}}`.
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

Construct a [`NativeVertex`](@ref) and record the accepting tolerance in metadata.

Set `metadata[:tolerance_used]` from the tolerance in `cfg` that matches `v_type`.
"""
function NativeVertex(
    cfg::HomotopyConfig{T},
    id::Int,
    coordinates::AbstractVector,
    v_type::VertexType;
    metadata::Dict{Symbol,Any} = Dict{Symbol,Any}(),
) where {T<:AbstractFloat}
    # Critical -> critical_point_tol; Boundary -> boundary_point_tol;
    # Singular -> singular_value_threshold; Artificial -> vertex_match_tol.
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

A path segment between two vertices with sampled points along it.

Raw samples come from `connect_the_dots!`; `sample_edge` refines spacing.

# Fields

- `id`: Unique edge identifier.
- `left_vertex_id`, `right_vertex_id`: Endpoint [`NativeVertex`](@ref) ids.
- `sampled_points`: Real-space samples as `Vector{Vector{T}}`.
- `is_singular`: Whether the edge passes near a singular locus.
"""
struct Edge{T<:AbstractFloat}
    id::Int
    left_vertex_id::Int
    right_vertex_id::Int
    # Vector{Vector{T}} keeps samples homogeneous; Vector{Vector{Any}} boxed coordinates
    # in the old prototype and defeated type inference downstream.
    sampled_points::Vector{Vector{T}}
    is_singular::Bool
end

"""
    Edge{T}(; id, left_vertex_id, right_vertex_id, sampled_points = Vector{T}[], is_singular = false)

Construct an [`Edge`](@ref) with concrete `Vector{Vector{T}}` samples.

Convert any nested numeric iterables to homogeneous `Vector{T}` rows.
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

A 2D decomposition cell from a mid-slice with boundary edges and mesh data.

# Fields

- `id`: Unique face identifier.
- `mid_slice_z`: Z-coordinate of the source mid-slice.
- `boundary_edges`: Ids of bounding [`Edge`](@ref)s.
- `mesh_vertices`: Mesh vertex coordinates, one row per vertex.
- `mesh_topology`: Mesh connectivity indices, one row per element.
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

Construct a [`Face`](@ref) with type-stable matrix fields.

Convert scalar and array inputs to concrete `T` and `Matrix` types.
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
