module Types
    export NativeVertex, Edge, Face, VertexType, Critical, Singular, Boundary

    @enum VertexType Critical Singular Boundary

    struct NativeVertex
        id::Int
        coordinates::Vector{ComplexF64} # Now [x, y, z]
        v_type::VertexType
    end

    struct Edge
        id::Int
        left_id::Int
        path::Vector{Vector{Float64}} # Points in 3D
        right_id::Int
    end

    struct Face
        id::Int
        mid_z::Float64
        # A face is bounded by a set of edges from the slices above and below
        boundary_edges::Vector{Int} 
        # For visualization, we will eventually store a mesh (Triangles)
        mesh_vertices::Matrix{Float64}
        mesh_faces::Matrix{Int}
    end
end