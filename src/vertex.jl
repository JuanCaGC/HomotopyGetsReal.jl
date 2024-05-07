module vertexModule

export Vertex, is_of_type, VertexType
# Define an Enum for VertexType


VertexType = Dict(
        "unset" => 0,
        "critical" => 1,
        "semicritical" => 2,
        "midpoint" => 4,
        "isolated" => 8,
        "new" => 16,
        "curve_sample_point" => 32,
        "surface_sample_point" => 64,
        "removed" => 128,
        "problematic" => 256,
        "singular" => 512
    )

# Define a Vertex object
mutable struct Vertex
    point::Any 
    input_filename_index::Int
    projection_value::Any
    type::Any

    # Constructor
    function Vertex(point, input_filename_index, projection_value, vertex_type::Int)
        new(point, input_filename_index, projection_value, vertex_type)
    end
end

# Representation of Vertex object (toString method)
Base.show(io::IO, v::Vertex) = print(io, "Vertex($(v.point), $(v.input_filename_index), $(v.type))")

# Function to check if a vertex matches a certain VertexType
function is_of_type(v::Any, type::Any)
    return (v.type & type) != 0
end

end