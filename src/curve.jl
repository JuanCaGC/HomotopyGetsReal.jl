module curveModule
export Curve
include("parser.jl")
using .parserModule

include("decomposition.jl")
using .decompositionModule

struct CurvePiece
    curve
    edge_indices
    directed_edges

    function CurvePiece(curve, edges, directed_edges)
        new(curve, edges, directed_edges)
    end

    function __repr__(self)
        s = "CurvePiece on curve $(self.curve.inputfilename), indices $(self.edge_indices), directed_edges $(self.directed_edges)\n"
        return s
    end

    function __str__(self)
        return __repr__(self)
    end

    function to_points(self)
        if isempty(self.directed_edges)
            error("Insert code memoizing / computing the directed edges")
        end

        vertices = self.curve.vertices
        points = Array{Float64}(undef, 0, 3)

        prev_point_index = -1
        for direct_edge in self.directed_edges
            edge_index = direct_edge.edge_index
            direction = direct_edge.direction
            edge = self.curve.edges[edge_index]
            if is_edge_degenerate(edge)
                continue
            end

            point_indices = length(self.curve.sampler_data) > 0 ? self.curve.sampler_data[edge_index] : edge

            if direction == EdgeDirection(backward)
                reverse!(point_indices)
            end

            list_of_points = []
            for ii in point_indices
                if ii != prev_point_index
                    push!(list_of_points, vertices[ii].point)
                    prev_point_index = ii
                end
            end

            points_this_edge = hcat(list_of_points...)'  # Transpose to match dimensions
            points = vcat(points, points_this_edge)
        end

        return points
    end
end

function EdgeDirection(direction)
    if direction == "forward"
        return 1
    elseif  direction == "backward"
        return 0
    else 
        return -1
    end
end

struct DirectedEdge
    edge_index
    direction

    function DirectedEdge(edge_index, direction)
        new(edge_index, direction)
    end
end

function is_edge_degenerate(e)
    return (e[1] == e[2]) || (e[2] == e[3])
end

mutable struct Curve 
    decomposition :: Decomposition
    directory::String
    is_embedded::Bool
    embedded_into::Union{Nothing, Decomposition}

    num_edges::Int
    edges
    sampler_data

    function Curve(directory; is_embedded=false, embedded_into=nothing)
        c = new()
        c.directory = directory
        c.is_embedded = is_embedded
        c.embedded_into = embedded_into
        c.num_edges = 0
        c.edges = []
        c.sampler_data = nothing

        c.decomposition = Decomposition(directory; is_embedded, embedded_into)
        parse_edge(c, directory)
        try
            parse_curve_sample(c, directory)
        catch e
            if occursin("no samples found for this surface", e.msg) && !c.is_embedded
                println("No samples to gather")
            else
                rethrow(e)
            end
        end

        return c
    end
end

function parse_edge(curve, directory)
    edge_data = parse_edges(directory)  
    curve.num_edges = edge_data["number of edges"]
    curve.edges = edge_data["edges"]
end

function parse_curve_sample(curve, directory)
    curve.sampler_data = parse_curve_samples(directory)
end

end