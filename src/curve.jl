include("parser.jl")
using .parser

include("decomposition.jl")
using .decomposition : Decomposition

module curve


struct DirectedEdge
    edge_index::Int
    direction::Int
end

@enum EdgeDirection forward=1 backward=0

mutable struct CurvePiece
    curve::Curve
    edge_indices::Vector{Int}
    directed_edges::Vector{DirectedEdge}
end

function CurvePiece(curve::Curve, edges::Vector{Int}, directed_edges::Vector{DirectedEdge})
    return CurvePiece(curve, edges, directed_edges)
end

function inputfilename(piece::CurvePiece)
    return piece.curve.inputfilename
end

function Base.repr(piece::CurvePiece)
    return "CurvePiece on curve $(piece.curve.inputfilename), indices $(piece.edge_indices), directed_edges $(piece.directed_edges)\n"
end

Base.show(io::IO, piece::CurvePiece) = print(io, repr(piece))

function to_points(piece::CurvePiece)
    if isempty(piece.directed_edges)
        throw(NotImplementedError("Insert code memoizing / computing the directed edges"))
    end
    
    vertices = piece.curve.vertices
    c = piece.curve
    
    points = zeros(Float64, 0, 3)
    prev_point_index = -1
    
    for (edge_index, direction) in piece.directed_edges
        if is_edge_degenerate(c.edges[edge_index])
            continue
        end
        
        point_indices = isempty(c.sampler_data) ? c.edges[edge_index] : c.sampler_data[edge_index]
        if direction == 0 # backward
            point_indices = reverse(point_indices)
        end
        
        list_of_points = []
        for ii in point_indices
            if ii != prev_point_index
                push!(list_of_points, real(vertices[ii].point))
                prev_point_index = ii
            end
        end
        
        points_this_edge = vcat(list_of_points...)
        points = vcat(points, points_this_edge)
    end
    
    return points
end

function is_edge_degenerate(e::Vector{Int})
    return (e[1] == e[2]) || (e[2] == e[3])
end

struct Curve <: Decomposition
    directory::String
    is_embedded::Bool
    embedded_into::Union{Nothing, Curve}
    num_edges::Int
    edges::Vector{Vector{Int}}
    sampler_data::Union{Nothing, Vector{Vector{Int}}}
    
    function Curve(directory::String; is_embedded::Bool=false, embedded_into::Union{Nothing, Curve}=nothing)
        num_edges = 0
        edges = []
        sampler_data = nothing
        new(directory, is_embedded, embedded_into, num_edges, edges, sampler_data)
    end
end

function Curve(directory::String, is_embedded::Bool=false, embedded_into::Union{Nothing, Curve}=nothing)
    curve = Curve(directory, is_embedded, embedded_into)
    parse_edge(curve, directory)
    try
        parse_curve_samples(curve, directory)
    catch e
        if !curve.is_embedded
            println("no samples to gather")
        end
    end
    return curve
end

function break_into_pieces(curve::Curve, edge_indices::Union{Nothing, Vector{Int}}=nothing)
    unsorted_edges = deepcopy(edge_indices)
    list_of_pieces = Vector{CurvePiece}()

    while !isempty(unsorted_edges)
        directed_edges = [DirectedEdge(popfirst!(unsorted_edges), EdgeDirection.forward)]
        added_edge_indicator = true

        while added_edge_indicator
            added_edge_indicator = false
            first_edge = curve.edges[directed_edges[1].edge_index]
            first_edge_direction = directed_edges[1].direction
            first_point_index = first_edge[first_edge_direction == EdgeDirection.forward ? 1 : end]

            last_edge = curve.edges[directed_edges[end].edge_index]
            last_edge_direction = directed_edges[end].direction
            last_point_index = last_edge[last_edge_direction == EdgeDirection.forward ? end : 1]

            used_edges_this = Set{Int}()

            for edge_ind in unsorted_edges
                if curve.edges[edge_ind][1] == last_point_index
                    push!(used_edges_this, edge_ind)
                    push!(directed_edges, DirectedEdge(edge_ind, EdgeDirection.forward))
                elseif curve.edges[edge_ind][1] == first_point_index
                    push!(used_edges_this, edge_ind)
                    unshift!(directed_edges, DirectedEdge(edge_ind, EdgeDirection.backward))
                elseif curve.edges[edge_ind][end] == last_point_index
                    push!(used_edges_this, edge_ind)
                    push!(directed_edges, DirectedEdge(edge_ind, EdgeDirection.backward))
                elseif curve.edges[edge_ind][end] == first_point_index
                    push!(used_edges_this, edge_ind)
                    unshift!(directed_edges, DirectedEdge(edge_ind, EdgeDirection.forward))
                end
            end

            added_edge_indicator = !isempty(used_edges_this)
            unsorted_edges = setdiff(unsorted_edges, used_edges_this)
        end

        push!(list_of_pieces, CurvePiece(curve, [e.edge_index for e in directed_edges], directed_edges))
    end

    return list_of_pieces
end

function parse_edge(curve::Curve, directory::String)
    edge_data = bertini_real.parse_edges(directory)
    curve.num_edges = edge_data["number of edges"]
    curve.edges = edge_data["edges"]
end

function parse_curve_samples(curve::Curve, directory::String)
    curve.sampler_data = bertini_real.parse_curve_samples(directory)
end

function Base.show(io::IO, curve::Curve)
    println(io, "curve with:")
    println(io, curve.num_edges, " edges")
end
end