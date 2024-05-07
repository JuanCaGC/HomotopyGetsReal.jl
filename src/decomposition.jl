module decompositionModule
    export Decomposition, parse_decomp, read_input, extract_points

include("parser.jl")
using .parserModule

include("data.jl")
using .dataModule 

mutable struct Decomposition
    directory::String
    is_embedded::Bool
    embedded_into::Union{Nothing, Decomposition}

    input::Union{Nothing, String}
    inputfilename::Union{Nothing, String}
    num_variables::Int
    pi::Array{Any, 1}
    num_patches::Int
    patch::Array{Any, 1}
    radius::Float64
    center_size::Int
    center::Array{Float64, 1}
    dimension::Int

    _memoized_data::Dict{String, Any}
    vertices::Array{Any, 1}
    filenames::Array{Any, 1}

    function Decomposition(directory::String ; is_embedded::Bool=false, embedded_into::Union{Nothing, Decomposition}=nothing)
        decomp = new()

        decomp.directory = directory
        decomp.is_embedded = is_embedded
        decomp.embedded_into = embedded_into

        decomp.input = nothing
        decomp.inputfilename = nothing
        decomp.num_variables = 0
        decomp.pi = []
        decomp.num_patches = 0
        decomp.patch = []
        decomp.radius = 0
        decomp.center_size = 0
        decomp.center = []
        decomp.dimension = 0

        parse_decomp(decomp, directory)
        read_input(decomp, directory)

        decomp._memoized_data = Dict()

        if !decomp.is_embedded
           
            decomp.vertices, decomp.filenames = gather_vertices(directory)
        else
            if decomp.embedded_into === nothing
                throw(EmbeddedIssue("parameter `embedded_into` cannot be unset if the decomposition is embedded"))
            end
            
            decomp.vertices = decomp.embedded_into.vertices
        end

        return decomp
    end
end

function parse_decomp(decomp, directory::String)
    
    decomposition_data = parse_decomposition(directory)
    decomp.inputfilename = decomposition_data["input file name"]
    decomp.pi = decomposition_data["pi info"]
    decomp.patch = decomposition_data["patch vectors"]
    decomp.radius = decomposition_data["radius"]
    decomp.center = decomposition_data["center"]
    decomp.num_patches = decomposition_data["num patches"]
    decomp.num_variables = decomposition_data["num_variables"]
    decomp.dimension = decomposition_data["dimension"]
end

function read_input(decomp, directory::String)
    filename = joinpath(directory, decomp.inputfilename)
    if !isfile(filename)
        println("Could not find input file in current directory: %s\n", directory)
    else
        decomp.input = read(filename, String)
    end
end

function extract_points(decomp; indices::Union{Nothing, Vector{Int}}=nothing)
    if indices === nothing
        indices = 1:length(decomp.vertices)
    end

    if decomp.is_embedded
        return decomp.embedded_into.extract_points()
    end

    if  "points" in keys(decomp._memoized_data)
        return decomp._memoized_data["points"]
    end

    points = []

    for ii in indices
        vertex = decomp.vertices[ii]

        point = fill(NaN, decomp.num_variables)

        for jj in 1:decomp.num_variables
            point[jj] = real(vertex.point[jj])
        end
        push!(points, point)
    end

    points = hcat(points...)

    decomp._memoized_data["points"] = points

    return points
end

end