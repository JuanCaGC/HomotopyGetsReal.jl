
module dataModule

export find_directory, gather_vertices 

include("vertex.jl")
using .vertexModule 

function find_directory(directory_path::String)
    # Find the directory of decomposition data based on directory path
    directories = split(directory_path, '/')
    directory = length(directories) > 1 ? directories[end] : directory_path
    return directory
end

function gather_vertices(directory::String)
    vertex_file_name = "V.vertex"
    if isfile(joinpath(directory, "V_samp.vertex"))
        vertex_file_name = "V_samp.vertex"
    end

    open(joinpath(directory, vertex_file_name), "r") do f
        line = split(readline(f), ' ')
        num_vertices = parse(Int, line[1])
        num_projections = parse(Int, line[2])
        num_natural_vars_incl_hom_coord = parse(Int, line[3])
        num_variables = num_natural_vars_incl_hom_coord - 1
        num_filenames = parse(Int, strip(line[4]))
        filenames = []

        for i in 1:(num_projections * num_natural_vars_incl_hom_coord)
            skip_this_line = readline(f)
            while strip(skip_this_line) == ""
                skip_this_line = readline(f)
            end
        end

        if strip(readline(f)) == ""
            readline(f)
        end

        for ii in 1:num_filenames
            push!(filenames, strip(readline(f)))
        end

        vertices = Vector(undef, num_vertices)

        for ii in 1:num_vertices
            line = readline(f)
            while strip(line) == ""
                line = readline(f)
            end
            number_of_variables = parse(Int, line)
            
            temporary_point = []
            for jj in 1:number_of_variables
                complex_num = split(readline(f), ' ')
                real_part = parse(Float64, complex_num[1])
                imaginary_part = parse(Float64, complex_num[2])

                push!(temporary_point,complex(real_part, imaginary_part))
            end

            point = dehomogenize(temporary_point[1:num_natural_vars_incl_hom_coord])
            
            line = readline(f)
            while strip(line) == ""
                line = readline(f)
            end

            num_projection_values = parse(Int, line)
            proj = []
            for jj in 1:num_projection_values
                complex_num = split(readline(f), ' ')
                real_part = parse(Float64, complex_num[1])
                imaginary_part = parse(Float64, complex_num[2])

                push!(proj,complex(real_part, imaginary_part))
            end

            line = readline(f)
            while strip(line) == ""
                line = readline(f)
            end

            input_in = parse(Float64, line)

            line = readline(f)
            while strip(line) == ""
                line = readline(f)
            end

            vertextype = parse(Int, line)
            

            v = Vertex(point,  input_in, proj, vertextype)
            
            vertices[ii] = v

        end
        
        return vertices, filenames
    end
end

function dehomogenize(points::Vector{Any}, index::Int=0)
    # Dehomogenizes points; assumes index 0 for dehomogenization and dimension 1
    new_points = []
    for i in 2:length(points)
        push!(new_points, points[i] / points[index + 1])  # Julia is 1-indexed, so `index + 1`
    end
    return new_points
end

end