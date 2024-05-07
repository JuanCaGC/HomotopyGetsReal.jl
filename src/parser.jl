module parserModule

export  parse_directory_name,
        parse_decomposition,
        parse_edges,
        parse_curve_samples

function parse_directory_name(directory_name)
    
    if !isfile(directory_name)
        println("The File does not exist")
        return nothing
    end

    open(directory_name, "r") do f
        directory = chomp(readline(f))
        MPtype = chomp(readline(f))
        dimension = chomp(readline(f))
        return [String(directory), String(MPtype), String(dimension)]
    end
end

function parse_decomposition(directory::String)
    """
    Lee datos del archivo decomp.

    :param directory: nombre del directorio
    :return: Lista que contiene los datos a almacenar en una instancia de la clase BRinfo [pi, patch_vectors, radius, center]
    """
    if !isfile(joinpath(directory, "decomp"))
        println("did not find decomp at $(pwd())")
        return Dict()
    end

    inputFileName = ""
    num_variables_and_dimension = []
    pi = []
    patch_vectors = []
    radius = 0.0
    center = []
    num_patches = 0
    num_variables = 0
    dimension = 0

    open(joinpath(directory, "decomp"), "r") do f
        inputFileName = chomp(readline(f))
        num_variables_and_dimension = split(chomp(readline(f)), " ")
        num_variables = parse(Int, num_variables_and_dimension[1])
        dimension = parse(Int, num_variables_and_dimension[2])
        
        pi = [[0.0, 0.0] for i in 1:(num_variables-1)]
        
        for ii in 1:dimension
            numVars = readline(f)
            
            while strip(numVars) == ""
                numVars = readline(f)
            end
            numVars = parse(Int, chomp(numVars))
            for jj in 1:numVars
                pi_nums = split(chomp(readline(f)), " ")
                
                if jj == 1
                    continue
                end
                pi[jj - 1][ii] = parse(Float64, pi_nums[1]) + parse(Float64, pi_nums[2]) * im
            end
        end
    
        num_patches = readline(f)
        
        while strip(num_patches) == ""
            num_patches = readline(f)
        end
        num_patches = parse(Int, chomp(num_patches))
        patch_vectors = []
        for ii in 1:num_patches
            push!(patch_vectors, [])
            patch_size = readline(f)
            while strip(patch_size) == ""
                patch_size = readline(f)
            end
            patch_size = parse(Int, chomp(patch_size))
            for jj in 1:patch_size
                patch_vectors_data = split(chomp(readline(f)))
                push!(patch_vectors[ii], parse(Float64, patch_vectors_data[1]) + parse(Float64, patch_vectors_data[2]) * im)
            end
        end
    
        radius = readline(f)
        while strip(radius) == ""
            radius = readline(f)
        end
        radius = parse(Float64,split(chomp(radius), " ")[1])
    
        centerSize = readline(f)
        while strip(centerSize) == ""
            centerSize = readline(f)
        end
        centerSize = parse(Int64, chomp(centerSize))
        center = []
        for ii in 1:centerSize
            center_data = split(chomp(readline(f)), " ")
            push!(center, parse(Float64, center_data[1]))
        end
        dicio =  Dict(
            "input file name" => inputFileName,
            "pi info" => pi,
            "patch vectors" => patch_vectors,
            "radius" =>  radius,
            "center" => center,
            "num patches" => num_patches,
            "num_variables" => num_variables - 1,
            "dimension" => dimension
        )
        return dicio
    end
end

function parse_edges(directory)
    """ Parse and store edges data

        :param directory: Directory of the edge folder
    """
    if !isfile(joinpath(directory, "E.edge"))
        println("E.edge file not found in current directory: ", pwd())
        return Dict("number of edges" => 0, "edges" => [])
    end
    curves = Dict()
    open(joinpath(directory, "E.edge"), "r") do f
        
        curves["number of edges"] = parse(Int, chomp(readline(f)))
        
        curves["edges"] = zeros(Int, curves["number of edges"], 3)
            for ii in 1:curves["number of edges"]
                edges = readline(f)
                while strip(edges) == ""
                    edges = readline(f)
                end
                edges = split(chomp(edges), " ")
                for jj in 1:3
                    curves["edges"][ii, jj] = parse(Int, edges[jj])
                end
            end
        
    end
    return curves
end

function parse_curve_samples(directory)
    """
    Parse and store curve samples data

    :param directory: Directory of the curve folder
    """
    filename = joinpath(directory, "samp.curvesamp")
    if !isfile(filename)
        error("no samples found for this surface")
    end

    open(filename, "r") do f
        num_edges = parse(Int, readline(f))
        readline(f)  # read blank line.
        sampler_data = []

        for ii in 1:num_edges
            num_samples = parse(Int, readline(f))
            temp = []
            thing = split(readline(f))
            for jj in thing
                push!(temp, parse(Int, jj))
            end
            push!(sampler_data, temp)
            readline(f)  # read blank line.
        end

        return sampler_data
    end
end


end