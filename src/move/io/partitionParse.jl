module PartitionParser

export partition_parse

function partition_parse(input_filepath::String, config_filepath::String, function_filepath::String)
    # Abrir el archivo de entrada y los archivos de salida
    input_file = open(input_filepath, "r")
    config_file = open(config_filepath, "w")
    function_file = open(function_filepath, "w")

    # Estado para saber dónde escribir
    in_config_section = false
    in_function_section = false

    try
        # Leer línea por línea del archivo de entrada
        for line in eachline(input_file)
            # Comprobar marcas de inicio y fin de secciones
            if occursin("CONFIG", line)
                in_config_section = true
                in_function_section = false
            elseif occursin("INPUT", line)
                in_config_section = false
                in_function_section = true
            elseif occursin("END", line)
                in_config_section = false
                in_function_section = false
            end

            # Escribir en el archivo correspondiente
            if in_config_section
                write(config_file, line * "\n")
            elseif in_function_section
                write(function_file, line * "\n")
            end
        end
    finally
        # Asegurarse de cerrar todos los archivos al final
        close(input_file)
        close(config_file)
        close(function_file)
    end
end

end
