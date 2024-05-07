module WitnessSetModule
export WitnessSet
abstract type Holder end

mutable struct PointHolder <: Holder
    points::Vector{Array{Float64, 1}}
end

mutable struct LinearHolder <: Holder
    linears::Vector{Array{Float64, 1}}
end

mutable struct PatchHolder <: Holder
    patches::Vector{Array{Float64, 1}}
end

mutable struct NameHolder <: Holder
    names::Vector{String}
end

mutable struct Function
    func::String
end

mutable struct WitnessSet 
    dim::Int
    comp_num::Int
    incid_num::Int
    num_vars::Int
    num_natty_vars::Int
    input_filename::String
    input_file::Function

    WitnessSet(nvar) = new(0, 0, 0, nvar, nvar, "unset_filename", Function(""))
end

function input_filename(ws::WitnessSet)
    return ws.input_filename
end

function set_input_filename!(ws::WitnessSet, new_input_filename)
    ws.input_filename = new_input_filename
end

function dimension(ws::WitnessSet)
    return ws.dim
end

function set_dimension!(ws::WitnessSet, new_dim)
    ws.dim = new_dim
end

function component_number(ws::WitnessSet)
    return ws.comp_num
end

function set_component_number!(ws::WitnessSet, new_comp_num)
    ws.comp_num = new_comp_num
end

function num_variables(ws::WitnessSet)
    return ws.num_vars
end

function set_num_variables!(ws::WitnessSet, new_num_vars)
    ws.num_vars = new_num_vars
end

function num_natural_variables(ws::WitnessSet)
    return ws.num_natty_vars
end

function set_num_natural_variables!(ws::WitnessSet, new_num_nat_vars)
    ws.num_natty_vars = new_num_nat_vars
end

function incidence_number(ws::WitnessSet)
    return ws.incid_num
end

function set_incidence_number!(ws::WitnessSet, new_incidence)
    ws.incid_num = new_incidence
end

function num_synth_vars(ws::WitnessSet)
    return ws.num_vars - ws.num_natty_vars
end

function only_natural_vars(ws::WitnessSet)
    # implementar según sea necesario
end

function only_first_vars(ws::WitnessSet, num_vars)
    # implementar según sea necesario
end

function sort_for_real(ws::WitnessSet, tol)
    # implementar según sea necesario
end

function sort_for_unique(ws::WitnessSet, tol)
    # implementar según sea necesario
end

function sort_for_inside_sphere(ws::WitnessSet, radius, center)
    # implementar según sea necesario
end

function RealifyPoint(ws::WitnessSet, ind, tol)
    # implementar según sea necesario
end

function Realify(ws::WitnessSet, tol)
    # implementar según sea necesario
end

function RescaleToPatch(ws::WitnessSet, patch)
    # implementar según sea necesario
end

function RealifyPatches(ws::WitnessSet)
    # implementar según sea necesario
end

function Parse(ws::WitnessSet, witness_set_file, num_vars)
    # implementar según sea necesario
end

function reset!(ws::WitnessSet)
    # implementar según sea necesario
end

function clear!(ws::WitnessSet)
    # implementar según sea necesario
end

function merge!(ws::WitnessSet, W_in, tol)
    # implementar según sea necesario
end

function print_to_screen(ws::WitnessSet, dehom_points=true, print_extras=false)
    # implementar según sea necesario
end

function print_to_file(ws::WitnessSet, filename)
    # implementar según sea necesario
end

function write_linears(ws::WitnessSet, filename)
    # implementar según sea necesario
end

function print_patches(ws::WitnessSet, filename)
    # implementar según sea necesario
end

function read_patches_from_file(ws::WitnessSet, filename)
    # implementar según sea necesario
end

function write_homogeneous_coordinates(ws::WitnessSet, filename)
    # implementar según sea necesario
end

function write_dehomogenized_coordinates(ws::WitnessSet, filename)
    # implementar según sea necesario
end

function write_dehomogenized_coordinates(ws::WitnessSet, filename, indices)
    # implementar según sea necesario
end


end