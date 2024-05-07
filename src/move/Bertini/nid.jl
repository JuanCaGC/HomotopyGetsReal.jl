using LinearAlgebra

mutable struct WitnessPointMetadata
    dimension::Int
    corank::Int
    typeflag::Int
    multiplicity::Int
    component_number::Int
    deflations_needed::Int
    condition_number::Float64
    smallest_nonzero_sing_value::Float64
    largest_zero_sing_value::Float64

    function WitnessPointMetadata(new_dim::Int)
        return new(new_dim, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0)
    end
end

mutable struct WitnessLinearMetadata
    dim::Int

    function WitnessLinearMetadata(new_dim::Int)
        return new(new_dim)
    end
end

mutable struct WitnessPatchMetadata
    dim::Int

    function WitnessPatchMetadata(new_dim::Int)
        return new(new_dim)
    end
end

mutable struct NumericalIrreducibleDecomposition
    point_metadata::Vector{WitnessPointMetadata}
    linear_metadata::Vector{WitnessLinearMetadata}
    patch_metadata::Vector{WitnessPatchMetadata}
    nonempty_dimensions::Vector{Int}
    dimension_component_counter::Dict{Int, Dict{Int, Int}}
    index_tracker::Dict{Int, Dict{Int, Vector{Int}}}
    homogenization_matrix_::Matrix{Int}
    num_variables::Int

    function NumericalIrreducibleDecomposition()
        return new(WitnessPointMetadata[], WitnessLinearMetadata[], WitnessPatchMetadata[],
                   Int[], Dict{Int, Dict{Int, Int}}(), Dict{Int, Dict{Int, Vector{Int}}}(),
                   Matrix{Int}[], 0)
    end

    function reset!(nid::NumericalIrreducibleDecomposition)
        reset_points(nid)
        reset_linears(nid)
        reset_patches(nid)
        nid.point_metadata = WitnessPointMetadata[]
        nid.linear_metadata = WitnessLinearMetadata[]
        nid.patch_metadata = WitnessPatchMetadata[]
    end

    function populate(nid::NumericalIrreducibleDecomposition, T::tracker_config_t)
        # implementar según sea necesario
    end

    function choose(nid::NumericalIrreducibleDecomposition, options::BertiniRealConfig)
        # implementar según sea necesario
    end

    function best_possible_automatic_set(nid::NumericalIrreducibleDecomposition, options::BertiniRealConfig)
        # implementar según sea necesario
    end

    function choose_set_interactive(nid::NumericalIrreducibleDecomposition, options::BertiniRealConfig)
        # implementar según sea necesario
    end

    function DisplayComponentsOfDim(nid::NumericalIrreducibleDecomposition, dim::Int)
        # implementar según sea necesario
    end

    function form_specific_witness_set(nid::NumericalIrreducibleDecomposition, dim::Int, comp::Int)
        # implementar según sea necesario
    end

    function print(nid::NumericalIrreducibleDecomposition)
        # implementar según sea necesario
    end

    function add_linear_w_meta(nid::NumericalIrreducibleDecomposition, lin::Vector{Float64}, meta::WitnessLinearMetadata)
        # implementar según sea necesario
    end

    function add_patch_w_meta(nid::NumericalIrreducibleDecomposition, pat::Vector{Float64}, meta::WitnessPatchMetadata)
        # implementar según sea necesario
    end

    function add_solution(nid::NumericalIrreducibleDecomposition, pt::Vector{Float64}, meta::WitnessPointMetadata)::Int
        # implementar según sea necesario
    end
end
