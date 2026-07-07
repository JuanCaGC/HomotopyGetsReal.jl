using HomotopyContinuation
using LinearAlgebra
using Statistics
using GeometryBasics
using GLMakie
using HomotopyGetsReal
using Test

const _TEST_OUTPUT = mkpath(joinpath(@__DIR__, "output"))

@testset "HomotopyGetsReal" begin
    include("test_types.jl")
    include("test_solver.jl")
    include("test_topology.jl")
    include("test_pathtracking.jl")
    include("test_surfacedecomposition.jl")
    include("test_visuals.jl")

    if get(ENV, "HOMOTOPYGETSREAL_RUN_SLOW_TESTS", "0") == "1"
        include("test_taubin.jl")
    else
        @info "Skipping test_taubin.jl (set HOMOTOPYGETSREAL_RUN_SLOW_TESTS=1 to run)"
    end
end
