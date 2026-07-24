using HomotopyContinuation
using LinearAlgebra
using Statistics
using GeometryBasics
using GLMakie
using HomotopyGetsReal
using Test

const _TEST_OUTPUT = mkpath(joinpath(@__DIR__, "output"))

@testset "HomotopyGetsReal" begin
    include("test_docstring_rendering.jl")
    include("test_types.jl")
    include("test_vertex_registry.jl")
    include("test_solver.jl")
    include("test_topology.jl")
    include("test_pathtracking.jl")
    include("test_surfacedecomposition.jl")
    include("test_projection.jl")
    include("test_incidence.jl")
    include("test_visuals.jl")

    if get(ENV, "HOMOTOPYGETSREAL_RUN_SLOW_TESTS", "0") == "1"
        include("test_taubin.jl")
    else
        @info "Skipping test_taubin.jl (set HOMOTOPYGETSREAL_RUN_SLOW_TESTS=1 to run)"
    end

    # Always included (matches test_visuals.jl's own fast/slow-within-one-file
    # pattern): only its Taubin-heart testsets are gated internally. Included
    # AFTER test_taubin.jl -- see test_isosingular_deflation.jl's own comment
    # on why its resolve_isosingular_dimension monkeypatch must run last.
    include("test_isosingular_deflation.jl")
end
