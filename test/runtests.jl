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
    # pattern): only its Taubin-heart testsets are gated internally.
    #
    # No ordering constraint here (2026-07): this file's Taubin verdict
    # testset used to permanently monkeypatch the real
    # resolve_isosingular_dimension, which required it to stay the LAST
    # include in this testset. It now uses that function's own private
    # _resolve_isosingular_trace ScopedValue hook instead (see
    # src/Solver.jl's docstring and this file's own comment) -- task-local,
    # restores automatically when its `with` block exits, cannot leak into
    # any other test regardless of include order.
    include("test_isosingular_deflation.jl")
end
