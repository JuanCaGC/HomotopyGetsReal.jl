function convert_test_file(path, testset_name)
    lines = readlines(path)
    out = String[]
    i = 1
    while i <= length(lines) && (startswith(lines[i], "#") || isempty(strip(lines[i])))
        i += 1
    end
    while i <= length(lines)
        l = strip(lines[i])
        skip = startswith(l, "using Test") || startswith(l, "using Pkg") ||
               startswith(l, "using HomotopyContinuation") || startswith(l, "using LinearAlgebra") ||
               startswith(l, "using Statistics") || startswith(l, "using GeometryBasics") ||
               startswith(l, "using GLMakie") || startswith(l, "include(joinpath") ||
               startswith(l, "using .HomotopyGetsReal") || startswith(l, "Pkg.activate")
        if skip
            i += 1
            continue
        end
        break
    end
    push!(out, "@testset \"$testset_name\" begin")
    while i <= length(lines)
        push!(out, lines[i])
        i += 1
    end
    while !isempty(out)
        l = strip(out[end])
        if startswith(l, "println") || occursin("All Phase", l)
            pop!(out)
            continue
        end
        break
    end
    push!(out, "end")
    write(path, join(out, "\n") * "\n")
end

for (f, name) in [
    ("test/test_solver.jl", "Solver (Phase 2)"),
    ("test/test_topology.jl", "Topology (Phase 3)"),
    ("test/test_pathtracking.jl", "PathTracking (Phase 4)"),
    ("test/test_surfacedecomposition.jl", "SurfaceDecomposition (Phase 5)"),
    ("test/test_taubin.jl", "Taubin heart (Phase 5 slow)"),
]
    convert_test_file(f, name)
end
