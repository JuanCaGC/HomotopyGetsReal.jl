# run_fixture.jl <fixture_name>
#
# One-shot driver: builds the System for <fixture_name>, runs
# decompose_1d_curve / decompose_3d_surface, computes residual statistics,
# writes dev/scratch/capability_survey/data/<fixture_name>.json and (if it
# got far enough) renders/<fixture_name>.png, then exits. Intended to be
# invoked as its own `julia --project=... run_fixture.jl <name>` subprocess
# per fixture (see capability_survey task spec) so a crash/hang in one
# fixture cannot affect any other.

using HomotopyContinuation
using HomotopyGetsReal
using GeometryBasics
using CairoMakie
using JSON
using Statistics
using LinearAlgebra
using Random
using HomotopyContinuation: evaluate

const SURVEY_DIR = "/Users/juancagc/HomotopyGetsReal/dev/scratch/capability_survey"
const DATA_DIR = joinpath(SURVEY_DIR, "data")
const RENDER_DIR = joinpath(SURVEY_DIR, "renders")

residual_at(F, point) = abs(only(evaluate(F.expressions, F.variables => point)))

function write_json(path, data)
    open(path, "w") do io
        JSON.print(io, data, 2)
    end
end

function stats_block(residuals::Vector{Float64})
    if isempty(residuals)
        return Dict{String,Any}(
            "mean" => nothing, "median" => nothing, "p90" => nothing,
            "p99" => nothing, "max" => nothing, "n_points" => 0,
        )
    end
    return Dict{String,Any}(
        "mean" => Statistics.mean(residuals),
        "median" => Statistics.median(residuals),
        "p90" => Statistics.quantile(residuals, 0.9),
        "p99" => Statistics.quantile(residuals, 0.99),
        "max" => maximum(residuals),
        "n_points" => length(residuals),
    )
end

function curve_residual_points(vertices, edges)
    pts = Vector{Vector{Float64}}()
    vmap = Dict(v.id => v for v in vertices)
    for e in edges
        if !isempty(e.sampled_points)
            for p in e.sampled_points
                push!(pts, Float64.(p))
            end
        else
            for vid in (e.left_vertex_id, e.right_vertex_id)
                if haskey(vmap, vid)
                    push!(pts, Float64.(real.(vmap[vid].coordinates)))
                end
            end
        end
    end
    if isempty(pts) && !isempty(vertices)
        for v in vertices
            push!(pts, Float64.(real.(v.coordinates)))
        end
    end
    return pts
end

function surface_residual_points(mesh)
    pts = GeometryBasics.coordinates(mesh)
    return [Float64.(p) for p in pts]
end

function vertex_counts(vertices)
    return Dict{String,Any}(
        "Critical" => count(v -> v.v_type == Critical, vertices),
        "Boundary" => count(v -> v.v_type == Boundary, vertices),
        "Singular" => count(v -> v.v_type == Singular, vertices),
        "Artificial" => count(v -> v.v_type == Artificial, vertices),
    )
end

# ----------------------------------------------------------------------
# Fixture registry: name -> (category, equation_str, F-builder, cfg,
# decompose kwargs NamedTuple, config-json-extra Dict)
# ----------------------------------------------------------------------

function build_fixture(name::String)
    if name == "circle"
        @var x y
        F = System([x^2 + y^2 - 1], variables = [x, y])
        return ("curve", "x^2 + y^2 - 1", F,
                HomotopyConfig{Float64}(edge_sample_density = 6), NamedTuple())
    elseif name == "astroid"
        @var x y
        F = System([(x^2 + y^2 - 1)^3 + 27 * x^2 * y^2], variables = [x, y])
        return ("curve", "(x^2+y^2-1)^3 + 27*x^2*y^2", F,
                HomotopyConfig{Float64}(edge_sample_density = 6), NamedTuple())
    elseif name == "nodal_cubic"
        @var x y
        F = System([y^2 - x^2 * (x + 1)], variables = [x, y])
        return ("curve", "y^2 - x^2*(x+1)", F,
                HomotopyConfig{Float64}(edge_sample_density = 6), NamedTuple())
    elseif name == "cuspidal_cubic"
        @var x y
        F = System([y^2 - x^3], variables = [x, y])
        return ("curve", "y^2 - x^3", F,
                HomotopyConfig{Float64}(edge_sample_density = 6), NamedTuple())
    elseif name == "squircle_quartic"
        @var x y
        F = System([x^4 + y^4 - 1], variables = [x, y])
        return ("curve", "x^4 + y^4 - 1", F,
                HomotopyConfig{Float64}(edge_sample_density = 6), NamedTuple())
    elseif name == "lemniscate_bernoulli"
        @var x y
        F = System([(x^2 + y^2)^2 - x^2 + y^2], variables = [x, y])
        return ("curve", "(x^2+y^2)^2 - x^2 + y^2", F,
                HomotopyConfig{Float64}(edge_sample_density = 6), NamedTuple())
    elseif name == "folium_descartes"
        @var x y
        F = System([x^3 + y^3 - 3 * x * y], variables = [x, y])
        return ("curve", "x^3 + y^3 - 3*x*y", F,
                HomotopyConfig{Float64}(edge_sample_density = 6), NamedTuple())
    elseif name == "cissoid_diocles"
        @var x y
        F = System([y^2 * (2 - x) - x^3], variables = [x, y])
        return ("curve", "y^2*(2-x) - x^3", F,
                HomotopyConfig{Float64}(edge_sample_density = 6), NamedTuple())
    elseif name == "tangent_parabolas_reducible"
        @var x y
        F = System([y^2 - x^4], variables = [x, y])
        return ("curve", "y^2 - x^4  [= (y-x^2)(y+x^2)]", F,
                HomotopyConfig{Float64}(edge_sample_density = 6), NamedTuple())
    elseif name == "three_concurrent_lines_reducible"
        @var x y
        F = System([x * (x - y) * (x + y)], variables = [x, y])
        return ("curve", "x*(x-y)*(x+y)", F,
                HomotopyConfig{Float64}(edge_sample_density = 6), NamedTuple())
    elseif name == "elliptic_two_components"
        @var x y
        F = System([y^2 - x * (x - 1) * (x + 1)], variables = [x, y])
        return ("curve", "y^2 - x*(x-1)*(x+1)", F,
                HomotopyConfig{Float64}(edge_sample_density = 6), NamedTuple())
    elseif name == "empty_curve"
        @var x y
        F = System([x^2 + y^2 + 1], variables = [x, y])
        return ("curve", "x^2 + y^2 + 1", F,
                HomotopyConfig{Float64}(edge_sample_density = 6), NamedTuple())

    elseif name == "sphere"
        @var x y z
        F = System([x^2 + y^2 + z^2 - 1], variables = [x, y, z])
        return ("surface", "x^2 + y^2 + z^2 - 1", F,
                HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8), NamedTuple())
    elseif name == "ellipsoid"
        @var x y z
        F = System([x^2 + 4 * y^2 + 9 * z^2 - 1], variables = [x, y, z])
        return ("surface", "x^2 + 4*y^2 + 9*z^2 - 1", F,
                HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8), NamedTuple())
    elseif name == "taubin_heart"
        @var x y z
        F = System([(x^2 + (1.2 * y)^2 + z^2 - 1)^3 - x^2 * z^3 - 0.1 * (1.2 * y)^2 * z^3], variables = [x, y, z])
        return ("surface", "(x^2+(1.2*y)^2+z^2-1)^3 - x^2*z^3 - 0.1*(1.2*y)^2*z^3", F,
                HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8), NamedTuple())
    elseif name == "torus"
        @var x y z
        F = System([(x^2 + y^2 + z^2 + 3)^2 - 16 * (x^2 + y^2)], variables = [x, y, z])
        return ("surface", "(x^2+y^2+z^2+3)^2 - 16*(x^2+y^2)", F,
                HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8),
                (projection = :random, rng = Xoshiro(42)))
    elseif name == "whitney_umbrella"
        @var x y z
        F = System([x^2 - y^2 * z], variables = [x, y, z])
        return ("surface", "x^2 - y^2*z", F,
                HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8), NamedTuple())
    elseif name == "horn_torus"
        @var x y z
        F = System([(x^2 + y^2 + z^2)^2 - 4 * (x^2 + y^2)], variables = [x, y, z])
        return ("surface", "(x^2+y^2+z^2)^2 - 4*(x^2+y^2)", F,
                HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8), NamedTuple())
    elseif name == "cone"
        @var x y z
        F = System([x^2 + y^2 - z^2], variables = [x, y, z])
        return ("surface", "x^2 + y^2 - z^2", F,
                HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8), NamedTuple())
    elseif name == "hyperboloid_one_sheet"
        @var x y z
        F = System([x^2 + y^2 - z^2 - 1], variables = [x, y, z])
        return ("surface", "x^2 + y^2 - z^2 - 1", F,
                HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8), NamedTuple())
    elseif name == "hyperboloid_two_sheets_reducible"
        @var x y z
        F = System([x^2 + y^2 - z^2 + 1], variables = [x, y, z])
        return ("surface", "x^2 + y^2 - z^2 + 1", F,
                HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8), NamedTuple())
    elseif name == "elliptic_paraboloid"
        @var x y z
        F = System([x^2 + y^2 - z], variables = [x, y, z])
        return ("surface", "x^2 + y^2 - z", F,
                HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8), NamedTuple())
    elseif name == "quartic_superellipsoid"
        @var x y z
        F = System([x^4 + y^4 + z^4 - 1], variables = [x, y, z])
        return ("surface", "x^4 + y^4 + z^4 - 1", F,
                HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8), NamedTuple())
    elseif name == "empty_surface"
        @var x y z
        F = System([x^2 + y^2 + z^2 + 1], variables = [x, y, z])
        return ("surface", "x^2 + y^2 + z^2 + 1", F,
                HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8), NamedTuple())
    else
        error("Unknown fixture name: $name")
    end
end

# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

function main()
    length(ARGS) >= 1 || error("usage: julia run_fixture.jl <fixture_name>")
    fixture_name = ARGS[1]
    println("=== Running fixture: $fixture_name ===")

    category, equation_str, F, cfg, extra_kwargs = build_fixture(fixture_name)

    config_json = Dict{String,Any}("edge_sample_density" => cfg.edge_sample_density)
    if category == "surface"
        config_json["midslice_sample_density"] = cfg.midslice_sample_density
    end
    if haskey(extra_kwargs, :projection)
        config_json["projection"] = string(extra_kwargs[:projection])
    end
    if haskey(extra_kwargs, :rng)
        config_json["rng"] = "Xoshiro(42)"
    end

    vertices = nothing
    edges = nothing
    faces = nothing
    mesh = nothing
    outcome = "clean_success"
    error_str = nothing
    notes = ""

    t0 = time()
    try
        if category == "curve"
            vertices, edges = decompose_1d_curve(F, cfg)
        else
            if haskey(extra_kwargs, :projection)
                vertices, edges, faces, mesh = decompose_3d_surface(
                    F, cfg; projection = extra_kwargs[:projection], rng = extra_kwargs[:rng],
                )
            else
                vertices, edges, faces, mesh = decompose_3d_surface(F, cfg)
            end
        end
    catch e
        outcome = "exception"
        io = IOBuffer()
        showerror(io, e, catch_backtrace())
        error_str = String(take!(io))
        println("EXCEPTION caught:")
        println(error_str)
    end
    wall_time = time() - t0

    vcounts = Dict{String,Any}("Critical" => nothing, "Boundary" => nothing, "Singular" => nothing, "Artificial" => nothing)
    edge_count = nothing
    face_count = nothing
    mesh_vertex_count = nothing
    mesh_triangle_count = nothing
    residuals_block = stats_block(Float64[])

    if outcome != "exception"
        vcounts = vertex_counts(vertices)
        edge_count = length(edges)

        if category == "curve"
            pts = curve_residual_points(vertices, edges)
            resids = Float64[residual_at(F, p) for p in pts]
            residuals_block = stats_block(resids)
            if isempty(vertices) && isempty(edges)
                notes *= "Empty result: 0 vertices, 0 edges. "
            end
        else
            face_count = length(faces)
            mesh_vertex_count = length(GeometryBasics.coordinates(mesh))
            mesh_triangle_count = length(GeometryBasics.faces(mesh))
            pts = surface_residual_points(mesh)
            resids = Float64[residual_at(F, p) for p in pts]
            residuals_block = stats_block(resids)
            if isempty(vertices) && isempty(faces)
                notes *= "Empty result: 0 vertices, 0 faces. "
            end
        end
    end

    # Plotting -- best-effort; failure here does not affect outcome/JSON.
    png_path = joinpath(RENDER_DIR, "$(fixture_name).png")
    if outcome != "exception"
        try
            if category == "curve"
                fig = plot_curve_decomposition(vertices, edges; cfg = cfg)
                CairoMakie.save(png_path, fig)
                println("Saved render: $png_path")
            else
                fig = plot_surface_decomposition(mesh; color_by = :z, cfg = cfg)
                CairoMakie.save(png_path, fig)
                println("Saved render: $png_path")
            end
        catch plot_e
            io = IOBuffer()
            showerror(io, plot_e, catch_backtrace())
            plot_err_str = String(take!(io))
            notes *= "PLOT FAILED: $(plot_err_str[1:min(end,500)]) "
            println("Plotting failed (non-fatal):")
            println(plot_err_str)
        end
    end

    result = Dict{String,Any}(
        "fixture_name" => fixture_name,
        "category" => category,
        "equation" => equation_str,
        "config" => config_json,
        "wall_time_seconds" => wall_time,
        "outcome" => outcome,
        "vertex_counts" => vcounts,
        "edge_count" => edge_count,
        "face_count" => face_count,
        "mesh_vertex_count" => mesh_vertex_count,
        "mesh_triangle_count" => mesh_triangle_count,
        "residuals" => residuals_block,
        "error" => error_str,
        "notes" => notes,
    )

    json_path = joinpath(DATA_DIR, "$(fixture_name).json")
    write_json(json_path, result)
    println("Saved JSON: $json_path")
    println("=== Done: $fixture_name (outcome=$outcome, wall_time=$(round(wall_time, digits=2))s) ===")
end

main()
