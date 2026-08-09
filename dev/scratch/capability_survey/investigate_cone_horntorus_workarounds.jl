# Exploratory, scratch-only: test two config-only workaround hypotheses
# for why cone/horn_torus return a completely empty decomposition.
# No src/ changes. 4 combinations: {cone, horn_torus} x {asymmetric bbox_z,
# projection=:random}.

using HomotopyContinuation
using HomotopyGetsReal
using GeometryBasics
using CairoMakie
using JSON
using Statistics
using Random
using HomotopyContinuation: evaluate

const SURVEY_DIR = "/Users/juancagc/HomotopyGetsReal/dev/scratch/capability_survey"
const DATA_DIR = joinpath(SURVEY_DIR, "data")
const RENDER_DIR = joinpath(SURVEY_DIR, "renders")

residual_at(F, point) = abs(only(evaluate(F.expressions, F.variables => point)))

function stats_block(residuals::Vector{Float64})
    isempty(residuals) && return Dict{String,Any}(
        "mean" => nothing, "median" => nothing, "p90" => nothing,
        "p99" => nothing, "max" => nothing, "n_points" => 0,
    )
    return Dict{String,Any}(
        "mean" => Statistics.mean(residuals), "median" => Statistics.median(residuals),
        "p90" => Statistics.quantile(residuals, 0.9), "p99" => Statistics.quantile(residuals, 0.99),
        "max" => maximum(residuals), "n_points" => length(residuals),
    )
end

function vertex_counts(vertices)
    return Dict{String,Any}(
        "Critical" => count(v -> v.v_type == Critical, vertices),
        "Boundary" => count(v -> v.v_type == Boundary, vertices),
        "Singular" => count(v -> v.v_type == Singular, vertices),
        "Artificial" => count(v -> v.v_type == Artificial, vertices),
    )
end

function write_json(path, data)
    open(path, "w") do io
        JSON.print(io, data, 2)
    end
end

@var x y z
fixtures = Dict(
    "cone" => (System([x^2 + y^2 - z^2], variables = [x, y, z]), "x^2 + y^2 - z^2"),
    "horn_torus" => (System([(x^2 + y^2 + z^2)^2 - 4 * (x^2 + y^2)], variables = [x, y, z]), "(x^2+y^2+z^2)^2 - 4*(x^2+y^2)"),
)

function try_decompose(label, F, eqn_str, cfg; extra_kwargs = NamedTuple())
    println("\n=== $label ===")
    outcome = "clean_success"
    error_str = nothing
    t0 = time()
    local vertices, edges, faces, mesh
    try
        if haskey(extra_kwargs, :projection)
            vertices, edges, faces, mesh = decompose_3d_surface(F, cfg; projection = extra_kwargs[:projection], rng = extra_kwargs[:rng])
        else
            vertices, edges, faces, mesh = decompose_3d_surface(F, cfg)
        end
    catch e
        outcome = "exception"
        io = IOBuffer()
        showerror(io, e, catch_backtrace())
        error_str = String(take!(io))
        println("  EXCEPTION: ", first(split(error_str, '\n')))
        wall = time() - t0
        result = Dict{String,Any}(
            "fixture_name" => label, "category" => "surface", "equation" => eqn_str,
            "config" => Dict("edge_sample_density" => cfg.edge_sample_density, "midslice_sample_density" => cfg.midslice_sample_density,
                              "bbox_z" => collect(cfg.bbox_z), "projection" => haskey(extra_kwargs, :projection) ? string(extra_kwargs[:projection]) : nothing),
            "wall_time_seconds" => wall, "outcome" => outcome,
            "vertex_counts" => Dict("Critical" => nothing, "Boundary" => nothing, "Singular" => nothing, "Artificial" => nothing),
            "edge_count" => nothing, "face_count" => nothing, "mesh_vertex_count" => nothing, "mesh_triangle_count" => nothing,
            "residuals" => stats_block(Float64[]), "error" => error_str, "notes" => "",
        )
        write_json(joinpath(DATA_DIR, "$(label).json"), result)
        return (:exception, nothing)
    end
    wall = time() - t0
    vcounts = vertex_counts(vertices)
    is_empty = isempty(vertices) && isempty(faces)
    println("  wall=$(round(wall,digits=2))s  vertices=$(length(vertices)) $vcounts  edges=$(length(edges))  faces=$(length(faces))  mesh_pts=$(length(GeometryBasics.coordinates(mesh)))  mesh_tris=$(length(GeometryBasics.faces(mesh)))")

    resids = Float64[]
    if !is_empty
        pts = GeometryBasics.coordinates(mesh)
        resids = Float64[residual_at(F, Float64.(p)) for p in pts]
        println("  residuals: ", stats_block(resids))
    end
    outcome = is_empty ? "empty_real_locus" : "clean_success"

    result = Dict{String,Any}(
        "fixture_name" => label, "category" => "surface", "equation" => eqn_str,
        "config" => Dict("edge_sample_density" => cfg.edge_sample_density, "midslice_sample_density" => cfg.midslice_sample_density,
                          "bbox_z" => collect(cfg.bbox_z), "projection" => haskey(extra_kwargs, :projection) ? string(extra_kwargs[:projection]) : nothing),
        "wall_time_seconds" => wall, "outcome" => outcome,
        "vertex_counts" => vcounts, "edge_count" => length(edges), "face_count" => length(faces),
        "mesh_vertex_count" => length(GeometryBasics.coordinates(mesh)), "mesh_triangle_count" => length(GeometryBasics.faces(mesh)),
        "residuals" => stats_block(resids), "error" => nothing, "notes" => "",
    )
    write_json(joinpath(DATA_DIR, "$(label).json"), result)

    if !is_empty
        try
            fig = plot_surface_decomposition(mesh; color_by = :z, cfg = cfg)
            png_path = joinpath(RENDER_DIR, "$(label).png")
            CairoMakie.save(png_path, fig)
            println("  Saved render: $png_path")
        catch pe
            println("  PLOT FAILED: ", first(split(sprint(showerror, pe), '\n')))
        end
    end

    return (is_empty ? :still_empty : :worked, mesh)
end

results = Dict{String,String}()

# --- Hypothesis A: asymmetric bbox_z = (-4.0, 4.3), coarse density, default bbox_x/y ---
cfg_A = HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8, bbox_z = (-4.0, 4.3))
for (name, (F, eqn)) in fixtures
    label = "$(name)_asymmetric_bbox"
    r, _ = try_decompose(label, F, eqn, cfg_A)
    results["A_$name"] = string(r)
end

# --- Hypothesis B: projection=:random, rng=Xoshiro(42), symmetric default bbox_z=(-4,4) ---
cfg_B = HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8)
for (name, (F, eqn)) in fixtures
    label = "$(name)_projection_random"
    r, _ = try_decompose(label, F, eqn, cfg_B; extra_kwargs = (projection = :random, rng = Xoshiro(42)))
    results["B_$name"] = string(r)
end

println("\n=== SUMMARY ===")
for k in ["A_cone", "A_horn_torus", "B_cone", "B_horn_torus"]
    println("  $k => $(results[k])")
end
println("ALL DONE")
