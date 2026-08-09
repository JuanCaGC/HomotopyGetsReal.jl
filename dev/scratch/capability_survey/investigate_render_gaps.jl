# Targeted follow-up: torus/taubin_heart naked-edge hypothesis, and
# hyperboloid_one_sheet density hypothesis. Not part of the original
# 24-fixture sweep -- reuses this directory's satellite environment.

using HomotopyContinuation
using HomotopyGetsReal
using GeometryBasics
using CairoMakie
using Random

const RENDER_DIR = "/Users/juancagc/HomotopyGetsReal/dev/scratch/capability_survey/renders"

naked_count(mesh) = length(HomotopyGetsReal._naked_mesh_edges(mesh))

function report(label, v, e, f, m, wall)
    nk = naked_count(m)
    println("$label: wall=$(round(wall,digits=2))s  vertices=$(length(v))  edges=$(length(e))  faces=$(length(f))  mesh_verts=$(length(GeometryBasics.coordinates(m)))  mesh_tris=$(length(GeometryBasics.faces(m)))  naked_edges=$nk")
    return nk
end

@var x y z
survey_cfg = HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8)

# ---------------------------------------------------------------------
# TORUS: bare vs incidence=true, same config as original survey run
# ---------------------------------------------------------------------
Ftorus = System([(x^2 + y^2 + z^2 + 3)^2 - 16 * (x^2 + y^2)], variables = [x, y, z])

println("\n=== TORUS bare (survey config, projection=:random, rng=Xoshiro(42)) ===")
t0 = time()
v, e, f, m = decompose_3d_surface(Ftorus, survey_cfg; projection = :random, rng = Xoshiro(42))
nk_torus_bare = report("torus_bare", v, e, f, m, time() - t0)
fig = plot_surface_decomposition(m; color_by = :z, cfg = survey_cfg)
CairoMakie.save(joinpath(RENDER_DIR, "torus_bare.png"), fig)

println("\n=== TORUS incidence=true (same config) ===")
t0 = time()
v, e, f, m, inc = decompose_3d_surface(Ftorus, survey_cfg; projection = :random, rng = Xoshiro(42), incidence = true)
nk_torus_inc = report("torus_incidence", v, e, f, m, time() - t0)
fig = plot_surface_decomposition(m; color_by = :z, cfg = survey_cfg)
CairoMakie.save(joinpath(RENDER_DIR, "torus_incidence.png"), fig)

# ---------------------------------------------------------------------
# TAUBIN HEART: bare vs incidence=true, same config as original survey run
# ---------------------------------------------------------------------
Ftaubin = System([(x^2 + (1.2 * y)^2 + z^2 - 1)^3 - x^2 * z^3 - 0.1 * (1.2 * y)^2 * z^3], variables = [x, y, z])

println("\n=== TAUBIN_HEART bare (survey config) ===")
t0 = time()
v, e, f, m = decompose_3d_surface(Ftaubin, survey_cfg)
nk_taubin_bare = report("taubin_heart_bare", v, e, f, m, time() - t0)
fig = plot_surface_decomposition(m; color_by = :z, cfg = survey_cfg)
CairoMakie.save(joinpath(RENDER_DIR, "taubin_heart_bare.png"), fig)

println("\n=== TAUBIN_HEART incidence=true (same config) ===")
t0 = time()
v, e, f, m, inc = decompose_3d_surface(Ftaubin, survey_cfg; incidence = true)
nk_taubin_inc = report("taubin_heart_incidence", v, e, f, m, time() - t0)
fig = plot_surface_decomposition(m; color_by = :z, cfg = survey_cfg)
CairoMakie.save(joinpath(RENDER_DIR, "taubin_heart_incidence.png"), fig)

# ---------------------------------------------------------------------
# HYPERBOLOID ONE SHEET: coarse (survey) density already rendered in the
# original survey (renders/hyperboloid_one_sheet.png) -- only need the
# default-density comparison run here.
# ---------------------------------------------------------------------
Fhyp = System([x^2 + y^2 - z^2 - 1], variables = [x, y, z])
default_cfg = HomotopyConfig{Float64}()  # edge_sample_density=50, midslice_sample_density=100

println("\n=== HYPERBOLOID_ONE_SHEET default density (edge_sample_density=$(default_cfg.edge_sample_density), midslice_sample_density=$(default_cfg.midslice_sample_density)) ===")
t0 = time()
v, e, f, m = decompose_3d_surface(Fhyp, default_cfg)
nk_hyp_default = report("hyperboloid_one_sheet_default_density", v, e, f, m, time() - t0)
fig = plot_surface_decomposition(m; color_by = :z, cfg = default_cfg)
CairoMakie.save(joinpath(RENDER_DIR, "hyperboloid_one_sheet_default_density.png"), fig)

println("\n=== SUMMARY ===")
println("torus:        naked_edges bare=$nk_torus_bare  incidence=$nk_torus_inc")
println("taubin_heart: naked_edges bare=$nk_taubin_bare  incidence=$nk_taubin_inc")
println("hyperboloid_one_sheet: naked_edges @ default density = $nk_hyp_default (survey/coarse density not re-measured here, only visually compared)")
println("ALL DONE")
