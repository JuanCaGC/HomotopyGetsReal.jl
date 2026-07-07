# scratch_phase1_check.jl
#
# Standalone Phase 1 sanity check (NOT part of test/ yet).
# Verifies HomotopyConfig{T}, NativeVertex{T}, Edge{T}, Face{T} all
# instantiate correctly for both Float64 and BigFloat, and that the
# outer constructors are type-stable via @inferred.
#
# Run with:
#   julia --project=. scratch_phase1_check.jl

using Test

include(joinpath(@__DIR__, "src", "Config.jl"))
include(joinpath(@__DIR__, "src", "Types.jl"))

println("=" ^ 60)
println("1. HomotopyConfig{T} instantiation")
println("=" ^ 60)

cfg64 = HomotopyConfig{Float64}()
cfgbig = HomotopyConfig{BigFloat}()

println("HomotopyConfig{Float64}():")
println("  critical_point_tol   = ", cfg64.critical_point_tol, " :: ", typeof(cfg64.critical_point_tol))
println("  boundary_point_tol   = ", cfg64.boundary_point_tol, " :: ", typeof(cfg64.boundary_point_tol))
println("  vertex_match_tol     = ", cfg64.vertex_match_tol, " :: ", typeof(cfg64.vertex_match_tol))
println("  bbox_x               = ", cfg64.bbox_x, " :: ", typeof(cfg64.bbox_x))
println("  max_path_steps       = ", cfg64.max_path_steps, " :: ", typeof(cfg64.max_path_steps))

println("\nHomotopyConfig{BigFloat}():")
println("  critical_point_tol   = ", cfgbig.critical_point_tol, " :: ", typeof(cfgbig.critical_point_tol))
println("  boundary_point_tol   = ", cfgbig.boundary_point_tol, " :: ", typeof(cfgbig.boundary_point_tol))
println("  vertex_match_tol     = ", cfgbig.vertex_match_tol, " :: ", typeof(cfgbig.vertex_match_tol))
println("  bbox_x               = ", cfgbig.bbox_x, " :: ", typeof(cfgbig.bbox_x))

@test cfg64 isa HomotopyConfig{Float64}
@test cfgbig isa HomotopyConfig{BigFloat}
@test cfg64.critical_point_tol isa Float64
@test cfgbig.critical_point_tol isa BigFloat
@test typeof(cfg64.bbox_x) == Tuple{Float64,Float64}
@test typeof(cfgbig.bbox_x) == Tuple{BigFloat,BigFloat}

# Custom override still respects T.
cfg_custom = HomotopyConfig{Float64}(critical_point_tol = 1e-3, max_path_steps = 500)
@test cfg_custom.critical_point_tol == 1e-3
@test cfg_custom.max_path_steps == 500

println("\nHomotopyConfig{T} checks passed.")

println()
println("=" ^ 60)
println("2. NativeVertex{Float64}, Edge{Float64}, Face{Float64}")
println("=" ^ 60)

v = NativeVertex{Float64}(
    id = 1,
    coordinates = [1.0 + 0.0im, 2.0 + 3.0im, 0.0 + 0.0im],
    v_type = Critical,
    metadata = Dict{Symbol,Any}(:jacobian_rank => 2, :tolerance_used => 1e-6),
)
println("v = ", v)
println("  typeof(v.coordinates) = ", typeof(v.coordinates))
println("  typeof(v.metadata)    = ", typeof(v.metadata))

e = Edge{Float64}(
    id = 1,
    left_vertex_id = 1,
    right_vertex_id = 2,
    sampled_points = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [1.0, 1.0, 0.0]],
    is_singular = false,
)
println("e = ", e)
println("  typeof(e.sampled_points)    = ", typeof(e.sampled_points))
println("  typeof(e.sampled_points[1]) = ", typeof(e.sampled_points[1]))

f = Face{Float64}(
    id = 1,
    mid_slice_z = 0.0,
    boundary_edges = [1, 2, 3],
    mesh_vertices = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0],
    mesh_topology = [1 2 3],
)
println("f = ", f)
println("  typeof(f.mesh_vertices)  = ", typeof(f.mesh_vertices))
println("  typeof(f.mesh_topology)  = ", typeof(f.mesh_topology))

@test v isa NativeVertex{Float64}
@test v.coordinates isa Vector{ComplexF64}
@test v.metadata isa Dict{Symbol,Any}

@test e isa Edge{Float64}
@test e.sampled_points isa Vector{Vector{Float64}}
@test eltype(e.sampled_points) == Vector{Float64}

@test f isa Face{Float64}
@test f.mesh_vertices isa Matrix{Float64}
@test f.mesh_topology isa Matrix{Int}

println("\nNativeVertex/Edge/Face{Float64} checks passed.")

println()
println("=" ^ 60)
println("3. @inferred type-stability checks")
println("=" ^ 60)

infer_vertex() = NativeVertex{Float64}(
    id = 1,
    coordinates = [1.0 + 0.0im, 2.0 + 3.0im],
    v_type = Critical,
)
infer_edge() = Edge{Float64}(
    id = 1,
    left_vertex_id = 1,
    right_vertex_id = 2,
    sampled_points = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
)
infer_face() = Face{Float64}(
    id = 1,
    mid_slice_z = 0.0,
    boundary_edges = [1, 2],
    mesh_vertices = [0.0 0.0 0.0; 1.0 0.0 0.0],
    mesh_topology = [1 2],
)
infer_config() = HomotopyConfig{Float64}()
infer_config_big() = HomotopyConfig{BigFloat}()

v2 = @inferred infer_vertex()
e2 = @inferred infer_edge()
f2 = @inferred infer_face()
c2 = @inferred infer_config()
c3 = @inferred infer_config_big()

println("@inferred NativeVertex{Float64} constructor -> ", typeof(v2), "  OK")
println("@inferred Edge{Float64} constructor         -> ", typeof(e2), "  OK")
println("@inferred Face{Float64} constructor         -> ", typeof(f2), "  OK")
println("@inferred HomotopyConfig{Float64}()         -> ", typeof(c2), "  OK")
println("@inferred HomotopyConfig{BigFloat}()        -> ", typeof(c3), "  OK")

println()
println("=" ^ 60)
println("4. Field-type red-flag scan (no Any leaking except metadata)")
println("=" ^ 60)

function scan_for_any(::Type{S}) where {S}
    flags = String[]
    for (fname, ftype) in zip(fieldnames(S), fieldtypes(S))
        if ftype === Any || ftype === Vector{Any} || (ftype isa UnionAll)
            if fname !== :metadata
                push!(flags, "$(S).$(fname) :: $(ftype)")
            end
        end
    end
    return flags
end

all_flags = String[]
append!(all_flags, scan_for_any(typeof(v)))
append!(all_flags, scan_for_any(typeof(e)))
append!(all_flags, scan_for_any(typeof(f)))
append!(all_flags, scan_for_any(typeof(cfg64)))

if isempty(all_flags)
    println("No red flags found (metadata::Dict{Symbol,Any} is the only intentional exception).")
else
    println("RED FLAGS FOUND:")
    foreach(println, all_flags)
end
@test isempty(all_flags)

println()
println("All Phase 1 sanity checks PASSED.")
