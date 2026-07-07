@testset "Types and Config (Phase 1)" begin
    cfg64 = HomotopyConfig{Float64}()
    cfgbig = HomotopyConfig{BigFloat}()

    @test cfg64 isa HomotopyConfig{Float64}
    @test cfgbig isa HomotopyConfig{BigFloat}
    @test cfg64.critical_point_tol isa Float64
    @test cfgbig.critical_point_tol isa BigFloat
    @test typeof(cfg64.bbox_x) == Tuple{Float64,Float64}
    @test typeof(cfgbig.bbox_x) == Tuple{BigFloat,BigFloat}

    cfg_custom = HomotopyConfig{Float64}(critical_point_tol = 1e-3, max_path_steps = 500)
    @test cfg_custom.critical_point_tol == 1e-3
    @test cfg_custom.max_path_steps == 500

    v = NativeVertex{Float64}(
        id = 1,
        coordinates = [1.0 + 0.0im, 2.0 + 3.0im, 0.0 + 0.0im],
        v_type = Critical,
        metadata = Dict{Symbol,Any}(:jacobian_rank => 2, :tolerance_used => 1e-6),
    )
    e = Edge{Float64}(
        id = 1,
        left_vertex_id = 1,
        right_vertex_id = 2,
        sampled_points = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [1.0, 1.0, 0.0]],
        is_singular = false,
    )
    f = Face{Float64}(
        id = 1,
        mid_slice_z = 0.0,
        boundary_edges = [1, 2, 3],
        mesh_vertices = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0],
        mesh_topology = [1 2 3],
    )

    @test v isa NativeVertex{Float64}
    @test v.coordinates isa Vector{ComplexF64}
    @test v.metadata isa Dict{Symbol,Any}
    @test e isa Edge{Float64}
    @test e.sampled_points isa Vector{Vector{Float64}}
    @test eltype(e.sampled_points) == Vector{Float64}
    @test f isa Face{Float64}
    @test f.mesh_vertices isa Matrix{Float64}
    @test f.mesh_topology isa Matrix{Int}

    infer_vertex() = NativeVertex{Float64}(
        id = 1, coordinates = [1.0 + 0.0im, 2.0 + 3.0im], v_type = Critical,
    )
    infer_edge() = Edge{Float64}(
        id = 1, left_vertex_id = 1, right_vertex_id = 2,
        sampled_points = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
    )
    infer_face() = Face{Float64}(
        id = 1, mid_slice_z = 0.0, boundary_edges = [1, 2],
        mesh_vertices = [0.0 0.0 0.0; 1.0 0.0 0.0], mesh_topology = [1 2],
    )
    infer_config() = HomotopyConfig{Float64}()
    infer_config_big() = HomotopyConfig{BigFloat}()

    @inferred infer_vertex()
    @inferred infer_edge()
    @inferred infer_face()
    @inferred infer_config()
    @inferred infer_config_big()

    function scan_for_any(::Type{S}) where {S}
        flags = String[]
        for (fname, ftype) in zip(fieldnames(S), fieldtypes(S))
            if ftype === Any || ftype === Vector{Any} || (ftype isa UnionAll)
                fname !== :metadata && push!(flags, "$(S).$(fname) :: $(ftype)")
            end
        end
        return flags
    end

    all_flags = String[]
    append!(all_flags, scan_for_any(typeof(v)))
    append!(all_flags, scan_for_any(typeof(e)))
    append!(all_flags, scan_for_any(typeof(f)))
    append!(all_flags, scan_for_any(typeof(cfg64)))
    @test isempty(all_flags)
end
