@testset "Visuals (Phase 6)" begin
    using LinearAlgebra

    outdir = isdefined(Main, :_TEST_OUTPUT) ? Main._TEST_OUTPUT : mktempdir()

    @testset "GLMakie headless smoke" begin
        fig = Figure()
        ax = Axis(fig[1, 1])
        scatter!(ax, [0.0], [0.0])
        GLMakie.save(joinpath(outdir, "00_smoke_test.png"), fig)
    end

    @var x y z
    f_sphere = x^2 + y^2 + z^2 - 1
    F_sphere = System([f_sphere], variables = [x, y, z])
    cfg = HomotopyConfig{Float64}(edge_sample_density = 6, midslice_sample_density = 8)

    @testset "plot_curve_decomposition" begin
        vertices2d, edges2d = slice_at_z(F_sphere, 0.0, cfg)
        fig1 = plot_curve_decomposition(vertices2d, edges2d; cfg = cfg, show_labels = true)
        @test fig1 isa Makie.Figure
        GLMakie.save(joinpath(outdir, "01_sphere_equator_curve.png"), fig1)
        fig1b = plot_curve_decomposition(vertices2d, edges2d; edge_color_by = :mono, show_vertices = false)
        @test fig1b isa Makie.Figure
    end

    all_vertices, all_edges, all_faces, mesh = decompose_3d_surface(F_sphere, cfg)

    @testset "plot_surface_decomposition(mesh)" begin
        z_warns, fig2 = Test.collect_test_logs() do
            plot_surface_decomposition(mesh; color_by = :z, show_wireframe = true, cfg = cfg, vertices = all_vertices)
        end
        @test fig2 isa Makie.Figure
        @test isempty([l for l in z_warns if l.level == Base.CoreLogging.Warn])

        radial_fn(px, py, pz) = sqrt(px^2 + py^2 + pz^2)
        radial_warns, fig2b = Test.collect_test_logs() do
            plot_surface_decomposition(mesh; color_by = radial_fn, show_colorbar = true)
        end
        @test fig2b isa Makie.Figure
        radial_warn_logs = [l for l in radial_warns if l.level == Base.CoreLogging.Warn]
        @test length(radial_warn_logs) == 1
        @test occursin("near-constant", radial_warn_logs[1].message)

        radial_warns2, _ = Test.collect_test_logs() do
            plot_surface_decomposition(mesh; color_by = radial_fn)
        end
        @test isempty([l for l in radial_warns2 if l.level == Base.CoreLogging.Warn])

        @test_throws ArgumentError plot_surface_decomposition(mesh; color_by = :bogus)
    end

    @testset "plot_surface_decomposition(faces)" begin
        warn_logs, fig3 = Test.collect_test_logs() do
            plot_surface_decomposition(all_faces; show_wireframe = true, cfg = cfg)
        end
        warn_logs1 = [l for l in warn_logs if l.level == Base.CoreLogging.Warn]
        @test length(warn_logs1) == 1
        @test occursin("winding correction", warn_logs1[1].message)
        @test fig3 isa Makie.Figure

        local_warns2 = Test.collect_test_logs() do
            plot_surface_decomposition(all_faces)
        end
        warn_logs2 = [l for l in local_warns2[1] if l.level == Base.CoreLogging.Warn]
        @test isempty(warn_logs2)
    end

    @testset "interactive_3d_viewer" begin
        fig4 = interactive_3d_viewer(mesh; color_by = :z, show_wireframe = false, cfg = cfg)
        @test fig4 isa Makie.Figure
    end

    @testset "ellipsoid mesh plot" begin
        @var xe ye ze
        f_ell = xe^2 + 4 * ye^2 + 9 * ze^2 - 1
        F_ell = System([f_ell], variables = [xe, ye, ze])
        _, _, _, emesh = decompose_3d_surface(F_ell, cfg)
        fig5 = plot_surface_decomposition(emesh; color_by = :z, cfg = cfg)
        @test fig5 isa Makie.Figure
    end

  if get(ENV, "HOMOTOPYGETSREAL_RUN_SLOW_TESTS", "0") == "1"
    @testset "Taubin heart visuals (slow)" begin
        @var xh yh zh
        f_heart = (xh^2 + (1.2 * yh)^2 + zh^2 - 1)^3 - xh^2 * zh^3 - 0.1 * (1.2 * yh)^2 * zh^3
        F_heart = System([f_heart], variables = [xh, yh, zh])
        cfg_heart = HomotopyConfig{Float64}(
            bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-1.3, 1.3),
            edge_sample_density = 8, midslice_sample_density = 8,
        )
        z_crits_heart = compute_critical_z_slices(F_heart, cfg_heart)
        sort!(z_crits_heart)
        narrow_i = argmin(diff(z_crits_heart))
        z_narrow_mid = (z_crits_heart[narrow_i] + z_crits_heart[narrow_i + 1]) / 2
        v_narrow, e_narrow = slice_at_z(F_heart, z_narrow_mid, cfg_heart)
        fig6 = plot_curve_decomposition(v_narrow, e_narrow; cfg = cfg_heart, edge_color_by = :cell)
        @test fig6 isa Makie.Figure
        _, _, hf, hmesh = decompose_3d_surface(F_heart, cfg_heart)
        fig6b = plot_surface_decomposition(hmesh; color_by = :z, show_wireframe = false, cfg = cfg_heart)
        @test fig6b isa Makie.Figure
        fig6c = plot_surface_decomposition(hf; cfg = cfg_heart)
        @test fig6c isa Makie.Figure
    end
  end

    @testset "investigation regressions" begin
        pts = GeometryBasics.coordinates(mesh)
        radial_fn(px, py, pz) = sqrt(px^2 + py^2 + pz^2)
        radvals = [radial_fn(p[1], p[2], p[3]) for p in pts]
        lo, hi = extrema(radvals)
        scale = max(abs(lo), abs(hi), 1.0)
        @test (hi - lo) / scale < 1e-4  # N1: near-constant radial_fn range

        mesh_pts = GeometryBasics.coordinates(mesh)
        mesh_tris = GeometryBasics.faces(mesh)
        n_inward = 0
        for t in mesh_tris
            i1, i2, i3 = Int(t[1]), Int(t[2]), Int(t[3])
            p1, p2, p3 = mesh_pts[i1], mesh_pts[i2], mesh_pts[i3]
            n = cross(Vector(p2 .- p1), Vector(p3 .- p1))
            nnorm = norm(n)
            nnorm < 1e-12 && continue
            n_hat = n ./ nnorm
            centroid = (Vector(p1) .+ Vector(p2) .+ Vector(p3)) ./ 3
            radial_hat = centroid ./ norm(centroid)
            dot(n_hat, radial_hat) <= 0 && (n_inward += 1)
        end
        @test n_inward == 0  # N2: all welded triangles outward-facing

        ms = methods(interactive_3d_viewer)
        @test length(ms) == 1  # N3: single dispatch target
        @test first(ms).sig.parameters[2] == GeometryBasics.Mesh

        fig_direct = plot_surface_decomposition(mesh; color_by = :z, show_wireframe = false, cfg = cfg)
        fig_viewer = interactive_3d_viewer(mesh; color_by = :z, show_wireframe = false, cfg = cfg)
        @test fig_direct isa Makie.Figure && fig_viewer isa Makie.Figure  # N4

        fixed = HomotopyGetsReal._near_constant_colorrange(fill(1.0, 10))
        @test fixed !== nothing  # N5
    end
end
