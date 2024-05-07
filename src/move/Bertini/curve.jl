using Base.Filesystem

function main(V::VertexSet, W_curve::WitnessSet, projections::Array{MPVec}, program_options::BertiniRealConfig, solve_options::SolverConfiguration)
    # Isosingular Deflation
    temp_path = W_curve.input_filename()
    dim_str = "_dim_$(W_curve.dimension())_comp_$(W_curve.component_number())_deflated"
    temp_path *= dim_str
    zeroonly = Set([0])
    # TODO write_dehomogenized_coordinates
    write_dehomogenized_coordinates(W_curve, "witness_points_dehomogenized", zeroonly)
    # TODO isosingular_deflation
    num_deflations, deflation_sequence = isosingular_deflation(program_options, W_curve.input_filename(), "witness_points_dehomogenized", temp_path, program_options.max_deflations())
    # TODO free
    free(deflation_sequence)
    # TODO set_input_deflated_filename
    program_options.set_input_deflated_filename(temp_path)
    # TODO set_input_filename
    W_curve.set_input_filename(temp_path)
    # TODO Decomposition::copy_data_from_witness_set
    copy_data_from_witness_set(W_curve)

    # TODO parse_input_file
    parse_input_file(W_curve.input_filename())
    # TODO preproc_data_clear
    preproc_data_clear(solve_options.PPD)
    # TODO parse_preproc_data
    parse_preproc_data("preproc_data", solve_options.PPD)

    self_conjugate = true
    if W_curve.num_synth_vars() == 0
        println("checking if component is self-conjugate")
        # TODO checkSelfConjugate
        self_conjugate = checkSelfConjugate(W_curve.point(1), program_options, program_options.input_filename())
        # TODO verify_projection_ok
        if verify_projection_ok(W_curve, projections, solve_options) == 1
            println("verified projection is ok")
        else
            println("the projection is invalid, in that the jacobian of the randomized system\nbecomes singular at a random point, when the projection is concatenated\n")
            br_exit(196)
        end
    end
    # TODO user_sphere
    if program_options.user_sphere()
        # TODO read_sphere
        read_sphere(program_options.bounding_sphere_filename())
    end
    # TODO add_projection
    add_projection(projections[1])

    if !self_conjugate
        # TODO computeCurveNotSelfConj
        computeCurveNotSelfConj(W_curve, V, num_variables(), program_options, solve_options)
    else
        # TODO computeCurveSelfConj
        computeCurveSelfConj(W_curve, projections, V, program_options, solve_options)
    end
end

function computeCurveSelfConj(W_curve::WitnessSet, projections::Array{MPVec}, V::VertexSet, program_options::BertiniRealConfig, solve_options::SolverConfiguration)
    println("computeCurveSelfConj")

    # Randomize down to N-1 equations to get a square system for the homotopies
    randomizer().setup(W_curve.num_variables() - W_curve.num_patches() - 1, solve_options.PPD.num_funcs)

    # Solve for critical conditions for random complex projection
    W_crit_real, W_singular = compute_critical_points(W_curve, projections, program_options, solve_options)
    interslice(W_curve, W_crit_real, projections, program_options, solve_options, V)

    if !IsEmbedded()
        V.add_type_to_points(W_singular, Singular)
    end
end

function compute_critical_points(W_curve::WitnessSet, projections::Array{MPVec}, program_options::BertiniRealConfig, solve_options::SolverConfiguration, W_crit_real::WitnessSet, W_singular::WitnessSet)
    println("compute_critical_points")

    if !randomizer().is_ready()
        throw(logic_error("randomizer is not setup at compute_critical_points"))
    end

    W_crit_real.set_input_filename(W_curve.input_filename())

    solve_out = SolverOutput()
    ns_config = NullspaceConfiguration()
    compute_crit_nullspace(solve_out, W_curve, randomizer(), projections, 1, 1, 1, program_options, solve_options, ns_config)
    ns_config.clear()

    solve_out.get_noninfinite_w_mult_full(W_crit_real)
    solve_out.get_sing(W_singular)

    W_crit_real.only_first_vars(W_curve.num_variables())
    W_crit_real.sort_for_real(solve_options.T.real_threshold)
    W_crit_real.sort_for_unique(program_options.same_point_tol())

    W_singular.copy_patches(W_crit_real)
    W_singular.only_first_vars(W_curve.num_variables())
    W_singular.sort_for_real(solve_options.T.real_threshold)
    W_singular.sort_for_unique(program_options.same_point_tol())

    if program_options.verbose_level() >= 2
        println("the critical points of the curve:\n\n")
        W_crit_real.print_to_screen()
    end

    if have_sphere()
        W_crit_real.sort_for_inside_sphere(sphere_radius(), sphere_center())
        W_singular.sort_for_inside_sphere(sphere_radius(), sphere_center())
    else
        println("computing sphere bounds...")
        compute_sphere_bounds(W_crit_real)
    end

    W_sphere_isect = WitnessSet()
    get_sphere_intersection_pts(W_sphere_isect, W_curve, program_options, solve_options)
    W_sphere_isect.sort_for_real(solve_options.T.real_threshold)
    W_sphere_isect.sort_for_unique(program_options.same_point_tol())

    if program_options.verbose_level() >= 2
        println("the sphere intersection points of the curve:\n\n")
        W_sphere_isect.print_to_screen()
    end

    W_crit_real.merge(W_sphere_isect, program_options.same_point_tol())

    return SUCCESSFUL
end

function get_sphere_intersection_pts(W_additional::WitnessSet, W_curve::WitnessSet, program_options::BertiniRealConfig, solve_options::SolverConfiguration)
    println("get_sphere_intersection_pts")

    if !randomizer().is_ready()
        throw(logic_error("randomizer is not ready to go at get_sphere_intersection_pts"))
    end

    if W_curve.num_linears() != 1
        error("the input witness set to get_additional_critpts had an incorrect number of linears: $(W_curve.num_linears())")
    end

    # Build up the start system
    parse_input_file(W_curve.input_filename())
    preproc_data_clear(solve_options.PPD)
    parse_preproc_data("preproc_data", solve_options.PPD)

    ml_config = MultilinConfiguration(solve_options, randomizer())

    multilin_linears = Array{MPVec}(undef, 1)
    init_vec_mp2(multilin_linears[1], W_curve.num_variables(), solve_options.T.AMP_max_prec)
    multilin_linears[1].size = W_curve.num_variables()

    W_sphere = deepcopy(W_curve)
    W_sphere.reset_points()
    W_sphere.reset_linears()
    W_sphere.reset_patches()

    sp_config = SphereConfiguration(randomizer())
    for jj in 1:W_curve.num_variables()
        set_zero_mp(multilin_linears[1].coord[jj])
    end

    for ii in 1:2
        for jj in 1:W_curve.num_natural_variables()
            get_comp_rand_mp(multilin_linears[1].coord[jj])
        end
        vec_cp_mp(sp_config.starting_linear[ii], multilin_linears[1])

        W_temp = WitnessSet()

        fillme = SolverOutput()
        multilin_solver_master_entry_point(W_curve, fillme, multilin_linears, ml_config, solve_options)
        fillme.get_noninfinite_w_mult(W_temp)

        merge(W_sphere, W_temp, program_options.same_point_tol())
    end

    clear_vec_mp(multilin_linears[1])
    free(multilin_linears)

    W_sphere.copy_patches(W_curve)

    if program_options.verbose_level() >= 1
        println("sphere intersection computation")
    end

    sp_config.set_memory(solve_options)
    sp_config.set_center(sphere_center())
    sp_config.set_radius(sphere_radius())

    fillme = SolverOutput()
    sphere_solver_master_entry_point(W_sphere, fillme, sp_config, solve_options)
    fillme.get_noninfinite_w_mult_full(W_additional)

    return SUCCESSFUL
end

function interslice(W_curve::WitnessSet, W_crit_real::WitnessSet, projections::Array{MPVec}, program_options::BertiniRealConfig, solve_options::SolverConfiguration, V::VertexSet)
    println("interslice")

    if !randomizer().is_ready()
        throw(logic_error("in interslice, randomizer is not set up properly."))
    end

    V.set_curr_projection(projections[1])
    V.set_curr_input(W_crit_real.input_filename())

    set_W(W_curve)
    copy_patches(W_curve)
    set_num_variables(W_crit_real.num_variables())
    set_input_filename(W_curve.input_filename())

    blabla = Cint(0)
    parse_input_file(W_curve.input_filename(), &blabla)
    solve_options.get_PPD()

    temp_vertex = Vertex()

    crit_point_counter = Dict{Int, Int}()

    W_canonicalized = WitnessSet()
    num_to_start = V.num_vertices()

    for ii in 1:W_crit_real.num_points()
        temp_vertex.set_point(W_crit_real.point(ii))
        temp_vertex.set_type(Critical)

        I = index_in_vertices_with_add(V, temp_vertex)
        crit_point_counter[I] = 0

        if program_options.verbose_level() >= 8
            println("using point $ii of $(W_crit_real.num_points()) from W_crit_real in VertexSet as point $I")
        end

        if I >= num_to_start
            add_point(W_canonicalized, W_crit_real.point(ii))
        else
            add_point(W_canonicalized, V.GetVertex(I).point())
        end
    end

    crit_downstairs = MPVec(0)
    mid_downstairs = MPVec(0)

    compute_downstairs_crit_midpts(V, W_canonicalized, crit_downstairs, mid_downstairs, projections[1], solve_options.T)

    SetCritSliceValues(crit_downstairs)

    if program_options.verbose_level() >= 0
        print_point_to_screen_matlab(crit_downstairs, "curve_interslice_crit_downstairs")
    end

    V.set_curr_input(W_curve.input_filename())

    num_midpoints = crit_downstairs.size

    edge_counter = 0
    midpoint_witness_sets = Vector{WitnessSet}(undef, num_midpoints)

    ml_config = MultilinConfiguration(solve_options, randomizer())
    particular_projection = MPVec(W_curve.num_variables())
    particular_projection.size = W_curve.num_variables()
    vec_cp_mp(particular_projection, projections[1])

    MidSlice(edge_counter, midpoint_witness_sets, ml_config, W_curve, particular_projection, mid_downstairs, program_options, solve_options)

    found_indices_crit = Vector{Set{Int}}(undef, num_midpoints)
    found_indices_mid = Vector{Set{Int}}(undef, num_midpoints)
    found_indices_left = Set{Int}()
    found_indices_right = Set{Int}()

    ConnectTheDots(found_indices_crit, found_indices_mid, found_indices_left, found_indices_right, crit_point_counter, V, crit_downstairs, mid_downstairs, particular_projection, midpoint_witness_sets, ml_config, program_options, solve_options)

    for ii in 1:num_midpoints
        bad_crit = assert_projection_value(V, found_indices_crit[ii], crit_downstairs.coord[ii])
        bad_mid = assert_projection_value(V, found_indices_mid[ii], mid_downstairs.coord[ii])
    end
    bad_crit = assert_projection_value(V, found_indices_crit[num_midpoints], crit_downstairs.coord[num_midpoints])

    crit_pt_iterator = crit_point_counter
    for (curr_index, _) in crit_pt_iterator
        AddEdge(Edge(curr_index, curr_index, curr_index), EdgeMetaData(0, 0))
    end

    if program_options.merge_edges
        Merge(midpoint_witness_sets[1], V, projections, program_options, solve_options)
    else
        for curr_index in found_indices_right
            if V[curr_index].is_type(New)
                V[curr_index].set_type(Semicritical)
                V[curr_index].remove_type(New)
            end
        end
        for curr_index in found_indices_left
            if V[curr_index].is_type(New)
                V[curr_index].set_type(Semicritical)
                V[curr_index].remove_type(New)
            end
        end
    end

    if program_options.verbose_level() >= 0
        println("num_edges = ", num_edges_)
    end

    clear_vec_mp(particular_projection)
    clear_vec_mp(crit_downstairs)
    clear_vec_mp(mid_downstairs)

    return SUCCESSFUL
end

function MidSlice(edge_counter::Int, midpoint_witness_sets::Vector{WitnessSet}, ml_config::MultilinConfiguration, W_curve::WitnessSet, particular_projection::MPVec, mid_downstairs::MPVec, program_options::BertiniRealConfig, solve_options::SolverConfiguration)
    println("MidSlice")

    num_midpoints = mid_downstairs.size

    for ii in 1:num_midpoints
        neg_mp(particular_projection.coord[1], mid_downstairs.coord[ii])
        real_threshold(particular_projection.coord[1], solve_options.T.real_threshold)

        if program_options.verbose_level() >= 2
            println("solving midpoints upstairs $ii, projection value ", mid_downstairs.coord[ii])
        end

        solve_options.backup_tracker_config("getting_midpoints_$(ii)")

        fillme = SolverOutput()
        multilin_solver_master_entry_point(W_curve, fillme, [particular_projection], ml_config, solve_options)
        fillme.get_noninfinite_w_mult_full(midpoint_witness_sets[ii])

        num_total_midslice_points = midpoint_witness_sets[ii].num_points()

        if program_options.verbose_level() >= 4
            midpoint_witness_sets[ii].print_to_screen()
            println("midpoint_downstairs $ii had $(midpoint_witness_sets[ii].num_points()) real and complex points total")
        end

        midpoint_witness_sets[ii].sort_for_unique(program_options.same_point_tol())
        num_unique_midslice_points = midpoint_witness_sets[ii].num_points()

        if num_total_midslice_points - num_unique_midslice_points > 0
            println(color_text("there were non-unique midslice points in interval $ii.", :red))
            println("trying to recover the failure by tightening tracking tolerances...")

            solve_options.T.endgameNumber = 2
            solve_options.T.basicNewtonTol *= 1e-2
            solve_options.T.endgameNewtonTol *= 1e-2
            println("new temporary tracktolBEFOREeg: ", solve_options.T.basicNewtonTol, " tracktolDURINGeg: ", solve_options.T.endgameNewtonTol)

            fillme2 = SolverOutput()
            multilin_solver_master_entry_point(W_curve, fillme2, [particular_projection], ml_config, solve_options)

            midpoint_witness_sets[ii].clear()

            fillme2.get_noninfinite_w_mult_full(midpoint_witness_sets[ii])
        end

        midpoint_witness_sets[ii].sort_for_unique(program_options.same_point_tol())
        num_unique_midslice_points = midpoint_witness_sets[ii].num_points()

        if num_total_midslice_points - num_unique_midslice_points > 0
            println(color_text("there were non-unique midslice points in interval $ii. your decomposition is possibly incorrect about the missed points, if the path crossings obscured real points", :red))
        end

        if program_options.verbose_level() >= 4
            midpoint_witness_sets[ii].print_to_screen()
            println("midpoint_downstairs $ii had $(midpoint_witness_sets[ii].num_points()) real and complex points total")
        end

        midpoint_witness_sets[ii].sort_for_real(solve_options.T.real_threshold)
        num_real_midslice_points = midpoint_witness_sets[ii].num_points()

        if program_options.verbose_level() >= 3
            if num_total_midslice_points - num_real_midslice_points > 0
                midpoint_witness_sets[ii].print_to_screen()
                println("midpoint_downstairs $ii had $(midpoint_witness_sets[ii].num_points()) real points total")
            else
                println("all midpoints real")
            end
        end

        if have_sphere()
            midpoint_witness_sets[ii].sort_for_inside_sphere(sphere_radius(), sphere_center())
        end

        num_real_interior_midslice_points = midpoint_witness_sets[ii].num_points()

        if program_options.verbose_level() >= 3
            if num_real_midslice_points - num_real_interior_midslice_points > 0
                midpoint_witness_sets[ii].print_to_screen()
                println("midpoint_downstairs $ii had $(midpoint_witness_sets[ii].num_points()) real points inside sphere of interest")
            else
                println("all real midpoints are inside sphere")
            end
        end

        edge_counter += midpoint_witness_sets[ii].num_points()

        solve_options.restore_tracker_config("getting_midpoints_$(ii)")
    end
end

function ConnectTheDots(found_indices_crit::Vector{Set{Int}}, found_indices_mid::Vector{Set{Int}}, found_indices_left::Set{Int}, found_indices_right::Set{Int}, crit_point_counter::Dict{Int, Int}, V::VertexSet, crit_downstairs::MPVec, mid_downstairs::MPVec, particular_projection::MPVec, midpoint_witness_sets::Vector{WitnessSet}, ml_config::MultilinConfiguration, program_options::BertiniRealConfig, solve_options::SolverConfiguration)
    num_midpoints = length(midpoint_witness_sets)

    Wleft = WitnessSet()
    Wright = WitnessSet()
    cycle_nums_left = Vector{Int}()
    cycle_nums_right = Vector{Int}()

    left_proj_val = init_mp()
    right_proj_val = init_mp()

    temp_vertex = Vertex()

    edge_occurence_tracker_left = Dict{Int, Vector{Int}}()
    edge_occurence_tracker_right = Dict{Int, Vector{Int}}()

    resize!(found_indices_mid, num_midpoints)
    resize!(found_indices_crit, num_midpoints + 1)

    for ii in 1:num_midpoints
        if (program_options.verbose_level() == 0 && solve_options.path_number_modulus != 0 && ii % solve_options.path_number_modulus == 0) || (program_options.verbose_level() >= 1)
            println(color_text("connecting midpoint downstairs, $ii of $num_midpoints", :brown))
        end

        cycle_nums_left = Int[]
        cycle_nums_right = Int[]

        solve_options.backup_tracker_config("midpoint_connect")

        try_again = true
        iterations = 0
        maxits = 2

        while try_again && iterations < maxits
            iterations += 1
            try_again = false

            if program_options.verbose_level() >= 2
                print_comp_matlab(crit_downstairs.coord[ii], "left_proj_val ")
                print_comp_matlab(crit_downstairs.coord[ii + 1], "right_proj_val ")
            end

            fillme0 = SolverOutput()
            neg_mp(particular_projection.coord[1], crit_downstairs.coord[ii])
            if program_options.Realify
                midpoint_witness_sets[ii].Realify(solve_options.T.real_threshold)
            end
            multilin_solver_master_entry_point(midpoint_witness_sets[ii], fillme0, [particular_projection], ml_config, solve_options)
            fillme0.get_noninfinite_w_mult_full(Wleft)
            cycle_nums_left = fillme0.get_cyclenums_noninfinite_w_mult()

            fillme0.reset()
            fillme0 = SolverOutput()
            neg_mp(particular_projection.coord[1], crit_downstairs.coord[ii + 1])
            multilin_solver_master_entry_point(midpoint_witness_sets[ii], fillme0, [particular_projection], ml_config, solve_options)
            fillme0.get_noninfinite_w_mult_full(Wright)
            cycle_nums_right = fillme0.get_cyclenums_noninfinite_w_mult()

            Wright_real = WitnessSet(Wright)
            Wleft_real = WitnessSet(Wleft)
            Wright_real.sort_for_real(solve_options.T.real_threshold)
            Wleft_real.sort_for_real(solve_options.T.real_threshold)

            if Wleft_real.num_points() != midpoint_witness_sets[ii].num_points()
                println(color_text("had a critical failure", :red))
                println("moving left was deficient ", midpoint_witness_sets[ii].num_points() - Wleft_real.num_points(), " points")
                try_again = true
            end

            if Wright_real.num_points() != midpoint_witness_sets[ii].num_points()
                println(color_text("had a critical failure", :red))
                println("moving right was deficient ", midpoint_witness_sets[ii].num_points() - Wright_real.num_points(), " points")
                try_again = true
            end

            if !try_again
                if iterations > 1
                    println(color_text("resolution successful", :green))
                end
                break
            elseif iterations < maxits
                Wleft.reset()
                Wright.reset()
                cycle_nums_left = Int[]
                cycle_nums_right = Int[]
                println("trying to recover the failure by tightening tolerances...")
                solve_options.T.endgameNumber = 2
                solve_options.T.basicNewtonTol *= 1e-2
                solve_options.T.endgameNewtonTol *= 1e-2
                println("tracktolBEFOREeg: ", solve_options.T.basicNewtonTol, " tracktolDURINGeg: ", solve_options.T.endgameNewtonTol)
                continue
            else
                Wleft.reset_points()
                Wright.reset_points()
                cycle_nums_left = Int[]
                cycle_nums_right = Int[]

                W_single = WitnessSet(midpoint_witness_sets[ii])
                W_single_sharpened = WitnessSet()

                W_midpoint_replacement = WitnessSet(midpoint_witness_sets[ii])
                W_midpoint_replacement.reset_points()

                for kk in 1:midpoint_witness_sets[ii].num_points()
                    W_single.reset_points()
                    W_single_sharpened.reset_points()
                    W_single_right = WitnessSet()
                    W_single_left = WitnessSet()

                    W_single.add_point(midpoint_witness_sets[ii].point(kk))

                    prev_sharpen_digits = solve_options.T.sharpenDigits
                    solve_options.T.sharpenDigits = min(4 * solve_options.T.sharpenDigits, 300)
                    neg_mp(particular_projection.coord[1], mid_downstairs.coord[ii])

                    fillme1 = SolverOutput()
                    multilin_solver_master_entry_point(W_single, fillme1, [particular_projection], ml_config, solve_options)
                    fillme1.get_noninfinite_w_mult_full(W_single_sharpened)
                    fillme1.reset()

                    if W_single_sharpened.num_points() == 0
                        println("sharpening failed, which sucks because the sharpened point was theoretically generic with respect to the system currently being used")
                    end

                    solve_options.T.sharpenDigits = prev_sharpen_digits

                    W_single_left.reset_points()
                    num_its = 0
                    while num_its < 2 && W_single_left.num_points() == 0
                        W_single_left.reset_points()
                        num_its += 1
                        println("$num_its th attempt, going left, midpoint $ii")

                        if num_its > 0
                            solve_options.T.maxNewtonIts = 2
                        end

                        fillme2 = SolverOutput()
                        multilin_solver_master_entry_point(W_single_sharpened, fillme2, [particular_projection], ml_config, solve_options)
                        fillme2.get_noninfinite_w_mult_full(W_single_left)
                        c1 = fillme2.get_cyclenums_noninfinite_w_mult()
                        fillme2.reset()

                        if length(c1) == 0
                            println(color_text("tracking left yielded a non-real point", :red))
                        end
                    end

                    W_single_right.reset_points()
                    num_its = 0
                    while num_its < 2 && W_single_right.num_points() == 0
                        W_single_right.reset_points()
                        num_its += 1
                        println("$num_its th attempt, going right, midpoint $ii")

                        if num_its > 0
                            solve_options.T.maxNewtonIts = 2
                        end

                        fillme2 = SolverOutput()
                        multilin_solver_master_entry_point(W_single_sharpened, fillme2, [particular_projection], ml_config, solve_options)
                        fillme2.get_noninfinite_w_mult_full(W_single_right)
                        c2 = fillme2.get_cyclenums_noninfinite_w_mult()
                        fillme2.reset()

                        if length(c2) == 0
                            println(color_text("tracking right yielded a non-real point", :red))
                        end
                    end

                    if length(c2) == 1 && length(c1) == 1
                        W_midpoint_replacement.add_point(midpoint_witness_sets[ii].point(kk))
                        Wleft.add_point(W_single_left.point(1))
                        Wright.add_point(W_single_right.point(1))
                        push!(cycle_nums_left, c1[end])
                        push!(cycle_nums_right, c2[end])
                    else
                        temp_vertex.set_point(midpoint_witness_sets[ii].point(kk))
                        temp_vertex.set_type(Problematic)
                        index_in_vertices_with_add(V, temp_vertex)
                    end
                end

                midpoint_witness_sets[ii].reset_points()
                midpoint_witness_sets[ii].copy_points(W_midpoint_replacement)
                break
            end
        end

        solve_options.restore_tracker_config("midpoint_connect")

        for kk in 1:midpoint_witness_sets[ii].num_points()
            temp_edge = Edge()
            temp_vertex.set_point(midpoint_witness_sets[ii].point(kk))
            temp_vertex.set_type(Midpoint)

            temp_edge.midpt(index_in_vertices_with_add(V, temp_vertex))

            temp_vertex.set_point(Wleft.point(kk))
            temp_vertex.set_type(New)
            temp_edge.left(index_in_vertices_with_add(V, temp_vertex))

            temp_vertex.set_point(Wright.point(kk))
            temp_vertex.set_type(New)
            temp_edge.right(index_in_vertices_with_add(V, temp_vertex))

            found_indices_left.add(temp_edge.left())
            found_indices_right.add(temp_edge.right())
            found_indices_crit[ii].add(temp_edge.left())
            found_indices_crit[ii + 1].add(temp_edge.right())
            found_indices_mid[ii].add(temp_edge.midpt())

            md = EdgeMetaData(cycle_nums_left[kk], cycle_nums_right[kk])
            edge_num = AddEdge(temp_edge, md)
            push!(get(edge_occurence_tracker_left, temp_edge.left(), []), edge_num)
            push!(get(edge_occurence_tracker_right, temp_edge.right(), []), edge_num)

            if !haskey(crit_point_counter, temp_edge.left())
                crit_point_counter[temp_edge.left()] = 1
            else
                crit_point_counter[temp_edge.left()] += 1
            end

            if !haskey(crit_point_counter, temp_edge.right())
                crit_point_counter[temp_edge.right()] = 1
            else
                crit_point_counter[temp_edge.right()] += 1
            end

            if program_options.verbose_level() >= 2
                println("done connecting upstairs midpoint $kk (downstairs midpoint $ii)")
                if program_options.verbose_level() >= 1
                    println("constructed edge: ", temp_edge, "\n")
                end
            end
        end

        Wleft.reset()
        Wright.reset()
    end

    clear_mp(left_proj_val)
    clear_mp(right_proj_val)
end

function GetMergeCandidates(V::VertexSet)
    println("curve::GetMergeCandidates")

    default_found_edges = [-1]

    # Looking for edges with the type New, by looking at the left endpoint
    for tentative_right_edge in 1:num_edges_
        if V[edges[tentative_right_edge].left()].type == New && V[edges[tentative_right_edge].right()].type != New
            # Found a starting point for the merges
            if edges[tentative_right_edge].is_degenerate()
                println(color_text("found a degenerate edge. the comment says this should never happen.", :red))
                continue
            end

            tentative_edge_list = [tentative_right_edge]

            while true
                tentative_left_edge = nondegenerate_edge_w_right(edges[tentative_edge_list[end]].left())

                if tentative_left_edge < 0
                    println(color_text("found that edge ", :red), tentative_edge_list[end], color_text(" with points ", :red), edges[tentative_edge_list[end]], color_text(" has NEW leftpoint, but exists edge w point ", :red), edges[tentative_edge_list[end]].left(), color_text(" as right point.", :red))
                    V.GetVertex(edges[tentative_edge_list[end]].left()).print()
                    break
                end

                tentative_edge_list = push!(tentative_edge_list, tentative_left_edge)

                if V[edges[tentative_left_edge].left()].type != New
                    break
                end
            end

            if length(tentative_edge_list) > 1
                return tentative_edge_list
            else
                continue
            end
        end
    end

    return default_found_edges
end

function Merge(W_midpt::WitnessSet, V::VertexSet, projections::Vector{vec_mp}, program_options::BertiniRealConfig, solve_options::SolverConfiguration)
    println("curve::merge")

    particular_projection = init_vec_mp(0)
    vec_cp_mp(particular_projection, projections[1])

    half = init_mp2(1024)
    temp = init_mp2(1024)
    temp2 = init_mp2(1024)
    mpf_set_str(half.r, "0.5", 10)
    mpf_set_str(half.i, "0.0", 10)

    ml_config = MultilinConfiguration(solve_options)

    edges_to_merge = GetMergeCandidates(V)

    while edges_to_merge[end] != -1
        rightmost_edge = edges_to_merge[1]
        leftmost_edge = edges_to_merge[end]
        moving_edge = edges_to_merge[div(length(edges_to_merge), 2) + 1]

        if solve_options.verbose_level >= 1
            println(color_text("merging edges: ", :cyan), reverse(edges_to_merge))
        end

        if edges_to_merge[end] < 0
            println("error: attempting to merge an edge with negative index!")
            break
        end

        W_temp = WitnessSet()

        projection_value_homogeneous_input!(particular_projection.coord[1], V[edges[moving_edge].midpt()].point, projections[1])
        neg_mp(particular_projection.coord[1], particular_projection.coord[1])

        reset_linears!(W_midpt)
        add_linear!(W_midpt, particular_projection)

        reset_points!(W_midpt)
        add_point!(W_midpt, V[edges[moving_edge].midpt()].point)

        if program_options.Realify
            Realify!(W_midpt, solve_options.T.real_threshold)
        end

        projection_value_homogeneous_input!(temp, V[edges[leftmost_edge].left()].point, projections[1])
        projection_value_homogeneous_input!(temp2, V[edges[rightmost_edge].right()].point, projections[1])

        add_mp!(new_proj_val, temp, temp2)
        mul_mp!(new_proj_val, new_proj_val, half)
        neg_mp(particular_projection.coord[1], new_proj_val)

        ml_config.set_randomizer(randomizer())
        fillme = SolverOutput()
        multilin_solver_master_entry_point(W_midpt, fillme, particular_projection, ml_config, solve_options)

        get_noninfinite_w_mult_full!(W_temp, fillme)

        if W_temp.num_points == 0
            println(color_text("merging multilin solver returned NO POINTS!!!  unable to continue merging, sorry.", :red))
            break
        end

        temp_vertex = Vertex()
        set_point!(temp_vertex, W_temp.point(1))
        set_type!(temp_vertex, Midpoint)

        temp_edge = Edge()
        left(temp_edge, edges[leftmost_edge].left())
        midpt(temp_edge, index_in_vertices_with_add(V, temp_vertex))
        right(temp_edge, edges[rightmost_edge].right())

        for zz in 1:length(edges_to_merge)
            merge_me_away = edges_to_merge[zz]

            for vec_iter in edges[merge_me_away].removed
                add_removed_point!(temp_edge, vec_iter)
            end

            if zz == 1
                add_removed_point!(temp_edge, edges[merge_me_away].left())
                add_removed_point!(temp_edge, edges[merge_me_away].midpt())
                V[edges[merge_me_away].midpt()].set_removed(true)
                V[edges[merge_me_away].left()].set_removed(true)
            elseif zz == length(edges_to_merge)
                add_removed_point!(temp_edge, edges[merge_me_away].midpt())
                V[edges[merge_me_away].midpt()].set_removed(true)
            else
                add_removed_point!(temp_edge, edges[merge_me_away].left())
                add_removed_point!(temp_edge, edges[merge_me_away].midpt())
                V[edges[merge_me_away].midpt()].set_removed(true)
                V[edges[merge_me_away].left()].set_removed(true)
            end
        end

        md = EdgeMetaData(edge_metadata[leftmost_edge].CycleNumLeft, edge_metadata[rightmost_edge].CycleNumRight)
        AddEdge(temp_edge, md)

        post_merge_edges = Edge[]
        post_merge_metadata = EdgeMetaData[]

        num_removed_edges = 0

        for ii in 1:num_edges_
            remove_flag = false

            for zz in 1:length(edges_to_merge)
                if edges_to_merge[zz] == ii
                    remove_flag = true
                    break
                end
            end

            if !remove_flag
                push!(post_merge_edges, edges[ii])
                push!(post_merge_metadata, edge_metadata[ii])
            else
                num_removed_edges += 1
            end
        end

        if num_removed_edges != length(edges_to_merge)
            error("claiming to have merged away $num_removed_edges edges, but had $(length(edges_to_merge)) in the list to merge.")
        end

        edges_ = post_merge_edges
        edge_metadata_ = post_merge_metadata
        num_edges_ = length(edges_)

        edges_to_merge = GetMergeCandidates(V)
    end

    clear_mp(half)
    clear_mp(temp)
    clear_mp(temp2)
    clear_mp(new_proj_val)
    clear_vec_mp(particular_projection)
end

function verify_projection_ok(W::WitnessSet, projection::Vector{vec_mp}, solve_options::SolverConfiguration)
    randomizer = SystemRandomizer()
    setup(randomizer, W.num_variables() - W.num_patches() - W.dimension(), solve_options.PPD.num_funcs)

    invalid_flag = verify_projection_ok(W, randomizer, projection, solve_options)

    return invalid_flag
end

function verify_projection_ok(W::WitnessSet, randomizer::SystemRandomizer, projection::Vector{vec_mp}, solve_options::SolverConfiguration)
    invalid_flag = 0

    parse_input_file(W.input_filename())

    temp_rand_point = init_vec_mp(W.num_variables(), 0)
    set_one_mp(temp_rand_point.coord[1])

    for ii in 2:W.num_variables()
        get_comp_rand_mp(temp_rand_point.coord[ii])
    end

    SLP = prog_t()
    setupProg(SLP, solve_options.T.Precision, 2)

    zerotime = init_mp()
    set_zero_mp(zerotime)

    ED = init_eval_struct_mp(0, 0, 0)
    evalProg_mp(ED.funcVals, ED.parVals, ED.parDer, ED.Jv, ED.Jp, temp_rand_point, zerotime, SLP)

    AtimesJ = init_mat_mp(1, 1)
    randomizer.randomize(temp_rand_point, AtimesJ, ED.funcVals, ED.Jv, temp_rand_point.coord[1])

    detme = init_mat_mp(W.num_variables() - 1, W.num_variables() - 1)

    for ii in 1:AtimesJ.rows
        for jj in 1:AtimesJ.cols - 1
            set_mp(detme.entry[ii, jj], AtimesJ.entry[ii, jj + 1])
        end
    end

    offset = W.num_variables() - 1 - W.dimension()

    for jj in 1:W.dimension()
        for ii in 1:W.num_variables() - 1
            set_mp(detme.entry[offset + jj, ii], projection[jj].coord[ii + 1])
        end
    end

    determinant = init_mp()
    take_determinant_mp(determinant, detme)

    if d_abs_mp(determinant) < 1e-2
        invalid_flag = 0
        println(color_text("determinant test revealed that your projection is probably inappropriate, and you should choose another one.", :red))
        println(d_abs_mp(determinant))
        print_matrix_to_screen_matlab(ED.Jv, "Jv")
        print_matrix_to_screen_matlab(detme, "detme")
    else
        invalid_flag = 1
    end

    clear_mat_mp(detme)
    clear_mat_mp(AtimesJ)
    clear_vec_mp(temp_rand_point)
    clear_mp(determinant)
    clear_mp(zerotime)
    clear_eval_struct_mp(ED)
    clearProg(SLP, solve_options.T.MPType, 1)

    return invalid_flag
end

function setup(containing_folder::String)
    setup_decomposition(containing_folder * "/decomp")
    setup_edges(containing_folder * "/E.edge")
    setup_cycle_numbers(containing_folder * "/curve.cnums")
    return 1
end

function setup_edges(INfile::String)
    println("curve::setup_edges")

    IN = safe_fopen_read(INfile)

    temp_num_edges = parse(Int, readline(IN))
    for ii in 1:temp_num_edges
        left, midpt, right = parse.(Int, split(readline(IN)))
        AddEdge(Edge(left, midpt, right), EdgeMetaData())
    end

    close(IN)
    return num_edges_()
end

function setup_cycle_numbers(INfile::String)
    println("curve::setup_cycle_numbers")

    if !isfile(INfile)
        return 0
    end

    IN = safe_fopen_read(INfile)

    temp_num_edges = parse(Int, readline(IN))

    if temp_num_edges != num_edges()
        throw(RuntimeError("mismatch in number of cycle number data, and number of edges in curve in file $INfile"))
    end

    for ii in 1:temp_num_edges
        left, right = parse.(Int, split(readline(IN)))
        edge_metadata_[ii] = EdgeMetaData(left, right)
    end

    close(IN)
    return num_edges_()
end

function print(base::String)
    println("curve::print")
    print_decomposition(base)
    edgefile = base * "/E.edge"
    print_edges(edgefile)
    cycle_num_file = base * "/curve.cnums"
    print_cycle_numbers(cycle_num_file)
end

function print_edges(outputfile::String)
    println("curve::print_edges")
    OUT = safe_fopen_write(outputfile)

    # output the number of vertices
    println(OUT, "$(num_edges_())\n")

    for ii in 1:num_edges_()
        println(OUT, "$(edges_[ii].left()) $(edges_[ii].midpt()) $(edges_[ii].right())")
    end

    close(OUT)
end

function print_cycle_numbers(outputfile::String)
    println("curve::print_cycle_numbers")
    OUT = safe_fopen_write(outputfile)

    # output the number of vertices
    println(OUT, "$(num_edges_())\n")

    for ii in 1:num_edges_()
        println(OUT, "$(edge_metadata_[ii].CycleNumLeft()) $(edge_metadata_[ii].CycleNumRight())")
    end

    close(OUT)
end


function computeCurveNotSelfConj(W_in::WitnessSet, V::VertexSet, num_vars::Int, program_options::BertiniRealConfig, solve_options::SolverConfiguration)
    println("curve::computeCurveNotSelfConj")

    IN = nothing

    declarations = nothing
    partition_parse(&declarations, W_in.input_filename(), "func_input_nsc", "config_nsc", 1)
    free(declarations)

    diag_homotopy_input_file("input_NSC", "func_input_nsc", "func_inputbar", "config_nsc", W_in.linear(0), num_vars - 1)
    diag_homotopy_start_file("start_NSC", W_in)

    copyfile("witness_data", "witness_data_0")

    command_line_options = ["input_NSC", "start_NSC"]

    program_options.call_for_help(BERTINI_MAIN)
    bertini_main_wrapper(command_line_options, program_options.num_procs(), 0, 0)

    rename("witness_data_0", "witness_data")

    IN = safe_fopen_read("real_solutions")

    num_sols = parse(Int, readline(IN))

    temp_vertex = Vertex()
    change_size_vec_mp(temp_vertex.point(), num_vars)
    temp_vertex.point().size = num_vars
    temp_vertex.set_type(Isolated)

    cur_sol = Array{mp}(undef, num_vars)
    init_vec_mp(cur_sol, num_vars)
    cur_sol[1] = parse(mp, "1")

    cur_sol_bar = Array{mp}(undef, num_vars)
    init_vec_mp(cur_sol_bar, num_vars)
    cur_sol_bar[1] = parse(mp, "1")

    for ii in 1:num_sols
        for jj in 1:num_vars - 1
            cur_sol[jj + 1] = parse(mp, readline(IN))
            cur_sol_bar[jj + 1] = parse(mp, readline(IN))
        end

        if isSamePoint_homogeneous_input(cur_sol, cur_sol_bar, solve_options.T.final_tol_times_mult)
            temp_vertex.set_point(cur_sol)
            index_in_vertices_with_add(V, temp_vertex)
        end
    end

    fclose(IN)

    clear_vec_mp(cur_sol)
    clear_vec_mp(cur_sol_bar)

    remove("func_input_nsc")
    remove("config_nsc")
    remove("func_inputbar")
    remove("var_names")
end

using Printf

function diag_homotopy_input_file(outputFile::String, funcInputx::String, funcInputy::String, configInput::String, L::Array{mp}, num_vars::Int)
    println("diag_homotopy_input_file")

    str = fill("", num_vars)
    for ii in 1:num_vars
        str[ii] = ""
    end

    OUT = safe_fopen_write(outputFile)

    # setup configurations in OUT
    println(OUT, "CONFIG")
    IN = safe_fopen_read(configInput)
    close(IN)
    println(OUT, "USERHOMOTOPY: 1;")
    println(OUT, "DeleteTempFiles: 0;")
    println(OUT, "END;")
    println(OUT, "INPUT")

    # setup variables in OUT
    IN = safe_fopen_read(funcInputx)
    write(OUT, read(IN, String))
    close(IN)

    # setup the function name in OUT
    IN = safe_fopen_read(funcInputy)
    write(OUT, read(IN, String))
    close(IN)

    IN = safe_fopen_read("var_names")
    ii = 1
    jj = 1
    while !eof(IN)
        ch = Char(read(IN, UInt8))
        if ch != '\n'
            str[ii] *= ch
        else
            ii += 1
            jj = 1
        end
    end
    close(IN)

    # setup the linear equations
    temp = Array{mp}(undef, L.size)
    for ii in 1:L.size
        println(OUT, "bertini_real_L$ii = ", mpf_out_str(L.coord[ii].r), "+I*", mpf_out_str(L.coord[ii].i), ";")
        println(OUT, "bertini_real_Lbar$ii = ", mpf_out_str(conjugate_mp(temp[ii], L.coord[ii]).r), "+I*", mpf_out_str(temp[ii].i), ";")
    end
    println(OUT)

    # Generate a random matrix A and output to input file.
    A = Array{mp}(undef, 2, num_vars)
    make_matrix_random_d(A, 2, num_vars)
    for ii in 1:2
        for jj in 1:num_vars
            println(OUT, "A$ii$jj = ", mpf_out_str(A[ii, jj].r), "+I*", mpf_out_str(A[ii, jj].i), ";")
        end
    end

    # setup the diagonal homotopy functions
    fmt = @sprintf("%%.%dlf+%%.%dlf*I", 15, 15)
    println(OUT, "bertini_real_L=t*(")
    for ii in 1:num_vars
        print(OUT, "bertini_real_L$ii*", str[ii], "+")
    end
    println(OUT, "-1)+(1-t)*(")
    for ii in 1:num_vars
        print(OUT, "A01*(", str[ii], "-$str[ii]bar)+")
    end
    println(OUT, "0);")
    println(OUT, "bertini_real_Lbar=t*(")
    for ii in 1:num_vars
        print(OUT, "bertini_real_Lbar$ii*", str[ii], "bar+")
    end
    println(OUT, "-1)+(1-t)*(")
    for ii in 1:num_vars
        print(OUT, "A1$ii*(", str[ii], "-$str[ii]bar)+")
    end
    println(OUT, "0);")
    println(OUT, "END;")

    close(OUT)

    return
end

using Printf

function diag_homotopy_start_file(startFile::String, W::WitnessSet)
    OUT = safe_fopen_write(startFile)

    # output the number of start points
    println(OUT, "${W.num_points()*W.num_points()}\n")

    temp = Array{mp}(undef)
    init_mp(temp)

    result = Vector{mp}(undef, W.num_variables() - 1)
    init_vec_mp(result, 0)

    result2 = Vector{mp}(undef, W.num_variables() - 1)
    init_vec_mp(result2, 0)

    for ii in 1:W.num_points()
        outer_point = W.point(ii)

        change_prec_vec_mp(result, outer_point.curr_prec)
        dehomogenize(result, outer_point)

        for jj in 1:W.num_points() # output {w \bar{w}}'
            inner_point = W.point(jj)

            change_prec_vec_mp(result, inner_point.curr_prec)
            dehomogenize(result2, inner_point)

            change_prec_mp(temp, inner_point.curr_prec)

            for kk in 1:W.num_variables() - 1
                print_mp(OUT, 0, result[kk])
                println(OUT)
                conjugate_mp(temp, result2[kk])
                print_mp(OUT, 0, temp)
                println(OUT)
            end
            println(OUT)
        end
    end

    clear_vec_mp(result2)
    clear_vec_mp(result)
    clear_mp(temp)

    close(OUT)
end

function ProjectionIntervalIndex(edge_index::Int, V::VertexSet)
    if edge_index >= num_edges()
        error("edge index $edge_index exceeds number of stored edges $(num_edges())")
    end

    temp = Array{mp}(undef)
    init_mp(temp)

    v = V.GetVertex(edges_[edge_index].left())
    ps = v.projection_values()

    minval = 1e200
    loc = -1

    for ii in 1:crit_slice_values.size
        sub_mp(temp, ps.coord[1], crit_slice_values.coord[ii])
        c = d_abs_mp(temp)
        if c < minval
            minval = c
            loc = ii
        end
    end

    if minval > 1e-5
        println(stdout, "$(color.red())returned index for projection interval index is almost certainly wrong\n\n$(color.console_default())")
    end

    clear_mp(temp)
    return loc
end
