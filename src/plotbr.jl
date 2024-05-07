module plotbrModule 

export plotbr, Options, Curve
include("decomposition.jl")
using .decompositionModule
include("vertex.jl")
using .vertexModule
include("curve.jl")
using .curveModule

using Plots

@enum ColorMode begin
    BY_CELL = 1
    MONO = 2
    BY_FUNCTION = 3
end

# StyleOptions
mutable struct StyleOptions
    line_thickness::Int
    colormap::Any
    colormode::Any
    mono_color::Union{Nothing, Color}
    color_function::Union{Nothing, Function}

    function StyleOptions()
        obj = new()
        set_defaults(obj)
        return obj
    end
end

function set_defaults(options::StyleOptions)
    options.line_thickness = 2
    options.colormap = cgrad(:viridis)  # not a function
    options.colormode = 1 # BY_CELL
    options.mono_color = nothing
    options.color_function = nothing
end

function set_color_function(options::StyleOptions, fun::Function)
    options.colormode = 3 # BY_FUNCTION
    options.color_function = fun
end

# VisibilityOptions
mutable struct VisibilityOptions
    vertices::Bool
    vertices_by_type::Dict{String, Bool}
    surface_samples::Bool
    surface_raw::Bool
    curve_samples::Bool
    curve_raw::Bool
    labels::Bool

    function VisibilityOptions()
        obj = new()
        set_defaults(obj)
        return obj
    end
end

function set_defaults(options::VisibilityOptions)
    options.vertices = false
    names = [string(split(string(t), '.')[1]) for t in keys(VertexType)]
    options.vertices_by_type = Dict(n => true for n in names)
    options.surface_samples = false
    options.surface_raw = false
    options.curve_samples = false
    options.curve_raw = false
    options.labels = false
end

function auto_adjust(options::VisibilityOptions, decomposition)
    if isa(decomposition, curveModule.Curve)
        adjust_for_curve(options, decomposition)    
    #elseif isa(decomposition, SurfacePiece)
    #    adjust_for_piece(options, decomposition)
    #elseif isa(decomposition, Surface)
    #    adjust_for_surface(options, decomposition)
    else
        error("NotImplementedError: cannot auto_adjust VisibilityOptions for dimension $(decomposition.dimension) components")
    end
end


function adjust_for_curve(options::VisibilityOptions, curve)
    if length(curve.decomposition.vertices) > 10000
        options.vertices = false
    end

    if curve.sampler_data === nothing
        options.curve_raw = true
        options.curve_samples = false
    else
        options.curve_raw = false
        options.curve_samples = true
    end
end

# RenderOptions
mutable struct RenderOptions
    vertices::Bool
    surface_samples::Bool
    surface_raw::Bool
    curve_samples::Bool
    curve_raw::Bool
    labels::Bool
    which_faces::Vector{Int}
    which_edges::Vector{Int}
    defer_show::Bool

    function RenderOptions()
        obj = new()
        set_defaults(obj)
        return obj
    end
end

function set_defaults(ro::RenderOptions)
    ro.vertices = true
    ro.surface_samples = true
    ro.surface_raw = true
    ro.curve_samples = true
    ro.curve_raw = true
    ro.labels = true
    ro.which_faces = Int[]
    ro.which_edges = Int[]
    ro.defer_show = false
end

function auto_adjust(ro::RenderOptions, decomposition)
    if isa(decomposition, Curve)
        adjust_for_curve(options, decomposition)    
    #elseif isa(decomposition, SurfacePiece)
    #    adjust_for_piece(options, decomposition)
    #elseif isa(decomposition, Surface)
    #    adjust_for_surface(options, decomposition)
    else
        error("NotImplementedError: cannot auto_adjust VisibilityOptions for dimension $(decomposition.dimension) components")
    end
end

function adjust_for_curve(ro::RenderOptions, curve)
    if isempty(ro.which_edges)
        ro.which_edges = collect(1:curve.num_edges)
    end

    if curve.sampler_data === nothing
        ro.curve_raw = true
        ro.curve_samples = false
    else
        ro.curve_raw = false
        ro.curve_samples = true
    end
end

# Options
mutable struct Options
    style::StyleOptions 
    visibility::VisibilityOptions  
    render::RenderOptions

    function Options()
        new(StyleOptions(), VisibilityOptions(), RenderOptions())
    end
end

# Plotter
mutable struct Plotter
    options::Any
    fig::Union{Nothing, Any}
    ax::Union{Nothing, Any}
    nondegen ::Union{Nothing, Any}
    plotted_decompositions::Vector
    widgets::Dict
    plot_results::Dict

    Plotter(options=Options()) = new(options, nothing, nothing, nothing, [], Dict("buttons" => Dict()), Dict("vertices" => Dict()))
end

function show(plotter::Plotter)
    display(plotter.fig)
      # Adjusted for Julia's display system
end


function plott(plotter::Plotter, decomposition)
    if plotter.fig === nothing
        make_new_figure(plotter)
        make_new_axes(plotter, decomposition)
        label_axes(plotter, decomposition)
        apply_title(plotter)
    end
    #TODO
    #auto_adjust_visibility(plotter.options, decomposition)
    #auto_adjust_render(plotter.options, decomposition)


    main(plotter, decomposition)

    if !plotter.options.render.defer_show
        #adjust_all_visibility(plotter)
        show(plotter)
    end 
end

function main(plotter::Plotter, decomposition)
    if isa(decomposition, Curve)
        plot_curve(plotter, decomposition) 
    #elseif isa(decomposition, AbstractArray) && all([isa(p, SurfacePiece) for p in decomposition])
    #    plot_pieces(plotter, decomposition)
    #elseif isa(decomposition, SurfacePiece)
    #    plot_piece(plotter, decomposition)
     
    #elseif isa(decomposition, Surface)
    #    plot_surface(plotter, decomposition)  # Assume `plot_surface` is defined elsewhere
    else
        print(Curve)
        error("NotImplementedError")
    end
end

# Plot functions
function make_new_figure(plotter::Plotter)
    plotter.fig = plot(size=(1000, 800))  # Size in pixels
end

function make_new_axes(plotter::Plotter, decomposition)
    if decomposition.decomposition.num_variables == 2
        plotter.fig = plotter.fig  
    else
        plotter.fig= plot(plotter.fig, projection=:`3d`)
    end

    try
        plotter.fig= plot(plotter.fig, aspect_ratio=:equal)
    catch e
        @warn "Setting aspect ratio failed, using `auto` instead" exception=(e)
        plotter.fig= plot(plotter.fig, aspect_ratio=:auto)
    end

    adjust_axis_bounds(plotter, decomposition)
end

function adjust_axis_bounds(plotter::Plotter, decomposition)
    xlims!(plotter.fig, (decomposition.decomposition.center[1] - decomposition.decomposition.radius, decomposition.decomposition.center[1] + decomposition.decomposition.radius))
    ylims!(plotter.fig, (decomposition.decomposition.center[2] - decomposition.decomposition.radius, decomposition.decomposition.center[2] + decomposition.decomposition.radius))

    if decomposition.decomposition.num_variables == 3
        zlims!(plotter.fig, (decomposition.decomposition.center[3] - decomposition.decomposition.radius, decomposition.decomposition.center[3] + decomposition.decomposition.radius))
    end
end

function apply_title(plotter::Plotter)
    title!(plotter.fig, split(pwd(), Base.Filesystem.path_separator)[end])
end

function label_axes(plotter::Plotter, decomposition)
    xlabel!(plotter.fig, "x")
    ylabel!(plotter.fig, "y")
    if decomposition.decomposition.dimension == 3
        zlabel!(plotter.fig, "z")
    end
end

# Vertices
function plot_vertices(plotter::Plotter, decomposition)
    xs, ys, zs = make_xyz(decomposition)

    # Define marker styles, assuming a way to map VertexType to markers is defined
    markers = Dict(
    "unset" => :circle,
    "critical" => :star,
    "semicritical" => :diamond,
    "midpoint" => :xcross,
    "isolated" => :square,
    "new" => :triangle,
    "curve_sample_point" => :star5,
    "surface_sample_point" => :star6,
    "removed" => :hexagon,
    "problematic" => :octagon,
    "singular" => :pentagon ) # Example mapping

    plotter.plot_results["vertices"] = Dict()

    for (T, m) in markers
        plot_these = [is_of_type(v, VertexType[T]) for v in decomposition.decomposition.vertices]
        if !any(plot_these)
            continue
        end

        if decomposition.decomposition.num_variables == 2
            scatter!(plotter.fig, xs[plot_these], ys[plot_these], marker=m)
        else
            scatter3d!(plotter.fig, xs[plot_these], ys[plot_these], zs[plot_these], markeralpha=1, marker=m)
        end

        #plotter.plot_results["vertices"][h] = T
    end

    #make_widgets_vertices!(plotter, decomposition) # Assuming this function exists
end

function make_xyz(decomposition)
    xs, ys, zs = Float64[], Float64[], Float64[]

    for v in decomposition.decomposition.vertices
        push!(xs, real(v.point[1]))
        push!(ys, real(v.point[2]))
        if decomposition.decomposition.num_variables > 2
            push!(zs, real(v.point[3]))
        end
    end

    return xs, ys, zs
end

# Curve
function plot_curve(plotter::Plotter, curve)
    #push!(plotter.plotted_decompositions, curve)
    determine_nondegen_edges(plotter, curve)

    if plotter.options.render.curve_raw
        plot_raw_edges(plotter, curve)
        #adjust_visibility(plotter, 'curveRaw')
    end

    if plotter.options.render.curve_samples
        #TODO
        #plot_edge_samples(plotter, curve)
        #adjust_visibility(plotter, 'curveSamples')
    end

    if plotter.options.render.vertices && !curve.is_embedded
        plot_vertices(plotter, curve)  
        #TODO
        #adjust_visibility(plotter, 'vertices')
    end

    if !curve.is_embedded
        #TODO
        #make_widgets_curve(plotter, curve)  
    end
end

function plot_raw_edges(plotter::Plotter, curve)
    num_nondegen_edges = length(plotter.nondegen)
    colormap = plotter.options.style.colormap
    color_list = [colormap[i] for i in range(0, stop=1, length=num_nondegen_edges)]

    for (index, edge_index) in enumerate(plotter.nondegen)
        color = color_list[index]
        xs, ys, zs = Float64[], Float64[], Float64[]
        inds = curve.edges[edge_index,:]
        for i in inds
            v = curve.decomposition.vertices[i+1]
            push!(xs, real(v.point[1]))
            push!(ys, real(v.point[2]))
            if curve.decomposition.num_variables > 2
                push!(zs, real(v.point[3]))
            end
        end
        
        if curve.decomposition.num_variables == 2
            plot!(plotter.fig, xs, ys, color=color, legend=false)
        else
            plot!(plotter.fig, xs, ys, zs, color=color, legend=false)
        end
        
    end
end

function determine_nondegen_edges(plotter::Plotter, decomposition)
    plotter.nondegen = Int[]
    for i in 1:decomposition.num_edges
        e = decomposition.edges[i,:]
        if e[1] != e[2] && e[2] != e[3] && e[1] != e[3]
            push!(plotter.nondegen, i)
        end
    end
end

# Finish
function plotbr(data, options=Options())
    plotter = Plotter(options)
    plott(plotter,data)
end

end