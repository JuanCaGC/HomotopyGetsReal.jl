# src/Visuals.jl
#
# Phase 6: plotting interface for Phase 1-5's ACTUAL output types, not
# the roadmap's originally-sketched signatures (which predate those
# types) -- confirmed directly against source before writing anything
# here (see the Phase 6 investigation): `decompose_1d_curve` returns
# `(vertices::Vector{NativeVertex{T}}, edges::Vector{Edge{T}})`
# (`src/Topology.jl`); `decompose_3d_surface` returns
# `(vertices, edges, faces::Vector{Face{T}}, mesh::GeometryBasics.Mesh)`
# (`src/SurfaceDecomposition.jl`), where `mesh` is `weld_mesh`'s output:
# ALREADY globally welded and ALREADY had its triangle winding corrected
# to align with `∇f` (`weld_mesh`'s own `fixed_triangles` step). `Face`'s
# own `mesh_topology` (`src/Types.jl`) is only locally, consistently
# wound WITHIN one face -- the winding correction is a `weld_mesh`-only
# step that needs `patch`/`F`, so a plotting function built directly
# from `Vector{Face{T}}` alone cannot reproduce it. This is why
# `mesh::GeometryBasics.Mesh` (not `faces::Vector{Face{T}}`) is the
# PRIMARY input for `plot_surface_decomposition` below, with the
# `faces`-based method kept as an explicitly-flagged secondary path
# (BertiniReal's own "color by cell" convention, see below).
#
# Inspired by, but not a literal port of, BertiniReal's Python plotting
# (`codigo_cplusplus/python/bertini_real/plot/__init__.py`):
#   - BertiniReal does NOT color-code vertices by type; it distinguishes
#     `VertexType` by MARKER SHAPE (cycling matplotlib's default marker
#     list) so a checkbox panel can toggle visibility per type
#     independently, and reserves COLOR for a separate per-cell/per-
#     function purpose (`ColorMode.BY_CELL`/`BY_FUNCTION`). We use BOTH
#     a fixed color AND a fixed marker per `VertexType` below (a cheap,
#     strictly-more-informative combination of BertiniReal's marker-only
#     convention and the roadmap's color-only sketch), and generalize
#     the roadmap's "colormap based on z-coordinate" to accept an
#     arbitrary `Function(x,y,z)`, matching BertiniReal's `BY_FUNCTION`.
#   - `VertexType` here is 4-valued (`Critical, Boundary, Singular,
#     Artificial`, `src/Types.jl`), not the 3 the roadmap's sketch names
#     -- `Artificial` is a real, load-bearing vertex type (from
#     `Topology._resolve_endpoint`'s fallback and `cluster_vertices`
#     merges), not a corner case, so it gets its own color/marker below.
#   - Deliberate deviation from the roadmap's literal color list:
#     `Singular => :orange`, not `:yellow` -- plain yellow scatter
#     points have poor contrast against Makie's default white
#     background.
#   - BertiniReal's `BY_CELL` default (every edge/face gets a distinct
#     color cycled from a colormap, so topologically distinct cells are
#     visually distinguishable) is carried over for both curve edges and
#     surface faces via `Makie.cgrad(colormap, n; categorical = true)`.
#
# Explicitly out of scope (see the Phase 6 investigation for why):
# raw-vs-smooth curve toggling (decompose_1d_curve never returns the
# pre-`sample_edge` data), live/interactive re-slicing at an
# interactively-chosen z (BertiniReal doesn't do this either -- a real
# performance/design question, not a rendering one), and
# `separate_into_nonsingular_pieces`/OBJ-STL export (a new topological
# feature, not plotting).
#
# GLMakie confirmed (empirically, not assumed) to render headlessly in
# this project's dev sandbox before writing anything below -- no
# CairoMakie fallback needed, `Project.toml`'s existing `GLMakie`
# dependency (already listed, never previously `using`d) is sufficient.

const _VERTEX_COLORS = Dict(
    Critical => :red,
    Boundary => :blue,
    Singular => :orange,
    Artificial => :gray,
)

const _VERTEX_MARKERS = Dict(
    Critical => :circle,
    Boundary => :utriangle,
    Singular => :diamond,
    Artificial => :xcross,
)

"""
    _plot_vertices_by_type!(ax, vertices::Vector{NativeVertex{T}}) where {T<:AbstractFloat}
        -> Dict{VertexType,Any}

Shared helper behind both [`plot_curve_decomposition`](@ref) (2D,
`length(coordinates) == 2`) and [`plot_surface_decomposition`](@ref)'s
`mesh` method (3D, `length(coordinates) == 3`, via its optional
`vertices` overlay keyword) -- one `scatter!` call per [`VertexType`](@ref)
ACTUALLY PRESENT (skipping empty types entirely, matching BertiniReal's
own `if not np.any(plot_these): continue`), colored/marked per
`_VERTEX_COLORS`/`_VERTEX_MARKERS`, labeled for `axislegend`. Returns the
per-type scatter handles (empty `Dict` if `vertices` is empty), so
callers can decide whether a legend is worth adding.
"""
function _plot_vertices_by_type!(ax, vertices::Vector{NativeVertex{T}}) where {T<:AbstractFloat}
    handles = Dict{VertexType,Any}()
    isempty(vertices) && return handles
    is_3d = length(vertices[1].coordinates) == 3

    for vtype in (Critical, Boundary, Singular, Artificial)
        these = filter(v -> v.v_type == vtype, vertices)
        isempty(these) && continue

        xs = T[real(v.coordinates[1]) for v in these]
        ys = T[real(v.coordinates[2]) for v in these]
        color = _VERTEX_COLORS[vtype]
        marker = _VERTEX_MARKERS[vtype]
        label = string(vtype)

        h = if is_3d
            zs = T[real(v.coordinates[3]) for v in these]
            scatter!(ax, xs, ys, zs; color = color, marker = marker, markersize = 12, label = label)
        else
            scatter!(ax, xs, ys; color = color, marker = marker, markersize = 12, label = label)
        end
        handles[vtype] = h
    end
    return handles
end

"""
    plot_curve_decomposition(
        vertices::Vector{NativeVertex{T}},
        edges::Vector{Edge{T}};
        cfg::Union{Nothing,HomotopyConfig} = nothing,
        show_vertices::Bool = true,
        show_labels::Bool = false,
        edge_color_by::Symbol = :cell,
        edge_color = :steelblue,
        colormap = :viridis,
    ) where {T<:AbstractFloat}
        -> Makie.Figure

Renders a 2D curve decomposition -- the direct output of
[`Topology.decompose_1d_curve`](@ref) or [`SurfaceDecomposition.slice_at_z`](@ref)
(both already return exactly `(vertices, edges)` in this shape; no
adapter needed).

- `cfg`, if given, sets axis limits from `cfg.bbox_x`/`cfg.bbox_y` rather
  than Makie's auto-fit-to-data, mirroring BertiniReal's own
  `_adjust_axis_bounds` (itself derived from the same kind of
  known-bounding-region config the decomposition was computed against).
- `edge_color_by`: `:cell` (default, BertiniReal's `BY_CELL`) gives every
  edge its own color cycled from `colormap`; `:mono` uses a flat
  `edge_color` for every edge.
- `show_vertices`: one legended `scatter!` per [`VertexType`](@ref)
  present, via [`_plot_vertices_by_type!`](@ref).
- `show_labels`: annotates each vertex with its `id` via `text!`.

Returns the `Figure` WITHOUT calling `display` on it -- composable with
a caller's own layout, or `save(path, fig)`, without forcing a window
open. Contrast [`interactive_3d_viewer`](@ref), whose entire contract IS
opening a window.
"""
function plot_curve_decomposition(
    vertices::Vector{NativeVertex{T}},
    edges::Vector{Edge{T}};
    cfg::Union{Nothing,HomotopyConfig} = nothing,
    show_vertices::Bool = true,
    show_labels::Bool = false,
    edge_color_by::Symbol = :cell,
    edge_color = :steelblue,
    colormap = :viridis,
) where {T<:AbstractFloat}
    edge_color_by in (:cell, :mono) || throw(ArgumentError(
        "plot_curve_decomposition: edge_color_by must be :cell or :mono, got $(edge_color_by)",
    ))

    fig = Figure()
    ax = Axis(fig[1, 1]; aspect = DataAspect(), xlabel = "x", ylabel = "y")

    if cfg !== nothing
        xlims!(ax, cfg.bbox_x...)
        ylims!(ax, cfg.bbox_y...)
    end

    n_edges = length(edges)
    edge_colors = edge_color_by == :cell ? Makie.cgrad(colormap, max(n_edges, 2); categorical = true) : nothing
    for (i, e) in enumerate(edges)
        isempty(e.sampled_points) && continue
        xs = T[p[1] for p in e.sampled_points]
        ys = T[p[2] for p in e.sampled_points]
        c = edge_color_by == :cell ? edge_colors[i] : edge_color
        lines!(ax, xs, ys; color = c)
    end

    if show_vertices
        handles = _plot_vertices_by_type!(ax, vertices)
        isempty(handles) || axislegend(ax)
    end

    if show_labels
        for v in vertices
            x, y = T(real(v.coordinates[1])), T(real(v.coordinates[2]))
            text!(ax, [x], [y]; text = [string(v.id)], fontsize = 10, offset = (4, 4))
        end
    end

    return fig
end

const _NEAR_CONSTANT_COLOR_WARNED = Ref(false)

# Below this relative spread, `color_by`'s range is treated as numerical
# noise rather than a meaningful signal -- see the `_near_constant_colorrange`
# docstring for the empirical motivation (found via the 02b investigation:
# a mathematically-exact-1.0 quantity's Float32 round-off, ~6e-8 relative,
# auto-scaled a colorbar into full-spectrum speckle).
const _COLOR_BY_MIN_REL_RANGE = 1e-4

"""
    _near_constant_colorrange(colorvals::Vector{Float64}) -> Union{Nothing,Tuple{Float64,Float64}}

Returns `nothing` if `colorvals`' spread is meaningful (safe to let Makie
auto-scale the colorbar to it, as before); otherwise returns a `colorrange`
tuple wide enough to stop Makie from stretching pure noise across the full
`colormap`.

Motivated directly by the Phase 6 `02b` investigation: `color_by`'s range
gets compared against `hi`/`lo`'s own magnitude (not an absolute
threshold), since "meaningful range" is scale-relative -- a surface
sitting at `z ~ 1e6` with a `1e-1`-wide range is fine, while one at
`z ~ 1.0` with a `6e-8`-wide range (exactly Float32 machine epsilon,
confirmed via `radial_fn`'s printed min/max on the unit sphere: literally
`0.99999994` to `1.0`) is not. When the ratio falls below
`_COLOR_BY_MIN_REL_RANGE`, emits a one-shot `@warn` (the same one-shot
`Ref{Bool}`-latch pattern as [`plot_surface_decomposition`](@ref)'s
`faces`-method winding warning) and returns a `colorrange` centered on the
data's mean, widened to `_COLOR_BY_MIN_REL_RANGE * scale` -- enough to
render as an (expected, correct) near-uniform flat color instead of
speckle, without silently lying about the data by clamping it to look
falsely varied.
"""
function _near_constant_colorrange(colorvals::Vector{Float64})
    lo, hi = extrema(colorvals)
    scale = max(abs(lo), abs(hi), 1.0)
    rel_range = (hi - lo) / scale
    rel_range >= _COLOR_BY_MIN_REL_RANGE && return nothing

    if !_NEAR_CONSTANT_COLOR_WARNED[]
        @warn "plot_surface_decomposition: color_by's range ($(lo) to $(hi), relative spread " *
              "$(rel_range)) is below the $(_COLOR_BY_MIN_REL_RANGE) noise-floor threshold -- " *
              "treating it as numerically near-constant and using a fixed colorrange instead of " *
              "Makie's auto-scaling, which would otherwise stretch floating-point round-off across " *
              "the full colormap and render as misleading speckle. If color_by is expected to vary " *
              "meaningfully here, this indicates a genuine upstream issue (e.g. a degenerate " *
              "surface), not a plotting bug."
        _NEAR_CONSTANT_COLOR_WARNED[] = true
    end
    mid = (lo + hi) / 2
    half_width = scale * _COLOR_BY_MIN_REL_RANGE / 2
    return (mid - half_width, mid + half_width)
end

"""
    plot_surface_decomposition(
        mesh::GeometryBasics.Mesh;
        color_by::Union{Symbol,Function} = :z,
        colormap = :viridis,
        show_wireframe::Bool = false,
        show_colorbar::Bool = true,
        vertices::Union{Nothing,Vector{<:NativeVertex}} = nothing,
        cfg::Union{Nothing,HomotopyConfig} = nothing,
    ) -> Makie.Figure

PRIMARY 3D surface plotting method: takes the ALREADY-WELDED,
ALREADY-correctly-wound `mesh` [`SurfaceDecomposition.decompose_3d_surface`](@ref)
returns (via [`SurfaceDecomposition.weld_mesh`](@ref)) -- not `faces`
alone. See this file's header for why the winding correction makes
`mesh` the right primary input; the `faces::Vector{Face{T}}` method
below is the explicitly-flagged secondary path.

- `color_by`: `:x`/`:y`/`:z` extract that coordinate per mesh point
  (`:z`, the default, matches the roadmap's own sketch); a `Function`
  `(x,y,z) -> Real` generalizes this to BertiniReal's `BY_FUNCTION`
  (e.g. distance from origin, or any other diagnostic of interest) at
  negligible extra cost.
- `show_colorbar`: adds a `Colorbar` for the (continuous) `color_by`
  scale.
- `show_wireframe`: overlays the mesh's own triangle edges.
- `vertices`: optional overlay (e.g. the surface's own critical/singular
  points from `decompose_3d_surface`'s first return value) via the SAME
  [`_plot_vertices_by_type!`](@ref) helper `plot_curve_decomposition`
  uses -- the 3D analogue of BertiniReal overlaying `critical_curve`/
  `singular_curves` on a raw surface.
- `cfg`, if given, sets axis limits from `cfg.bbox_x`/`cfg.bbox_y`/
  `cfg.bbox_z`, mirroring `plot_curve_decomposition`'s own convention.

If `color_by`'s range across `mesh` is near-constant relative to its own
magnitude (see [`_near_constant_colorrange`](@ref)), a fixed `colorrange`
is used instead of Makie's default auto-scaling, and a one-shot `@warn` is
emitted -- prevents floating-point round-off on an otherwise-correct mesh
from rendering as misleading full-spectrum speckle (see this file's own
Phase 6 follow-up investigation for the motivating case).

Returns the `Figure` without calling `display` (see
[`plot_curve_decomposition`](@ref)'s docstring for why).
"""
function plot_surface_decomposition(
    mesh::GeometryBasics.Mesh;
    color_by::Union{Symbol,Function} = :z,
    colormap = :viridis,
    show_wireframe::Bool = false,
    show_colorbar::Bool = true,
    vertices::Union{Nothing,Vector{<:NativeVertex}} = nothing,
    cfg::Union{Nothing,HomotopyConfig} = nothing,
)
    pts = GeometryBasics.coordinates(mesh)

    colorvals = if color_by isa Function
        [Float64(color_by(p[1], p[2], p[3])) for p in pts]
    elseif color_by == :x
        [Float64(p[1]) for p in pts]
    elseif color_by == :y
        [Float64(p[2]) for p in pts]
    elseif color_by == :z
        [Float64(p[3]) for p in pts]
    else
        throw(ArgumentError(
            "plot_surface_decomposition: color_by must be :x, :y, :z, or a Function(x,y,z)->Real, got $(color_by)",
        ))
    end

    fig = Figure()
    ax = Axis3(fig[1, 1]; aspect = :data, xlabel = "x", ylabel = "y", zlabel = "z")

    if cfg !== nothing
        xlims!(ax, cfg.bbox_x...)
        ylims!(ax, cfg.bbox_y...)
        zlims!(ax, cfg.bbox_z...)
    end

    fixed_range = _near_constant_colorrange(colorvals)
    mp = if fixed_range === nothing
        mesh!(ax, mesh; color = colorvals, colormap = colormap)
    else
        mesh!(ax, mesh; color = colorvals, colormap = colormap, colorrange = fixed_range)
    end
    show_colorbar && Colorbar(fig[1, 2], mp)
    show_wireframe && wireframe!(ax, mesh; color = :black, linewidth = 0.5)

    if vertices !== nothing
        handles = _plot_vertices_by_type!(ax, vertices)
        isempty(handles) || axislegend(ax)
    end

    return fig
end

const _FACES_WINDING_WARNED = Ref(false)

"""
    plot_surface_decomposition(
        faces::Vector{Face{T}};
        colormap = :viridis,
        show_wireframe::Bool = false,
        cfg::Union{Nothing,HomotopyConfig} = nothing,
    ) where {T<:AbstractFloat} -> Makie.Figure

SECONDARY 3D surface plotting method (multiple dispatch on the same
exported name, not a keyword flag -- see this file's header): BertiniReal's
"color by cell" view, one solid, distinct color per `Face` (cycled from
`colormap` via `Makie.cgrad(colormap, n; categorical = true)`, exactly
like [`plot_curve_decomposition`](@ref)'s `edge_color_by = :cell`
default), useful for visually distinguishing which triangles belong to
which topological cell -- information [`weld_mesh`](@ref)'s globally-welded
output has already discarded.

**Known, deliberately-flagged limitation**: unlike the `mesh` method
above, this path builds each face's local mesh directly from
`face.mesh_vertices`/`face.mesh_topology` and does NOT reproduce
`weld_mesh`'s gradient-based winding correction (that correction needs
`patch`/`F`, which this method never receives -- see this file's
header). Some faces may therefore render with inverted/inward-facing
normals under directional lighting. This is flagged via an explicit,
one-time `@warn` on first call (not merely a docstring note, and not a
per-call warning either -- a module-level `Ref{Bool}` latch), since a
visually-inverted face could otherwise look like a rendering bug rather
than a known, documented limitation of this specific path.
"""
function plot_surface_decomposition(
    faces::Vector{Face{T}};
    colormap = :viridis,
    show_wireframe::Bool = false,
    cfg::Union{Nothing,HomotopyConfig} = nothing,
) where {T<:AbstractFloat}
    if !_FACES_WINDING_WARNED[]
        @warn "plot_surface_decomposition(::Vector{Face}) does not apply weld_mesh's " *
              "gradient-based winding correction -- some faces may render with inverted/" *
              "inward-facing normals under directional lighting. This is a known, documented " *
              "limitation of this per-cell coloring path (see weld_mesh's own docstring), not a " *
              "rendering bug. Prefer plot_surface_decomposition(mesh::GeometryBasics.Mesh) (the " *
              "output of decompose_3d_surface) for correctly-oriented shading."
        _FACES_WINDING_WARNED[] = true
    end

    fig = Figure()
    ax = Axis3(fig[1, 1]; aspect = :data, xlabel = "x", ylabel = "y", zlabel = "z")

    if cfg !== nothing
        xlims!(ax, cfg.bbox_x...)
        ylims!(ax, cfg.bbox_y...)
        zlims!(ax, cfg.bbox_z...)
    end

    n_faces = length(faces)
    face_colors = Makie.cgrad(colormap, max(n_faces, 2); categorical = true)

    for (i, face) in enumerate(faces)
        size(face.mesh_vertices, 1) == 0 && continue
        points3 = [
            GeometryBasics.Point3f(face.mesh_vertices[r, 1], face.mesh_vertices[r, 2], face.mesh_vertices[r, 3])
            for r in 1:size(face.mesh_vertices, 1)
        ]
        tris = [
            GeometryBasics.TriangleFace{Int}(face.mesh_topology[r, 1], face.mesh_topology[r, 2], face.mesh_topology[r, 3])
            for r in 1:size(face.mesh_topology, 1)
        ]
        face_mesh = GeometryBasics.Mesh(points3, tris)

        mesh!(ax, face_mesh; color = face_colors[i])
        show_wireframe && wireframe!(ax, face_mesh; color = :black, linewidth = 0.5)
    end

    return fig
end

"""
    interactive_3d_viewer(
        mesh::GeometryBasics.Mesh;
        color_by::Union{Symbol,Function} = :z,
        colormap = :viridis,
        show_wireframe::Bool = false,
        show_colorbar::Bool = true,
        vertices::Union{Nothing,Vector{<:NativeVertex}} = nothing,
        cfg::Union{Nothing,HomotopyConfig} = nothing,
    ) -> Makie.Figure

Thin wrapper around [`plot_surface_decomposition`](@ref)'s `mesh` method
(reused internally, not duplicated) that activates the `GLMakie` backend
and `display`s the resulting `Figure` before returning it -- unlike
[`plot_surface_decomposition`](@ref)/[`plot_curve_decomposition`](@ref),
whose entire contract is composability (return a `Figure`, never force a
window), THIS function's entire contract is "open an interactive
window", so forcing `display` here is a deliberate, documented asymmetry.

Rotation/zoom come for free from `GLMakie`'s default `Axis3` camera --
no custom trackball code needed (unlike BertiniReal's separate,
hand-written `glumpy`/OpenGL viewer, which duplicates its own matplotlib
renderer's data-extraction logic just to get free camera controls).

A live/interactive re-slicing control (dragging a z-value to recompute a
NEW cross-section on the fly) is explicitly OUT of scope here -- see this
file's header; BertiniReal's own "slice exploration" only toggles
visibility of already-computed curves, never recomputes one.
"""
function interactive_3d_viewer(
    mesh::GeometryBasics.Mesh;
    color_by::Union{Symbol,Function} = :z,
    colormap = :viridis,
    show_wireframe::Bool = false,
    show_colorbar::Bool = true,
    vertices::Union{Nothing,Vector{<:NativeVertex}} = nothing,
    cfg::Union{Nothing,HomotopyConfig} = nothing,
)
    GLMakie.activate!()
    fig = plot_surface_decomposition(
        mesh;
        color_by = color_by,
        colormap = colormap,
        show_wireframe = show_wireframe,
        show_colorbar = show_colorbar,
        vertices = vertices,
        cfg = cfg,
    )
    display(fig)
    return fig
end
