# src/HomotopyGetsReal.jl
#
# Top-level module for the HomotopyGetsReal.jl rebuild. Following the
# pattern established in Phase 1, this is a single flat module: the
# files below are `include`d directly (none of them declare their own
# `module ... end` block), so `Config.jl`, `Types.jl`, `Clustering.jl`,
# `Solver.jl`, `PathTracking.jl`, `Topology.jl`, `FaceTracking.jl`, and
# `SurfaceDecomposition.jl` all share one namespace. Include order
# matters: `Clustering.jl` must come before `Solver.jl`/`Topology.jl`
# since both call `cluster_vertices`/`cluster_scalars` directly;
# `Solver.jl` must come before `PathTracking.jl`/`Topology.jl` since both
# call `jacobian_rank_info`/`_newton_polish` (Topology) directly;
# `PathTracking.jl` must come before `Topology.jl` since
# `connect_the_dots!` calls `build_tracker`/`track_bidirectional`
# directly; `Topology.jl` must come before `FaceTracking.jl`/
# `SurfaceDecomposition.jl` since `slice_at_z` calls
# `decompose_1d_curve` directly, and `track_dense_path` reuses
# `PathTracking._track_path_segment!` directly; `FaceTracking.jl` must
# come before `SurfaceDecomposition.jl` since `decompose_3d_surface`
# calls `build_patch_system`/`track_face` directly and `weld_mesh` calls
# `_gradient_at` directly; `Projection.jl` must come before
# `SurfaceDecomposition.jl` since `decompose_3d_surface`'s projection
# branch calls `_resolve_projection`/`_rotate_system`/
# `_verify_projection_ok`/`_chart_config`/`_map_to_world` directly;
# `SurfaceDecomposition.jl` must come before
# `Visuals.jl` since `plot_surface_decomposition`'s methods take
# `decompose_3d_surface`'s/`weld_mesh`'s own return types directly (no
# `using` needed within a flat module).

module HomotopyGetsReal

using HomotopyContinuation
using LinearAlgebra
using GenericLinearAlgebra
using Parameters
using GeometryBasics
using GLMakie
using Random
using Combinatorics

include("Config.jl")
include("Types.jl")
include("Clustering.jl")
include("Solver.jl")
include("PathTracking.jl")
include("Topology.jl")
include("FaceTracking.jl")
include("Projection.jl")
include("SurfaceDecomposition.jl")
include("Visuals.jl")

# Config.jl
export HomotopyConfig

# Types.jl
export VertexType, Critical, Boundary, Singular, Artificial
export NativeVertex, Edge, Face

# Clustering.jl
export cluster_vertices, cluster_scalars, cluster_points_indexed
export VertexRegistry, register!

# Solver.jl
export jacobian_rank_info, compute_critical_points, intersect_bounding_object
export estimate_corank, deflate_once
export IsosingularVerdict, Verified, NotTerminal, Inconclusive
export VerifyResult, verify_isosingular_dimension

# PathTracking.jl
export is_near_singular, build_tracker, track_path, track_bidirectional

# Topology.jl
export compute_midslice, connect_the_dots!, sample_edge, decompose_1d_curve

# FaceTracking.jl
export build_patch_system, patch_direction, build_face_tracker, track_dense_path,
       sweep_face_bidirectional, track_face

# Projection.jl
export random_orthogonal_matrix

# SurfaceDecomposition.jl
export compute_critical_z_slices, slice_at_z, decompose_3d_surface, weld_mesh
export CritSlice, ColumnLanding, SurfaceIncidence

# Visuals.jl
export plot_curve_decomposition, plot_surface_decomposition, interactive_3d_viewer

end # module HomotopyGetsReal
