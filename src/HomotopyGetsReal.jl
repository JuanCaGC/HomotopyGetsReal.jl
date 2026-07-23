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
#
# Convention (2026-07, isosingular deflation Stage 4a): because this is one
# flat namespace, two `@enum` blocks anywhere in `src/` can collide on a
# member name with no compile error -- Julia's `==` across two different
# enum types just silently returns `false` rather than erroring, so a
# collision doesn't announce itself; it quietly breaks whichever comparison
# assumed the wrong one was in scope. This nearly happened between
# `IsosingularVerdict.Inconclusive` (Solver.jl, Stage 3) and a first draft of
# `ResolveVerdict`'s own `Inconclusive` (Stage 4a) -- caught only by
# re-deriving the docstring's own example by hand, not by any check. Before
# adding any new `@enum` to this module: grep the whole `src/` tree for
# every candidate member name first, and prefix multi-value "verdict"-shaped
# enums with a short tag tied to their own owning function (e.g.
# `ResolveResolved`/`ResolveAmbiguous`/`ResolveExhausted`, not bare
# `Resolved`/`Ambiguous`/`Exhausted`) rather than relying on a word simply
# sounding unlikely to collide -- that's exactly what `Inconclusive` sounded
# like the first time too. Not applied retroactively to
# `IsosingularVerdict`/`ResolveVerdict`'s own already-shipped, exported
# members; applies going forward.
#
# Note (2026-07, isosingular deflation Stage 4c validation): established-normal
# round counts for `resolve_isosingular_dimension`, for future reference if a
# deeper or unexplained round count ever shows up elsewhere and someone needs
# to know whether that's within observed range or new. Every ground-truth case
# through Stage 3/4a resolved in 1-2 rounds. The Taubin heart fixture's full
# Stage 4c validation (individually inspecting all 22 real firings during a
# `decompose_3d_surface(...; deflate=true)` run, not just the ones surviving to
# final output) found 17 of 22 at 1 round, but 5 genuine outliers at 4 rounds
# (`corank_seq=[2,1,1,1,0]`), all at the slice-level critical points x=(±1,0)
# -- the heart's left/right-most z-slice extremes. All 22 still resolved
# cleanly (`Resolved`, zero `Ambiguous`/`Exhausted`/`attempts>=15`), so 4
# rounds at a genuine local degeneracy is confirmed normal behavior, not a
# red flag by itself -- but it's a real, non-spurious deviation from the [1,2]
# range seen everywhere else, worth knowing about before assuming a similar
# count elsewhere is a bug. Separately, every one of those 22 firings resolved
# to `isosingular_dimension=0`, and every `d=0` resolution has `attempts=0` by
# construction (`verify_isosingular_dimension` needs zero hyperplanes when
# `d=0`, so its retry mechanism never runs) -- the soft-flag attempts range of
# [1,4] used during this validation was only ever scoped to `d>=1` cases that
# actually exercise that retry loop; `attempts=0` at `d=0` is not itself
# evidence of anything and should not be treated as a flag in future runs.

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
export ResolveVerdict, Resolved, Ambiguous, Exhausted
export ResolveResult, resolve_isosingular_dimension

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
