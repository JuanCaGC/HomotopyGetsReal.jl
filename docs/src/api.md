```@meta
CurrentModule = HomotopyGetsReal
```

# API reference

Exported symbols from `HomotopyGetsReal`. Docstrings are taken from the package sources (Phase 7 docstring audit); no additional narrative is defined here.

## Configuration and types

```@docs
HomotopyConfig
VertexType
NativeVertex
Edge
Face
```

## Clustering

```@docs
cluster_vertices
cluster_scalars
cluster_points_indexed
VertexRegistry
register!
```

## Critical points and bounding objects

```@docs
jacobian_rank_info
compute_critical_points
intersect_bounding_object
```

## Isosingular deflation

```@docs
estimate_corank
deflate_once
IsosingularVerdict
VerifyResult
verify_isosingular_dimension
ResolveVerdict
ResolveResult
resolve_isosingular_dimension
```

## Path tracking

```@docs
is_near_singular
build_tracker
track_path
track_bidirectional
```

## Curve topology

```@docs
compute_midslice
connect_the_dots!
sample_edge
decompose_1d_curve
```

## Face tracking

```@docs
build_patch_system
patch_direction
build_face_tracker
track_dense_path
track_face
sweep_face_bidirectional
```

## Projection

```@docs
random_orthogonal_matrix
```

## Surface decomposition

```@docs
compute_critical_z_slices
slice_at_z
decompose_3d_surface
weld_mesh
CritSlice
ColumnLanding
SurfaceIncidence
```

## Visualization

```@docs
plot_curve_decomposition
plot_surface_decomposition
interactive_3d_viewer
```
