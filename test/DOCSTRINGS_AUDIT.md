# Docstrings audit — exported API

Automated check: every symbol in `src/HomotopyGetsReal.jl` `export` list has a docstring on its definition (or parent type for enum members).

| Export | Docstring | Location |
|--------|-----------|----------|
| `HomotopyConfig` | HAS | `src/Config.jl` |
| `VertexType` | HAS | `src/Types.jl` (documents `Critical`, `Boundary`, `Singular`, `Artificial`) |
| `NativeVertex` | HAS | `src/Types.jl` |
| `Edge` | HAS | `src/Types.jl` |
| `Face` | HAS | `src/Types.jl` |
| `cluster_vertices` | HAS | `src/Clustering.jl` |
| `cluster_scalars` | HAS | `src/Clustering.jl` |
| `cluster_points_indexed` | HAS | `src/Clustering.jl` |
| `jacobian_rank_info` | HAS | `src/Solver.jl` |
| `compute_critical_points` | HAS | `src/Solver.jl` |
| `intersect_bounding_object` | HAS | `src/Solver.jl` |
| `is_near_singular` | HAS | `src/PathTracking.jl` |
| `build_tracker` | HAS | `src/PathTracking.jl` |
| `track_path` | HAS | `src/PathTracking.jl` |
| `track_bidirectional` | HAS | `src/PathTracking.jl` |
| `compute_midslice` | HAS | `src/Topology.jl` |
| `connect_the_dots!` | HAS | `src/Topology.jl` |
| `sample_edge` | HAS | `src/Topology.jl` |
| `decompose_1d_curve` | HAS | `src/Topology.jl` |
| `build_patch_system` | HAS | `src/FaceTracking.jl` |
| `patch_direction` | HAS | `src/FaceTracking.jl` |
| `build_face_tracker` | HAS | `src/FaceTracking.jl` |
| `track_dense_path` | HAS | `src/FaceTracking.jl` |
| `track_face` | HAS | `src/FaceTracking.jl` |
| `sweep_face_bidirectional` | HAS | `src/FaceTracking.jl` |
| `compute_critical_z_slices` | HAS | `src/SurfaceDecomposition.jl` |
| `slice_at_z` | HAS | `src/SurfaceDecomposition.jl` |
| `decompose_3d_surface` | HAS | `src/SurfaceDecomposition.jl` |
| `weld_mesh` | HAS | `src/SurfaceDecomposition.jl` |
| `plot_curve_decomposition` | HAS | `src/Visuals.jl` |
| `plot_surface_decomposition` | HAS | `src/Visuals.jl` (both methods) |
| `interactive_3d_viewer` | HAS | `src/Visuals.jl` |

**Result: 0 missing** (enum value exports inherit `VertexType` docstring).
