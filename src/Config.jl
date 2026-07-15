# src/Config.jl
#
# Phase 1: central, precision-parametric configuration object for the
# HomotopyGetsReal pipeline. Every numerical tolerance used by the six
# algorithmic steps (compute_critical_points, intersect_bounding_object,
# interslice/MidSlice!, ConnectTheDots!, Merge/GetMergeCandidates,
# sample_edge) should be threaded through a `HomotopyConfig{T}` instance
# instead of being hard-coded, so that the whole pipeline can be re-run
# in BigFloat precision without touching algorithm code.

using Parameters

"""
    HomotopyConfig{T<:AbstractFloat}

Central configuration for numerical tolerances and algorithm knobs.

All pipeline steps read settings from a single `HomotopyConfig{T}` so
they can be tuned in one place and instantiated at a chosen floating-point
precision `T` (e.g. `Float64` or `BigFloat`). Numeric defaults use `T(...)`
so `HomotopyConfig{BigFloat}()` gets full-precision literals.

# Fields

## Tolerances

- `critical_point_tol`: Residual tolerance for accepting a homotopy endpoint as a critical point.
- `boundary_point_tol`: Distance tolerance for classifying points on the bounding-box boundary.
- `vertex_match_tol`: Distance tolerance for merging duplicate vertex coordinates.
- `jacobian_rank_tol`: Singular-value cutoff when computing numerical Jacobian rank.
- `singular_value_threshold`: Threshold for classifying a vertex as `Singular`.
- `path_tracker_precision`: Requested precision for HomotopyContinuation path tracking.
- `patch_transversality_cos_tol`: Cosine threshold for re-anchoring skewed z-sweep patches.

## Knobs

- `max_path_steps`: Maximum path-tracker steps before declaring failure.
- `bbox_x`, `bbox_y`, `bbox_z`: Bounding-box `(min, max)` extents along each axis.
- `edge_sample_density`: Number of points sampled per edge in `sample_edge`.
- `midslice_sample_density`: Sample count when building a mid-slice mesh.
- `z_mid_retry_frac`: Fraction of slab width used to perturb `z_mid` on retry.
- `max_z_mid_retries`: Maximum perturbed-`z_mid` attempts in `_robust_slice_at_z`.
- `z_mid_gradient_ratio_tol`: Minimum patch-direction strength ratio for accepting a slice.
"""
@with_kw struct HomotopyConfig{T<:AbstractFloat}
    # critical_point_tol: solution quality at compute_critical_points ("did the solver converge?").
    # boundary_point_tol: geometric containment at intersect_bounding_object ("is this on the box?").
    # vertex_match_tol: deduplication at ConnectTheDots!/Merge ("are these the same vertex?").
    # Mixing these three is a common source of subtle bugs — each bounds "close enough" at a
    # different pipeline stage (residual vs containment vs identity).
    critical_point_tol::T = T(1e-6)
    boundary_point_tol::T = T(1e-5)
    vertex_match_tol::T = T(1e-4)
    jacobian_rank_tol::T = T(1e-8)
    singular_value_threshold::T = T(1e-6)
    path_tracker_precision::T = T(1e-10)
    # patch_transversality_cos_tol: dimensionless cosine bound for FaceTracking._sweep_direction
    # re-anchoring. Measures rotation of the patch anchor gradient vs the local surface gradient
    # (in-slice x,y components), not absolute Jacobian singular values — reusing
    # singular_value_threshold for this fired one hop too late on asymmetric-ellipsoid regression
    # because that scale tracks equation magnitude, not patch skew. Default 0.9 ≈ 26° rotation.
    patch_transversality_cos_tol::T = T(0.9)
    max_path_steps::Int = 1000
    bbox_x::Tuple{T,T} = (T(-4.0), T(4.0))
    bbox_y::Tuple{T,T} = (T(-4.0), T(4.0))
    bbox_z::Tuple{T,T} = (T(-4.0), T(4.0))
    edge_sample_density::Int = 50
    midslice_sample_density::Int = 100
    # z_mid_retry_frac: perturbation step as a fraction of slab width (z_top - z_bottom), not an
    # absolute offset. Used by SurfaceDecomposition._robust_slice_at_z when a slab's naive z_mid
    # yields both an Artificial :endpoint_fallback vertex and a Singular vertex (e.g. Taubin heart
    # at z_mid=0). Default 0.01 (1%) is empirically generous for that case.
    z_mid_retry_frac::T = T(0.01)
    # max_z_mid_retries: attempt cap for perturbed-z_mid retries before _robust_slice_at_z throws.
    max_z_mid_retries::Int = 8
    # z_mid_gradient_ratio_tol: second gate in _robust_slice_at_z — minimum per-edge patch_direction
    # magnitude at candidate z_mid vs cross-z reference at slab quarter-points. Catches topologically
    # clean but numerically ill-conditioned slices (Taubin heart [-1,1] slab). Dimensionless ratio;
    # raw residual and |∇_xy f|/|f_z| thresholds were tried and rejected (see _robust_slice_at_z).
    z_mid_gradient_ratio_tol::T = T(0.01)
end
