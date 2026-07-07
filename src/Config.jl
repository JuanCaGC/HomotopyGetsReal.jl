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

Central configuration object for the HomotopyGetsReal pipeline. All
tolerances and numerical knobs used across the six algorithmic steps
live here so they can be tuned in one place, and so the entire pipeline
can be instantiated in a different floating-point precision `T` (e.g.
`Float64` for speed, `BigFloat` for certified/high-precision runs)
simply by writing `HomotopyConfig{BigFloat}()`.

All numeric defaults are wrapped in `T(...)` so that constructing
`HomotopyConfig{BigFloat}()` actually gets `BigFloat`-precision defaults
instead of `Float64` literals silently truncating the precision.

# Fields

## Tolerances (the three easy-to-confuse ones)

These three tolerances all bound "how close is close enough", but they
are used at very different stages of the pipeline and mixing them up is
a classic source of subtle bugs:

- `critical_point_tol::T`: Used during `compute_critical_points`. Bounds
  how close the path tracker's endpoint must be to an actual solution of
  the critical-point system (i.e. residual/Jacobian-based stopping
  tolerance for Newton correction) before a numerically-found point is
  accepted as a genuine critical point of the algebraic variety. This is
  about *solution quality* of a single homotopy path.

- `boundary_point_tol::T`: Used during `intersect_bounding_object`.
  Bounds how close a tracked point must be to the bounding box's faces
  (or other bounding object) before it is classified as lying *on* the
  boundary rather than strictly inside/outside it. This is about
  *geometric containment*, not path-tracking residual.

- `vertex_match_tol::T`: Used during `ConnectTheDots!` / `Merge` /
  `GetMergeCandidates`. Bounds how close two independently-computed
  `NativeVertex` coordinates must be to be considered *the same
  geometric vertex* (e.g. a critical point recomputed from two adjacent
  slices, or a vertex found via two different homotopy paths). This is
  about *deduplication/identity*, not solution quality or containment.

In short: `critical_point_tol` asks "did the solver converge?",
`boundary_point_tol` asks "is this point on the box?", and
`vertex_match_tol` asks "are these two points actually the same vertex?".

## Other tolerances

- `jacobian_rank_tol::T`: Threshold below which a singular value of the
  Jacobian is treated as zero when computing its numerical rank (used to
  detect singular points / rank-deficient loci).
- `singular_value_threshold::T`: Threshold used when classifying a
  vertex as `Singular` from the smallest singular value(s) of the
  Jacobian at that point (distinct from `jacobian_rank_tol`, which is
  used purely for rank computation; this one drives vertex
  classification/metadata).
- `path_tracker_precision::T`: Requested precision passed to
  HomotopyContinuation.jl's path tracker (controls adaptive
  step-size/Newton corrector precision along a homotopy path).
- `patch_transversality_cos_tol::T`: Used by
  `FaceTracking._sweep_direction`'s adaptive re-anchoring. Bounds the
  cosine of the angle between a z-sweep patch anchor's ORIGINAL fixed
  gradient direction and the CURRENT local surface gradient (both
  restricted to their in-slice `(x,y)` components) before the patch is
  considered too skewed relative to the curve's drifting tangent and is
  rebuilt at the current point. This is a genuinely different kind of
  quantity from `jacobian_rank_tol`/`singular_value_threshold` above
  (which bound absolute Jacobian singular values of the *system actually
  being solved*, at whatever scale that system's own equations happen to
  carry) -- this one is a dimensionless, scale-free RATIO (a cosine),
  because the quantity of interest here is purely the *rotation* of the
  gradient direction, not any solve residual or absolute singular value.
  Reusing `singular_value_threshold` for this was tried first and found
  empirically unusable: on an asymmetric-ellipsoid regression case, the
  augmented system's smallest singular value stayed well above
  `singular_value_threshold`'s default (`1e-6`) for every sample point
  right up until the hop where transversality had *already* collapsed
  and the tracked point had already landed measurably off-surface (the
  singular value's absolute scale is set by the surface's own gradient
  magnitude at that point, not by how close the patch is to becoming
  tangent to the curve) -- i.e. that check reliably fired one hop too
  late to actually prevent the failure it was meant to catch. The cosine
  threshold instead directly measures the geometric quantity that
  determines patch validity (see `_sweep_direction`'s docstring for the
  full derivation of why `H_sys`'s Jacobian determinant is exactly this
  dot product), so a single default (`0.9`, i.e. re-anchor once the
  gradient has rotated more than ~26° from the anchor) is meaningful
  independent of the surface's scale or `midslice_sample_density`.

## Numerical/algorithmic knobs

- `max_path_steps::Int`: Maximum number of steps the path tracker is
  allowed to take before a path is declared a failure/truncation.
- `bbox_x::Tuple{T,T}`, `bbox_y::Tuple{T,T}`, `bbox_z::Tuple{T,T}`: The
  `(min, max)` extents of the bounding box used by
  `intersect_bounding_object` along each axis.
- `edge_sample_density::Int`: Number of points sampled along each `Edge`
  in `sample_edge`.
- `midslice_sample_density::Int`: Number of sample points used when
  building a mid-slice mesh in `interslice`/`MidSlice!`.
- `z_mid_retry_frac::T`: Used by
  `SurfaceDecomposition._robust_slice_at_z`. When a slab's naive
  `z_mid = (z_bottom+z_top)/2` produces a
  [`Topology.decompose_1d_curve`](@ref) result containing BOTH an
  `Artificial` vertex tagged `metadata[:origin] == :endpoint_fallback`
  (i.e. `Topology._resolve_endpoint` had to fabricate a vertex because a
  tracked path failed to close back onto any known vertex) AND a
  `Singular`-typed vertex (the co-occurrence, not `:endpoint_fallback`
  alone, is what actually distinguishes a genuinely non-reduced/
  repeated-factor curve -- e.g. the Taubin heart surface's `z_mid=0`
  slice -- from a merely topologically-complex-but-well-conditioned one;
  see `_robust_slice_at_z`'s own docstring for the empirical comparison
  that ruled out both "`:endpoint_fallback` alone" and a raw
  residual-magnitude threshold), `_robust_slice_at_z` retries at a
  perturbed `z_mid`. This field is the perturbation step size, expressed
  as a FRACTION OF THE SLAB WIDTH `(z_top - z_bottom)` rather than an
  absolute offset, so the same default behaves sensibly regardless of a
  given surface's own z-extent. A new field, not a reuse of
  `vertex_match_tol`/`critical_point_tol`/etc.: none of those tolerances
  mean "how far to nudge a slicing plane", and reusing one anyway would
  tie an unrelated numerical knob's tuning to this different job (the
  same reasoning that justified `patch_transversality_cos_tol` as its
  own field above). Default `0.01` (1% of slab width) is empirically
  generous: on the Taubin heart's degenerate `z_mid=0` case, even a
  0.4%-of-width perturbation already produced a clean, non-`Artificial`
  decomposition (`scratch_phase5_taubin_check.jl` section 6's retry
  report).
- `max_z_mid_retries::Int`: Companion to `z_mid_retry_frac`, analogous to
  `max_path_steps` (an attempt-count cap, not a tolerance): the number of
  perturbed-`z_mid` attempts `_robust_slice_at_z` makes before giving up
  and throwing a loud `ErrorException` naming the offending slab, rather
  than silently proceeding with (or looping forever trying to escape) a
  slice already known to contain an unexplained endpoint-fallback vertex.
  Default `8`.
- `z_mid_gradient_ratio_tol::T`: Second, independent gate in
  `_robust_slice_at_z`, added after the vertex-type gate above was found
  (Taubin heart `[-1,1]` slab, retry investigation) to accept a
  topologically-clean candidate `z_mid` that was nonetheless still
  numerically ill-conditioned -- `decompose_1d_curve` classified its
  curve cleanly (no stray `Artificial`/`Singular` vertices), but the
  surface's true gradient near that `z_mid` was orders of magnitude
  smaller than elsewhere on the same slab (a genuine near-degenerate
  NEIGHBORHOOD around the excluded exact critical/repeated-factor
  z-value, not just a single pathological point), which wrecked
  `FaceTracking.track_face`'s Newton-based patch tracking downstream
  even though the raw 2D decomposition looked fine. A raw
  residual-magnitude threshold was tried first and rejected (see
  `_robust_slice_at_z`'s own docstring): `sample_edge`'s linear
  interpolation error has no clean scale-free gap to split. A same-point
  `|∇_xy f| / |f_z|` ratio was tried next and also rejected: it does not
  correlate with sweep quality at all (measured ratios for KNOWN-GOOD
  reference slabs, `0.60`-`1.16`, sit BELOW several KNOWN-BAD candidates'
  ratios, `2.08`-`4.29`, making it actively misleading, not merely
  noisy). What DOES work, reusing `FaceTracking.patch_direction` (the
  exact `(a,b)` pair `track_face` itself seeds each sweep with, per this
  field's own investigation) rather than any new gradient-computation
  logic: comparing a candidate `z_mid`'s minimum per-edge anchor
  `hypot(patch_direction(...)...)` against a CROSS-Z reference computed
  at two FIXED, retry-schedule-independent locations well inside the
  same slab (`z_bottom + 0.25*(z_top-z_bottom)` and
  `z_top - 0.25*(z_top-z_bottom)`, i.e. the slab's own quarter-points,
  taking the larger of the two magnitudes found there). This ratio is
  genuinely dimensionless (both numerator and denominator are the same
  `hypot(patch_direction(...)...)` quantity evaluated at different
  `z`-values, so it is invariant to the surface's own gradient scale) and
  showed a clean, wide gap on the Taubin heart: the two already-healthy
  reference slabs measured `0.82`/`0.98` (i.e. ~as strong as the slab's
  own quarter-point baseline), while the `[-1,1]` slab's bad candidates
  measured `0.0014`-`0.0057` and its eventually-accepted good candidate
  measured `0.013`-`0.023` -- nearly two orders of magnitude of margin on
  both sides of the default `0.01` (1%) chosen here, not a threshold
  sitting uncomfortably close to either boundary. A new field, for the
  same reason `z_mid_retry_frac` is: no existing tolerance means "is the
  local gradient strong enough relative to the rest of this slab",
  and reusing one anyway would tie an unrelated knob's tuning to this
  different job. Known limitation, confirmed and left as-is rather than
  tuned around (see `_robust_slice_at_z`'s own docstring): a curve with
  branches of very different local conditioning (e.g. a thin inner loop
  next to a strong outer loop) can trigger one harmless extra retry,
  since the gate compares a candidate's WEAKEST branch against the
  reference's STRONGEST branch.
"""
@with_kw struct HomotopyConfig{T<:AbstractFloat}
    critical_point_tol::T = T(1e-6)
    boundary_point_tol::T = T(1e-5)
    vertex_match_tol::T = T(1e-4)
    jacobian_rank_tol::T = T(1e-8)
    singular_value_threshold::T = T(1e-6)
    path_tracker_precision::T = T(1e-10)
    patch_transversality_cos_tol::T = T(0.9)
    max_path_steps::Int = 1000
    bbox_x::Tuple{T,T} = (T(-4.0), T(4.0))
    bbox_y::Tuple{T,T} = (T(-4.0), T(4.0))
    bbox_z::Tuple{T,T} = (T(-4.0), T(4.0))
    edge_sample_density::Int = 50
    midslice_sample_density::Int = 100
    z_mid_retry_frac::T = T(0.01)
    max_z_mid_retries::Int = 8
    z_mid_gradient_ratio_tol::T = T(0.01)
end
