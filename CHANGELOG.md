# Changelog

## [Unreleased]

## [0.2.1] - 2026-08-08

### Fixed
- `plot_surface_decomposition` no longer throws an uncaught `ArgumentError`
  when `decompose_3d_surface` returns a completely empty mesh (0 points, 0
  triangles — a known outcome for some fixtures, e.g. a surface whose only
  critical z lands exactly on an undetected degenerate point). Now renders
  an empty axes frame (axis limits from `cfg` if given) with a one-shot
  `@warn` instead of crashing.

## [0.2.0]

### Breaking
- `sample_edge` now requires an explicit `F::System` argument
  (`sample_edge(F, edge, cfg)`, was `sample_edge(edge, cfg)`) so it can
  re-project samples onto the true curve via Newton iteration instead of
  interpolating a straight chord between endpoints.

### Added
- Isosingular deflation pipeline for classifying `Singular` vertices at
  non-isolated/higher-multiplicity points: `estimate_corank`, `deflate_once`,
  `verify_isosingular_dimension`, `resolve_isosingular_dimension`, plus the
  `IsosingularVerdict`/`ResolveVerdict` result types. Opt-in via
  `decompose_3d_surface(...; deflate = true)`; diagnostic-only (adds
  `metadata[:isosingular_verdict]` etc. to affected vertices; does not
  change mesh geometry).
- Generic projection support for `decompose_3d_surface` via a new
  `projection` keyword (`:random`, a user-supplied orthogonal matrix, or
  `nothing` for the prior z-aligned behavior), plus `random_orthogonal_matrix`.
- Face/edge incidence tracking (`incidence = true` on `decompose_3d_surface`)
  and cross-call persistent vertex identity via `VertexRegistry`/`register!`,
  substantially reducing mesh naked-edge counts through boundary
  snap-unification and coordinated triangle lofting in `weld_mesh`.
- `CritSlice`, `ColumnLanding`, `SurfaceIncidence` exported types supporting
  the above.

### Fixed
- `sample_edge` no longer approximates curved edges with a straight chord
  (see Breaking, above) — this was a real geometric inaccuracy for edges
  with significant curvature between sampled endpoints.
- Several docstring-rendering and cross-reference corruption issues found
  during two internal documentation audits.
- A flaky winding-check assertion at genuine surface singularities.

### Changed
- `Combinatorics` and `Random` are now properly declared with `compat`
  bounds in `Project.toml` (previously present as implicit/undeclared
  dependencies).
- `HomotopyConfig` gained several new fields with defaults (all
  keyword-constructed; no existing call sites are affected):
  `projection_orthonormality_tol`, `min_slab_width`,
  `incidence_snap_tol_ratio`, `isosingular_verify_retries`,
  `max_deflations`.
