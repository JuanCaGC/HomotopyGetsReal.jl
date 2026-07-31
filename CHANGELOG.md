# Changelog

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
  `decompose_3d_surface(...; deflate = true)`; diagnostic-only.
- Generic projection support for `decompose_3d_surface` via a new
  `projection` keyword (`:random`, a user-supplied orthogonal matrix, or
  `nothing`), plus `random_orthogonal_matrix`.
- Face/edge incidence tracking (`incidence = true`) and cross-call
  persistent vertex identity via `VertexRegistry`/`register!`, reducing
  mesh naked-edge counts through boundary snap-unification and coordinated
  triangle lofting in `weld_mesh`.
- `CritSlice`, `ColumnLanding`, `SurfaceIncidence` exported types.

### Fixed
- `sample_edge` no longer approximates curved edges with a straight chord.
- Several docstring-rendering/cross-reference issues from two internal
  documentation audits.
- A flaky winding-check assertion at genuine surface singularities.

### Changed
- `Combinatorics`/`Random` now properly declared with `compat` bounds.
- `HomotopyConfig` gained 5 new defaulted fields (all keyword-constructed;
  no existing call sites affected).
