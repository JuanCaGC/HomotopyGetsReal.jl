# Assertion Audit — Phase 7 Test Conversion

Baseline counts from scratch scripts (method: `^\s*@test[^_a-zA-Z]`, `@test_throws`, `= @inferred `).

## Taubin runtime (measured in this environment, 2026-07-07)

| Measurement | Wall time |
|-------------|-----------|
| Full `scratch_phase5_taubin_check.jl` | **77.4 s** |
| `compute_critical_z_slices` (Taubin heart, low-res cfg) | **18.52 s** |
| `decompose_3d_surface` (Taubin heart, low-res cfg) | **27.54 s** |
| Core two calls combined (no section-6 extras) | **46.06 s** |

**→ FLAG FOR MANUAL REVIEW:** skip-gate design (env-var vs CI tag). Measured numbers above; implementation uses `HOMOTOPYGETSREAL_RUN_SLOW_TESTS=1` gate pending your decision.

---

## Scratch → test file mapping

| Scratch source | Test file | `@test` | `@test_throws` | `@inferred` | Expected total |
|----------------|-----------|---------|----------------|-------------|----------------|
| `scratch_phase1_check.jl` | `test_types.jl` | 18 | 0 | 5 | **23** |
| `scratch_phase2_check.jl` | `test_solver.jl` | 19 | 0 | 3 | **22** |
| `scratch_phase3_check.jl` | `test_topology.jl` | 33 | 0 | 5 | **38** |
| `scratch_phase4_check.jl` | `test_pathtracking.jl` | 16 | 0 | 3 | **19** |
| `scratch_phase5_check.jl` | `test_surfacedecomposition.jl` | 50 | 0 | 6 | **56** |
| `scratch_phase5_taubin_check.jl` | `test_taubin.jl` | 13 | 0 | 0 | **13** |
| `scratch_phase6_check.jl` | `test_visuals.jl` | 17 | 1 | 0 | **18** |
| `scratch_phase6_investigate.jl` | `test_visuals.jl` (new) | 0 → **5** | 0 | 0 | **+5** |
| | **Grand total** | **166 + 5** | **1** | **22** | **194** |

---

## New tests from `scratch_phase6_investigate.jl` → `test_visuals.jl`

These are **named by assertion**, not generic placeholders:

| # | `@testset` name | Assertion (exact) | Captures |
|---|-----------------|-------------------|----------|
| N1 | `investigation regressions / radial_fn relative range` | `(maximum(radvals) - minimum(radvals)) / max(abs(maximum(radvals)), abs(minimum(radvals)), 1.0) < 1e-4` on unit-sphere welded mesh | Near-constant `color_by` range (Float32 round-off ~6e-8), motivates `_near_constant_colorrange` |
| N2 | `investigation regressions / welded mesh outward normals` | `n_inward == 0` after per-triangle `dot(normal, radial) > 0` check on all `GeometryBasics.faces(mesh)` | Corrected mesh path has no inverted triangles |
| N3 | `investigation regressions / interactive_3d_viewer dispatch` | `length(methods(interactive_3d_viewer)) == 1` AND sole method's first parameter type is `GeometryBasics.Mesh` | `interactive_3d_viewer` calls mesh path, not faces path |
| N4 | `investigation regressions / interactive_3d_viewer reuses mesh plotter` | `interactive_3d_viewer(mesh; ...)` returns `Makie.Figure` with same type as `plot_surface_decomposition(mesh; ...)` (both succeed; optional: same mesh vertex count in figure data — type check only to avoid brittle pixel compare) | No accidental wrong-method dispatch |
| N5 | `investigation regressions / near-constant colorrange helper` | `_near_constant_colorrange(fill(1.0, 10)) !== nothing` (returns fixed range tuple) | Helper fires on genuinely near-constant data (unit test of guard itself) |

Note: N5 tests `_near_constant_colorrange` directly; radial_fn warning latch is already covered by scratch_phase6_check assertions (not duplicated).

---

## Per-file completion status

| test_file | expected | actual @test | actual @test_throws | actual @inferred | total (inventory) | Test.jl Pass (fast) | status |
|-----------|----------|--------------|---------------------|------------------|-------------------|---------------------|--------|
| `test_types.jl` | 23 | 18 | 0 | 5 | 23 | 23 | PASS |
| `test_solver.jl` | 22 | 19 | 0 | 3 | 22 | 19 | PASS* |
| `test_topology.jl` | 38 | 33 | 0 | 5 | 38 | 33 | PASS* |
| `test_pathtracking.jl` | 19 | 16 | 0 | 3 | 19 | 16 | PASS* |
| `test_surfacedecomposition.jl` | 56 | 50 | 0 | 6 | 56 | 50 | PASS* |
| `test_taubin.jl` | 13 | 13 | 0 | 0 | 13 | 13 (slow only) | PASS |
| `test_visuals.jl` | 23 (18 scratch + 5 investigate) | 23 | 1 | 0 | 24 | 21 fast / 24 slow | PASS |

\* `Test.jl` summary counts 162 Pass in fast mode. Inventory totals 179 executable assertions in fast mode (including assignment-form `@inferred`); the 17 assignment-form `@inferred` forms run and would fail on regression but do not increment the Pass counter in Julia 1.12's summary (observed behavior, not lost coverage).

**Fast suite observed:** 162 Pass, ~70s wall time.

**Full suite observed (slow):** 178 Pass, ~135s wall time (`HOMOTOPYGETSREAL_RUN_SLOW_TESTS=1`).

**Taubin heart visuals (3 `@test` in `test_visuals.jl`)** gated behind the same `HOMOTOPYGETSREAL_RUN_SLOW_TESTS=1` as `test_taubin.jl` — not a silent drop; they run only in full/slow mode.

### Taubin skip-gate — FLAG FOR MANUAL REVIEW

Measured in this environment (2026-07-07):

| Benchmark | Time |
|-----------|------|
| Full `scratch_phase5_taubin_check.jl` | **77.4 s** |
| `compute_critical_z_slices` only | **18.52 s** |
| `decompose_3d_surface` only | **27.54 s** |
| Core calls combined | **46.06 s** |

Implementation uses `HOMOTOPYGETSREAL_RUN_SLOW_TESTS=1` env gate. **Your call** on whether this is the right default vs CI tags.
