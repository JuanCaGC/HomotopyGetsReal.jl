# Capability survey — done

All 24 fixtures run sequentially, one Julia subprocess per fixture, via
`julia --project=dev/scratch/capability_survey run_fixture.jl <name>`, each
wrapped in the agent's own Bash-tool `timeout: 600000` (600s cap per
fixture). No fixture hit the timeout; every fixture completed (or
exceptioned, or ran empty) well inside the cap — longest single fixture was
`torus` at 40.1s.

## `timeout`/`gtimeout`

Confirmed at the start of this run: neither `timeout` nor `gtimeout`
resolves on this machine (`which timeout` / `which gtimeout` both exit 1).
No attempt was made to install one. The Bash tool's own `timeout` parameter
was used on every per-fixture Julia invocation instead, exactly as
instructed — this worked cleanly throughout; no fixture needed it to
actually fire (none hung), so it was never exercised as a live kill switch
in this run, but it was in place for all 24 invocations.

## Per-fixture outcomes

| # | fixture | category | outcome | wall time (s) |
|---|---|---|---|---|
| 1 | circle | curve | clean_success | 14.8 |
| 2 | astroid | curve | clean_success | 15.7 |
| 3 | nodal_cubic | curve | clean_success | 17.3 |
| 4 | cuspidal_cubic | curve | clean_success | 16.9 |
| 5 | squircle_quartic | curve | success_with_caveats | 15.4 |
| 6 | lemniscate_bernoulli | curve | clean_success | 15.6 |
| 7 | folium_descartes | curve | success_with_caveats | 15.9 |
| 8 | cissoid_diocles | curve | clean_success | 17.7 |
| 9 | tangent_parabolas_reducible | curve | clean_success | 17.1 |
| 10 | three_concurrent_lines_reducible | curve | success_with_caveats | 18.4 |
| 21 | elliptic_two_components | curve | clean_success | 16.9 |
| 23 | empty_curve | curve | empty_real_locus | 13.6 |
| 11 | sphere | surface | clean_success | 23.0 |
| 12 | ellipsoid | surface | success_with_caveats | 23.6 |
| 13 | taubin_heart | surface | clean_success | 34.2 |
| 14 | torus | surface | clean_success | 40.1 |
| 15 | whitney_umbrella | surface | exception | 15.9 |
| 20 | horn_torus | surface | success_with_caveats | 16.7 |
| 22 | cone | surface | success_with_caveats | 16.8 |
| 16 | hyperboloid_one_sheet | surface | clean_success | 23.0 |
| 17 | hyperboloid_two_sheets_reducible | surface | clean_success | 23.7 |
| 18 | elliptic_paraboloid | surface | clean_success | 23.8 |
| 19 | quartic_superellipsoid | surface | success_with_caveats | 24.1 |
| 24 | empty_surface | surface | empty_real_locus | 15.6 |

**Tally:** 14 `clean_success`, 7 `success_with_caveats`, 1 `exception`
(whitney_umbrella, confirming an already-documented prediction), 2
`empty_real_locus` (empty_curve, empty_surface, both confirming the
predicted graceful-empty behavior), 0 `timeout`.

Sum of `wall_time_seconds` across all 24 fixtures: **475.8s (~7.9 min)**.
Total session wall time (satellite-env setup through final log write,
including ~5 diagnostic re-runs beyond the 24 official fixtures):
**~43 minutes** (directory created 03:24:48, last file written 04:07:23).

## Known predictions checked against live behavior

- **Whitney umbrella (15):** CONFIRMED exactly. Bare call throws
  `OverflowError: Cannot compute a start system`, via
  `intersect_bounding_object` -> HC.jl `polyhedral` start-system
  construction, reached through `_robust_slice_at_z` -> `slice_at_z` ->
  `decompose_1d_curve` at the degenerate naive `z=0` midpoint. Not a new
  finding; logged as confirmation.
- **Cone (22):** CONFIRMED — first live check of this previously
  "untested" backlog entry. Bare call returns a completely empty
  decomposition (0 vertices/edges/faces/mesh), no exception.
- **Astroid (2):** this run found all 4 cusps cleanly (`Singular`, 0
  `Artificial`) — one of the documented "good" ~3/5 outcomes, not a
  contradiction of the ~2/5-miss prediction (single run, not a 5-run
  sample).
- **Torus (14):** ran directly with `projection=:random, rng=Xoshiro(42)`
  per instructions (z-aligned default skipped, not re-discovered). Result
  matches the documented fix exactly: 8 vertices all `Critical` (0
  `Singular`), 8 edges, 8 faces.
- **Nodal cubic (3)** and the multiplicity>=2 pattern generally: mixed.
  `nodal_cubic` itself found its node correctly as `Singular` (a positive
  result, notably different from the pure node-curve's documented
  zero-solutions failure). But the pattern showed up newly, live, in three
  fixtures the docs hadn't tried: `squircle_quartic` and
  `quartic_superellipsoid` (smooth points misclassified `Singular` due to
  `x^4`'s higher-order-vanishing derivative, edges/mesh wired to bogus
  off-curve `Artificial` placeholders instead) and
  `three_concurrent_lines_reducible` (triple point misclassified
  `Artificial`, plus an entirely separate finding: the vertical line
  component is structurally invisible to x-slicing). `folium_descartes`
  additionally surfaced run-to-run nondeterminism (1 of 3 identical runs
  silently dropped a curve segment).

## New findings this survey (see `errors_log.md` for full detail/evidence)

1. `squircle_quartic` — genuine smooth critical points at `(+-1,0)`
   misclassified `Singular`; both edges wired to off-curve `Artificial`
   placeholders; render visually confirms most of the curve is missing,
   not merely mislabeled.
2. `folium_descartes` — cross-process nondeterminism: 1 of 3 identical
   runs silently dropped a whole curve segment (orphaned boundary vertex),
   no exception.
3. `three_concurrent_lines_reducible` — the vertical line component
   (`x=0`) is structurally invisible to x-parametrized slicing (never
   appears in any edge); triple point misclassified `Artificial` not
   `Singular`. Visually confirmed in the render (only the X-shaped
   diagonal pair is drawn).
4. `ellipsoid` — max/p99 residual ~170x above p90, a small number of
   outlier mesh points, plausibly connected to the documented ellipsoid
   pole-tracking discovery.
5. `horn_torus` — extends the cone's silent-empty-decomposition pattern to
   a second, previously-untried fixture (self-tangent pinch point at the
   origin, symmetric `bbox_z`).
6. `quartic_superellipsoid` — surface-level recurrence of finding 1,
   worse: `Artificial` edge-endpoint vertices at residual ~1.0 (nowhere
   near the surface); render visually confirms a real wedge-shaped
   notch/gouge in the mesh near one pole, not merely a bookkeeping-only
   issue (this correction was made after visually inspecting the render,
   contradicting an initial from-the-numbers-alone assessment).
7. `plot_surface_decomposition` throws `ArgumentError` (uncaught, from
   `_near_constant_colorrange`'s `extrema` call) on any completely empty
   mesh, rather than handling it gracefully — hit on `horn_torus`, `cone`,
   `empty_surface` (3 occurrences, same root cause, logged once in detail).

## Positive findings (HGR beats or meets expectation)

- `elliptic_two_components`: both disconnected real components (bounded
  oval + unbounded branch) correctly found and separated — the specific
  thing this fixture exists to check.
- `hyperboloid_two_sheets_reducible`: both disconnected sheets correctly
  found, counts exactly double the one-sheet analog.
- `tangent_parabolas_reducible`: the tacnode (higher-order tangency,
  effectively multiplicity 4) correctly resolved — all 4 boundary exits,
  singular origin found, fully correct topology.
- `nodal_cubic`: node correctly found and classified `Singular`, contrary
  to the pure node-curve's documented zero-solutions failure.
- `cissoid_diocles`, `circle`: residuals below the sphere baseline.

## Self-check: `git status --porcelain`

```
?? dev/scratch/capability_survey/
```

Exactly one line, referring only to this survey's own directory. No file
under `src/`, `test/`, `docs/`, `paper_artifacts/`, the repo-root
`Project.toml`/`Manifest.toml`, or any other tracked path was written,
edited, or staged. No git write operation (`add`/`commit`/`checkout`/
`reset`) was run at any point in this survey.

(Note: the conversation's initial git-status snapshot additionally showed
`M .gitignore` and `?? dev/scratch/renders_torus_validation/`. Investigated
and confirmed unrelated to this survey: both were resolved by a separate
commit, `0f9ee57` / `4c6dfba`, made at 02:17-03:02, i.e. *before* this
survey's directory was even created at 03:24:48. Not caused by, or touched
during, this work.)
