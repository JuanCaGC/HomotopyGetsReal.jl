# Capability Survey — Phase 2 Evaluator Summary

Read-only synthesis of Phase 1's 24-fixture survey of `decompose_1d_curve` /
`decompose_3d_surface` at coarse test-density config
(`edge_sample_density=6`, surfaces also `midslice_sample_density=8`).
Every number below was independently re-read from
`dev/scratch/capability_survey/data/*.json` (all 24 files, parsed directly,
not copied from Phase 1's own summaries), cross-checked against
`docs/DESIGN_NOTES.md`'s actual text (grepped and read in full, not
inferred), and in five cases (`squircle_quartic`, `quartic_superellipsoid`,
`three_concurrent_lines_reducible`, `sphere`, `empty_curve`) against the
actual PNG renders. No decomposition or Julia code was re-run; no file
other than this one was intended to be written.

**Outcome tally, independently counted from the 24 JSON files' own
`outcome` fields:** 14 `clean_success`, 7 `success_with_caveats`, 2
`empty_real_locus`, 1 `exception`. This matches `SURVEYOR_DONE.md`'s
tally exactly — confirmed, not just repeated. Sum of `wall_time_seconds`
across all 24 files, independently recomputed: **475.82s**, matching
`SURVEYOR_DONE.md`'s "475.8s" to the reported precision.

**Render count: 20 of 24.** The 4 missing — `whitney_umbrella` (exception
before the plotting step; no PNG expected), `cone`/`horn_torus`/`empty_surface`
(all three: empty 0/0/0/0 decomposition, `plot_surface_decomposition`
throws `ArgumentError` on the empty mesh via `_near_constant_colorrange`,
`src/Visuals.jl:210`) — are fully accounted for in each fixture's own
`data/<fixture>.json` `notes` field; this was a documentation-completeness
gap (neither this report nor `SURVEYOR_DONE.md` stated the 20/24 figure
explicitly), not a new or unlogged discrepancy. See the 2026-08-06 entry
in `errors_log.md` for the full accounting.

**2026-08-06 addendum**: a targeted follow-up (not a re-run of the
24-fixture survey) investigated visible defects in `torus.png`,
`taubin_heart.png`, and `hyperboloid_one_sheet.png` that the original
survey never flagged as findings. Summary: torus/taubin_heart's visible
gaps are CONFIRMED as the known no-`incidence=true` naked-edge mechanism
(`docs/DESIGN_NOTES.md`'s "Watertightness measurements" entry); the
mechanism reproduces at this survey's own coarse density (torus: 80→8
naked edges, taubin_heart: 132→29, both via `incidence=true`), though
`incidence=true` substantially reduces rather than eliminates the visible
gaps — matching that entry's own "multi-face edge-type boundaries reduce
without fully closing" caveat exactly, not a clean before/after fix.
hyperboloid_one_sheet's jagged rim is a DIFFERENT, unrelated mechanism —
confirmed as a coarse-sampling artifact specific to its bbox-clipped
boundary (smooths out almost completely at default `edge_sample_density=50`
vs. this survey's coarse `6`), not the naked-edge/incidence issue above.
Full detail, numbers, and render paths: `errors_log.md`'s 2026-08-06
entry. §1's table notes and §3 below are updated accordingly.

---

## 1. Results table

Vertex-type shorthand: `C`=Critical, `S`=Singular, `B`=Boundary,
`A`=Artificial (fallback). Omitted letters = 0 of that type. "Gap
targeted" is inferred from each fixture's evident role in the suite, not
a literal JSON field — flagged inline wherever direct verification changed
the obvious reading (see §2.1 for `hyperboloid_two_sheets_reducible`).

### Curves (`decompose_1d_curve`)

| Fixture | Category / gap targeted | Outcome | Wall (s) | Vertices | Edges | Median resid. | Max resid. |
|---|---|---|---|---|---|---|---|
| circle | baseline (smooth closed curve) | clean_success | 14.8 | C2 | 2 | 9.03e-13 | 2.91e-08 |
| astroid | baseline (4 cusps; documented mult.≥2 reliability case) | clean_success | 15.7 | S4 | 4 | 2.31e-09 | 5.89e-08 |
| nodal_cubic | new — node singular structure, contrasted vs. pure-node-curve | clean_success | 17.3 | C1/S1/B2 | 4 | 2.29e-10 | 5.49e-07 |
| cuspidal_cubic | new — single generic cusp (simpler than astroid's 4) | clean_success | 16.9 | S1/B2 | 2 | 1.32e-13 | 3.43e-07 |
| squircle_quartic | new — smooth-critical-point misclassification | success_with_caveats | 15.4 | S2/A4 | 2 | 3.81e-10 | 9.60e-07 |
| lemniscate_bernoulli | new — self-crossing node + turning point | clean_success | 15.6 | C2/S1 | 4 | 2.84e-08 | 7.67e-07 |
| folium_descartes | new — node + unbounded branch; surfaced nondeterminism | success_with_caveats | 15.9 | C1/S1/B2/A1 | 5 | 1.55e-11 | 7.51e-07 |
| cissoid_diocles | new — cusp + unbounded branch | clean_success | 17.7 | S1/B2 | 2 | 2.32e-11 | 9.16e-10 |
| tangent_parabolas_reducible | new — genuinely reducible (2 factors) + tacnode | clean_success | 17.1 | S1/B4 | 4 | 3.09e-11 | 7.51e-09 |
| three_concurrent_lines_reducible | new — genuinely reducible (3 factors) + triple pt. + axis-aligned component | success_with_caveats | 18.4 | B6/A1 | 4 | 0.0 | 0.0 |
| elliptic_two_components | new — disconnected-but-irreducible | clean_success | 16.9 | C3/B2 | 4 | 7.70e-10 | 5.20e-07 |
| empty_curve | new — empty real locus (2D) | empty_real_locus | 13.6 | (none) | 0 | n/a (0 pts) | n/a |

### Surfaces (`decompose_3d_surface`)

| Fixture | Category / gap targeted | Outcome | Wall (s) | Vertices | E/F | Mesh V/Tri | Median resid. | Max resid. |
|---|---|---|---|---|---|---|---|---|
| sphere | baseline | clean_success | 23.0 | C2 | 2/2 | 152/300 | 2.14e-08 | 5.61e-08 |
| ellipsoid | baseline (documented "ellipsoid discovery" fixture) | success_with_caveats | 23.6 | C2 | 2/2 | 160/316 | 2.57e-08 | 1.00e-05 |
| taubin_heart | baseline (Phase 5/8/9 fixture; singular notch) [†] | clean_success | 34.2 | C10/A4 | 14/14 | 1174/2236 | 2.32e-08 | 2.85e-05 |
| torus | baseline (`projection=:random` fixture, seed 42) [†] | clean_success | 40.1 | C8 | 8/8 | 680/1280 | 1.44e-06 | 7.68e-04 |
| whitney_umbrella | baseline (documented crash prediction) | exception | 15.9 | — | — | — | — | — |
| horn_torus | new — degenerate topology (self-tangent pinch, `z=0`) | success_with_caveats | 16.7 | all 0 | 0/0 | 0/0 | n/a | n/a |
| cone | baseline (apex; explicitly-"untested" prediction) | success_with_caveats | 16.8 | all 0 | 0/0 | 0/0 | n/a | n/a |
| hyperboloid_one_sheet | new — unbounded, connected topology [‡] | clean_success | 23.0 | C2 | 2/2 | 170/320 | 1.45e-07 | 8.22e-07 |
| hyperboloid_two_sheets_reducible | new — disconnected-but-irreducible (mislabeled "reducible" — see §2.1) | clean_success | 23.7 | C4 | 4/4 | 322/620 | 1.71e-07 | 1.09e-06 |
| elliptic_paraboloid | new — unbounded topology | clean_success | 23.8 | C2 | 2/2 | 161/310 | 4.83e-08 | 2.72e-07 |
| quartic_superellipsoid | new — surface analog of `squircle_quartic` misclassification | success_with_caveats | 24.1 | S2/A4 | 2/2 | 204/320 | 5.26e-08 | 1.38e-06 |
| empty_surface | new — empty real locus (3D) | empty_real_locus | 15.6 | all 0 | 0/0 | 0/0 | n/a | n/a |

`whitney_umbrella`'s JSON has every numeric field (`vertex_counts`,
`edge_count`, `face_count`, `residuals`) set to `null`/`None` — confirmed
directly, not assumed — because `decompose_3d_surface` threw before
producing any output; only `error` (the full stacktrace) is populated.

**[†] taubin_heart / torus** — both `outcome=clean_success` per their JSON
(no exception, no logged discrepancy at time of the original survey), but
their renders visibly show mid-surface gaps that were not flagged as
findings until a 2026-08-06 follow-up. CONFIRMED as the known
no-`incidence=true` naked-edge mechanism (`docs/DESIGN_NOTES.md`'s
"Watertightness measurements"): bare→incidence=true naked-edge counts,
torus 80→8 (10.0x), taubin_heart 132→29 (4.6x); `incidence=true`
substantially reduces but does NOT eliminate the visible gaps for either
fixture, matching that entry's own documented caveat exactly. Full detail
in `errors_log.md`'s 2026-08-06 entry; the `clean_success` tag itself is
not being revised (both fixtures behave exactly as this already-known,
already-documented limitation predicts — this is not a new bug), but the
visual gap should not have gone unremarked in this report's original
version.

**[‡] hyperboloid_one_sheet** — `outcome=clean_success`, but its render
shows a visibly jagged/faceted rim at the bbox-clipped boundary
(`z=±4`), also unflagged in the original pass. CONFIRMED as a coarse-
sampling artifact specific to bbox-clipped unbounded surfaces (a
DIFFERENT mechanism from taubin_heart/torus above, not the same
naked-edge issue) — the rim smooths out almost completely at default
`edge_sample_density=50` (this survey used coarse `6` throughout). Full
detail in `errors_log.md`'s 2026-08-06 entry.

---

## 2. Patterns by category

### 2.1 Reducible vs. disconnected-but-irreducible — NOT the same class, and one fixture is mislabeled

**A naming problem found first, because it changes the answer.**
`hyperboloid_two_sheets_reducible`'s equation is `x^2+y^2-z^2+1`. Factoring
all four "reducible"/"disconnected" fixtures directly (independent of
anything Phase 1 claimed):

```
tangent_parabolas_reducible:          -x**4+y**2  -> -(x**2-y)*(x**2+y)   [factors — genuinely reducible]
three_concurrent_lines_reducible:     x**3-x*y**2 -> x*(x-y)*(x+y)        [factors — genuinely reducible]
hyperboloid_two_sheets_reducible:     x**2+y**2-z**2+1 -> unchanged        [does NOT factor]
elliptic_two_components:              -x**3+x+y**2 -> unchanged           [does NOT factor]
```

`x^2+y^2-z^2+1` is a rank-4 quadratic form (homogenizing to
`x²+y²-z²+w²`, matrix `diag(1,1,-1,1)`, full rank) — a product of two
linear polynomials always has associated-form rank ≤ 2, so a rank-4 form
can never factor that way, over any field. This isn't a marginal
distinction: **`hyperboloid_two_sheets_reducible` is a single irreducible
quadric whose real locus happens to have two disconnected sheets — the
same phenomenon as `elliptic_two_components`, not the same phenomenon as
the two fixtures that actually do factor.** The fixture's own name
(`..._reducible`) is a misnomer against its actual algebra. This is
independent of anything Phase 1's log claims (neither `errors_log.md` nor
`SURVEYOR_DONE.md` addresses this), so it is reported here as a fresh
cross-check, not attributed to either.

With that correction, the real comparison is:

| Class | Members | Outcome |
|---|---|---|
| **Genuinely reducible** (multiple polynomial factors) | `tangent_parabolas_reducible`, `three_concurrent_lines_reducible` | 1 clean, 1 caveats |
| **Genuinely disconnected-but-irreducible** (one factor, split real locus) | `elliptic_two_components`, `hyperboloid_two_sheets_reducible` | 2/2 clean |

The disconnected-but-irreducible class is uniformly clean in this survey.
The reducible class is mixed — but the failure in
`three_concurrent_lines_reducible` is traceable to two specific,
independent mechanisms (an axis-aligned vertical-line component invisible
to x-parametrized slicing, and a triple point misclassified `Artificial`)
that have nothing to do with "reducibility" as an abstract property:
`tangent_parabolas_reducible` is equally reducible (2 factors) and
resolved perfectly (tacnode correctly `Singular`, all 4 boundary exits
wired). Neither pipeline (`decompose_1d_curve` for both curve pairs,
`decompose_3d_surface` for the hyperboloid) has any code path that treats
"reducible" and "disconnected-irreducible" differently — both classes go
through the same generic real-locus slicing. **Conclusion: the survey
data does not show reducibility itself driving different behavior from
disconnectedness; it shows the SAME generic pipeline hitting unrelated,
fixture-specific edge cases (axis alignment, singularity order) that
happen to land more often in the small "reducible" sample here.**

### 2.2 Cone vs. the astroid/Whitney-umbrella multiplicity issue — same backlog entry, three distinct manifestations

`docs/DESIGN_NOTES.md`'s "HC.jl polyhedral-solve reliability on
multiplicity≥2 / reducible critical-point systems" backlog entry
consolidates four prior cases under one root cause (HC.jl's default
polyhedral `solve` struggling when the critical-point system has a
multiplicity-≥2 solution): node curve, astroid cusps, Whitney umbrella,
and — explicitly flagged **"untested"** — the cone apex. This survey's
`cone.json` (`decompose_3d_surface` on `x^2+y^2-z^2`, bare/default call)
returns a completely empty decomposition (0/0/0/0), no exception, in
16.8s — **this is the first live confirmation of that specific
"untested" sub-entry**, matching its prediction exactly (empty
`compute_critical_z_slices`, symmetric `bbox_z` landing the naive
midpoint on the undetected apex).

But "same backlog entry" does not mean "same failure signature." Within
that one root cause, this survey's own data shows three genuinely
different observable behaviors:

- **Astroid** (this survey, `clean_success`): all 4 cusps found as
  `Singular`, 0 `Artificial` — one of the "good" runs of a documented
  ~3/5-good, ~2/5-miss-1-or-2-cusps nondeterministic pattern. Failure
  mode when it occurs: **partial, nondeterministic miss** (some cusps
  silently become `Artificial` fallback vertices, others fine, varies
  run to run).
- **Whitney umbrella** (this survey, `exception`): `compute_critical_z_slices`
  returns empty AND the resulting `z=0` midpoint itself throws inside HC.jl's
  `polyhedral` start-system construction. Failure mode: **deterministic
  hard crash**, `OverflowError`, loud.
- **Cone** (this survey, `success_with_caveats`, functionally empty):
  `compute_critical_z_slices` also returns empty, but the `z=0` midpoint
  does *not* throw — `decompose_3d_surface` just silently returns nothing.
  Failure mode: **deterministic silent total emptiness**, no error at all.

So the cone's confirmed behavior is *not* a new failure signature and
*not* the astroid's partial/nondeterministic pattern — it is exactly the
already-predicted "no crash this time, just nothing, which is arguably
worse to notice" case docs item 4 describes, now backed by live data for
the first time.

**A separate mechanism that must not be conflated with the above:**
`squircle_quartic` and its surface analog `quartic_superellipsoid` show a
*different* documented mechanism — `docs/DESIGN_NOTES.md`'s "Stage 4c —
`_deflation_applicable`" entry (a point can be `Faug`-singular while
`F_original`'s own gradient is nonzero, using `f=x-y^3` as its abstract
example). Here HC.jl's solve *succeeds* at finding the point (unlike the
multiplicity≥2 backlog cases, where solve fails to find it at all); the
bug is downstream classification (`Singular` instead of `Critical`) and,
worse and undocumented anywhere, edge construction bypassing the
correctly-found points entirely in favor of off-curve `Artificial`
fallbacks. This is confirmed as visually real, not just a bookkeeping
label, by direct inspection of both renders (see §4).

### 2.3 Empty-locus fixtures — graceful decomposition in both dimensions, but plotting diverges 2D vs. 3D

`empty_curve` (`x^2+y^2+1`) and `empty_surface` (`x^2+y^2+z^2+1`) both
return `outcome=empty_real_locus`, all counts 0, `error=null` — confirmed
directly from both JSONs. Decomposition itself is graceful and identical
in kind across dimensions.

**Plotting is not identical.** `empty_curve.png` is a genuine,
successfully-produced blank `[-4,4]×[-4,4]` axes frame — no crash.
`empty_surface.json`'s own `notes` field records a `PLOT FAILED`
`ArgumentError: reducing over an empty collection is not allowed` from
`plot_surface_decomposition` → `_near_constant_colorrange`
(`src/Visuals.jl:210`/`:282`) — the same crash also hit by `horn_torus`
and `cone` (3 occurrences total, all from an empty *mesh*, since the 2D
curve-plotting path apparently has no equivalent `extrema`-over-empty-
collection call). So: **decomposition behaves identically 2D vs. 3D
(graceful); plotting does not (2D tolerates empty, 3D crashes).**

One correction to how Phase 1 framed this: both `errors_log.md`'s absence
of an empty-locus entry and `SURVEYOR_DONE.md`'s phrase "confirms the
predicted graceful-empty behavior" imply this was an existing documented
prediction. Grepping `docs/DESIGN_NOTES.md` in full for "empty": the only
4 hits are all about `compute_critical_z_slices` returning an empty
*critical-value set* (the cone/Whitney/torus-fold discussions), never
about a fixture's real locus being empty outright. **There is no
`docs/DESIGN_NOTES.md` entry documenting expected behavior for a wholly
empty real variety** — so "graceful handling of an empty real locus" is,
as far as this repo's own design-notes file goes, a new finding (a
positive one), not a confirmation of anything written there.

---

## 3. Timing outliers and residual anomalies

**Wall time.** Curve fixtures cluster tightly: 13.6s (`empty_curve`,
fastest) to 18.4s (`three_concurrent_lines_reducible`), no real outliers.
Surfaces split into three bands:
- **Successful, non-trivial surfaces** (`sphere`, `ellipsoid`,
  `hyperboloid_one_sheet`, `hyperboloid_two_sheets_reducible`,
  `elliptic_paraboloid`): tight band, 22.96–24.11s.
- **Outliers, both slower**: `taubin_heart` 34.2s (~1.45x that band) and
  `torus` 40.1s (~1.71x, the single slowest fixture in the whole survey).
  Both are explicable by real structural cost (taubin: 14 edges/faces vs.
  the band's 2–4; torus: `projection=:random` plus 8 edges/faces), not
  flagged as anomalous, but worth naming since they're the two furthest
  outliers by wall time in the entire 24-fixture set.
- **The four `success_with_caveats`/`exception`/`empty_real_locus`
  surfaces that terminate early** (`whitney_umbrella` 15.9s, `horn_torus`
  16.7s, `cone` 16.8s, `empty_surface` 15.6s): all *faster* than the
  healthy-surface band, ~65–73% of it. This is an artifact worth flagging
  explicitly so it isn't misread as "fine, ran fast" — these fixtures are
  fast *because* they short-circuit into a crash or an empty return
  before ever reaching the expensive mesh/face-tracking stage, not
  because anything about them is cheap to correctly decompose.

**Residual anomalies relative to the sphere/ellipsoid baseline (median
~2e-8, this survey's own sphere: median 2.136e-08, max 5.610e-08):**

- **`ellipsoid`**: median/p90 sit in baseline range, but p99/max jump to
  1.0003e-05 — ~170x above its own p90, ~180x above sphere's own max.
  `n_points=160`, so this is ~1–2 outlier mesh points, not a broad
  degradation.
- **`torus`**: median 1.44e-06 (~67x sphere's median), max 7.68e-04
  (~13,700x sphere's median, ~347x sphere's own max). The fixture's own
  `notes` field attributes this to the deliberately coarse
  `edge_sample_density=6` used by this survey vs. a denser validation
  script elsewhere, and does not treat it as a discrepancy. Plausible,
  but no docs entry states an expected torus residual at this specific
  coarse density.
- **`taubin_heart`**: median 2.32e-08 (baseline), but max 2.85e-05 and
  p99 2.27e-06. `docs/DESIGN_NOTES.md`'s "Validation of the refined
  gates" entry (Phase 8, `_robust_slice_at_z`) states, for the *same*
  fixed-axis (non-rotated) Taubin heart, "full fixed-axis decompose max
  `|f| = 2.4e-6`" as the current post-fix baseline (`DESIGN_NOTES.md:734`,
  independently re-confirmed for this report). This survey's max
  (2.85e-05) is **~12x higher** than that documented figure. The
  fixture's own JSON `notes` field attributes the elevation to a
  *different* docs entry — "Known limitation: generic projections over
  singular curves" — but that entry is specifically about the **rotated**
  Taubin heart under `projection=:random`; this survey's `taubin_heart`
  config carries no `projection` key, i.e. it's the fixed-axis run. That
  makes the fixture's own self-diagnosis a mismatch against the docs text
  it cites — the entry it should be compared against is the fixed-axis
  "2.4e-6" figure, not the rotated-projection one, and against that
  figure there's a real, unexplained 12x gap. Plausibly the same
  coarse-density effect seen in `torus` (methodology, not confirmed
  identical), but this specific comparison is not addressed anywhere in
  Phase 1's own log.
  **2026-08-06 update**: both `torus`'s and `taubin_heart`'s elevated max
  residuals were investigated alongside their visible render gaps (see
  §1's [†] footnote and `errors_log.md`'s 2026-08-06 entry) — both are
  CONFIRMED to be run without `incidence=true`, and both show substantial
  naked-edge counts in that state (torus 80, taubin_heart 132, vs. 8/29
  once `incidence=true` is added). Elevated max/p99 residuals sitting at
  unstitched mesh-crack boundaries is a plausible connection to the same
  underlying no-incidence limitation, but this was NOT directly measured
  (i.e., residuals were not recomputed with `incidence=true` in this
  follow-up) — flagged as a plausible link, not a confirmed one.
- **`nodal_cubic`** and **`lemniscate_bernoulli`**: both self-intersecting
  curves, max residuals 5.49e-07 and 7.67e-07 respectively (~19x/~26x
  circle's own max of 2.91e-08), but medians (2.29e-10, 2.84e-08) stay in
  baseline range — consistent with locally tighter Newton conditioning
  near a node, a plausible and minor effect, not flagged as a discrepancy.
- **`hyperboloid_one_sheet`/`hyperboloid_two_sheets_reducible`**: medians
  1.45e-07/1.71e-07, ~7–8x sphere's median but still same order of
  magnitude — mild, plausibly a bbox-clipped-unbounded-surface effect, not
  previously quantified anywhere in `docs/DESIGN_NOTES.md` for either
  fixture (neither is mentioned there at all).

---

## 4. Confirms-documented vs. new-finding ledger

Every finding below is tagged explicitly, per-item, against
`docs/DESIGN_NOTES.md`'s actual text.

1. **`squircle_quartic` — smooth points at `(±1,0)` misclassified
   `Singular`.** *Mechanism* — **Confirms `docs/DESIGN_NOTES.md`'s "Stage
   4c — `_deflation_applicable`" entry**, applied concretely to a new
   fixture: the entry's own abstract illustration (`f=x-y^3`, `Faug`
   singular while `F_original` is regular) is exactly what's measured here
   (`jacobian_rank=1`, `singular_values=[4,0]`, `F_original` gradient
   `(4,0)` full rank 1). *Edge-graph wiring to 4 off-curve `Artificial`
   fallbacks (residual up to 1.0 on the surface analog) instead of the
   correctly-found points* — **new finding, not covered by the Stage 4c
   entry or anywhere else in `docs/DESIGN_NOTES.md`.** Visually confirmed
   by direct inspection of `renders/squircle_quartic.png`: only two tiny
   disconnected arcs are drawn near the corner regions; the two `Singular`
   markers at `(±1,0)` are isolated points with no curve segment reaching
   them at all — most of the true closed quartic boundary is simply
   absent from the render, not merely mislabeled.
2. **`quartic_superellipsoid` — surface analog of #1.** Same split:
   misclassification mechanism **confirms Stage 4c**; `Artificial`
   edge-endpoints at residual ~1.0 and a **visible wedge-shaped
   notch/gouge in the mesh** (confirmed directly by rendering
   `quartic_superellipsoid.png` and comparing against `sphere.png`: the
   sphere is a fully closed smooth ball, the superellipsoid clearly has a
   flat cut wedge removed near one pole) — **new finding.**
3. **`folium_descartes` — 1 of 3 identical runs silently drops a curve
   segment.** The *general* mechanism (cross-process HC.jl solver jitter)
   **confirms** patterns already documented elsewhere in
   `docs/DESIGN_NOTES.md` (the torus naked-edge-count spread under
   "Watertightness measurements," and the full-suite 536-vs-537 flake
   under "Full-suite test timing and the 536-vs-537 assertion-count
   variance"). But this *specific consequence* — `decompose_1d_curve`
   silently orphaning a genuine `Boundary` vertex and dropping a whole
   curve segment with zero error signal — is **a new finding**; no prior
   entry documents jitter manifesting this way on this pipeline.
4. **`three_concurrent_lines_reducible` — vertical line component
   (`x=0`) entirely absent from the edge graph.** Grepped
   `docs/DESIGN_NOTES.md` for "vertical line" and related terms: zero
   hits. **New finding** — a structural limitation (x-parametrized
   slicing cannot represent a component existing at exactly one x-value)
   with no precedent in the docs.
5. **`three_concurrent_lines_reducible` — triple point misclassified
   `Artificial`, not `Singular`.** The underlying mechanism (a genuine
   multiplicity-3 zero defeating HC.jl's polyhedral solve reliability) is
   the concrete, first-ever-measured instance of a case the "HC.jl
   polyhedral-solve reliability" backlog entry's own introductory
   sentence explicitly names in the abstract ("structurally guaranteed at
   cusps, **reducible crossing points**, and conical apexes") but had
   never previously backed with data. **Confirms the abstract
   `docs/DESIGN_NOTES.md` "reducible crossing points" case** with its
   first concrete measurement — not a wholly new mechanism, but the first
   evidence for an already-named-but-previously-untested sub-case.
   Visually confirmed via `renders/three_concurrent_lines_reducible.png`:
   only the X-shaped diagonal pair (`x=y`, `x=-y`) is drawn; the vertical
   line and its two `Boundary` markers at `(0,±4)` sit disconnected with
   no line reaching them; the origin is plotted as a gray `Artificial` X,
   not a `Singular` marker.
6. **`ellipsoid` — max/p99 residual ~170x above p90.** The *mechanism*
   (fixed patch-line transversality loss nearer a pole, asymmetric
   surface) is **the same story as `docs/DESIGN_NOTES.md`'s "Adaptive
   re-anchoring"/"The ellipsoid discovery" entry** (Face tracking, Phase
   5), which used the ellipsoid to originally discover the problem. But
   that entry's own validation *measurements* (the two-gate mechanism's
   before/after numbers) are all taken on the Taubin heart, not the
   ellipsoid itself — there is no prior documented residual figure for the
   ellipsoid fixture post-fix. **This specific 170x/1.0e-05 quantification
   is a new data point**, plausibly consistent with, but not previously
   recorded by, the existing narrative.
7. **`whitney_umbrella` — `OverflowError: Cannot compute a start
   system.`** **Confirms `docs/DESIGN_NOTES.md`'s "Uncaught exception
   when `compute_critical_z_slices` finds zero critical z-values and the
   naive bbox midpoint is itself degenerate" entry** exactly — stacktrace
   terminal frame independently verified in this survey's own JSON to
   match `Solver.jl:873` → `intersect_bounding_object`, as the docs entry
   describes.
8. **`horn_torus` — completely empty decomposition, no exception.**
   The fixture itself is untested territory (`docs/DESIGN_NOTES.md` never
   mentions "horn torus," confirmed via grep), so at the fixture level
   this is a **new finding**. But the mechanism it demonstrates is **not
   new** — it is the identical already-documented pattern as the cone
   (item 4 of "HC.jl polyhedral-solve reliability"): empty
   `compute_critical_z_slices`, symmetric `bbox_z`, naive midpoint landing
   on a genuine degenerate point (here a self-tangency/pinch at the
   origin rather than a cone apex). Reported as: new fixture, confirming
   mechanism — do not read as either a pure confirmation or a pure novel
   bug.
9. **`cone` — completely empty decomposition, no exception.**
   **Confirms `docs/DESIGN_NOTES.md`'s backlog item 4** — explicitly
   flagged there as "untested" prior to this survey. This is the first
   live data for that specific sub-entry.
10. **`plot_surface_decomposition` throws `ArgumentError` on any
    completely empty mesh** (hit 3x: `horn_torus`, `cone`,
    `empty_surface`, same stack each time, from `_near_constant_colorrange`
    at `src/Visuals.jl:210`). `docs/DESIGN_NOTES.md` does have an entry
    for `_near_constant_colorrange` ("Phase 6 '02b' discovery"), but it
    documents a completely different failure — Float32 round-off on a
    *non-empty*, near-constant colorbar range causing full-spectrum
    speckle, not an empty-collection `extrema` crash. **New finding**,
    same function, unrelated bug.
11. **`empty_curve`/`empty_surface` — graceful empty-locus handling.**
    Grepped `docs/DESIGN_NOTES.md` for "empty" in full: all 4 hits concern
    an empty *critical-value set* mid-pipeline (cone/Whitney/torus-fold
    discussions), never a wholly empty real variety. **No
    `docs/DESIGN_NOTES.md` entry documents this scenario at all** — so,
    contrary to how `SURVEYOR_DONE.md`'s phrasing ("confirms the predicted
    graceful-empty behavior") reads, this is **a new finding** as far as
    `docs/DESIGN_NOTES.md` is concerned, not a confirmation of anything
    written there.
12. **`taubin_heart` — max residual 2.85e-05 vs. `docs/DESIGN_NOTES.md`'s
    own quoted fixed-axis full-decompose figure of `2.4e-6`.** See §3 for
    the full number. **New finding** — a ~12x gap against a specific,
    citable documented number for the same (fixed-axis, non-rotated)
    configuration, not resolved by the rotated-Taubin entry the fixture's
    own `notes` field cites instead.
13. **`nodal_cubic` node correctly resolved as `Singular`.** Genuinely
    **contrasts with, rather than confirms**, `docs/DESIGN_NOTES.md`'s
    documented pure-node-curve case (`y²-x²`, "returns zero solutions"):
    this is a *different* fixture (`y²-x²(x+1)`, a nodal cubic, not the
    bare two-crossing-lines curve) that happens to resolve cleanly where
    the documented pure-node case fails. Correctly described this way in
    the fixture's own `notes` field — verified accurate, not a
    discrepancy.
14. **`hyperboloid_one_sheet`, `hyperboloid_two_sheets_reducible`,
    `elliptic_paraboloid`** (all clean, unbounded-topology fixtures) and
    **`tangent_parabolas_reducible`, `cissoid_diocles`,
    `lemniscate_bernoulli`, `elliptic_two_components`, `cuspidal_cubic`**
    (all clean, singular/multi-component fixtures): none of these
    fixtures or their specific numbers appear anywhere in
    `docs/DESIGN_NOTES.md` (confirmed via grep — "hyperboloid," "folium,"
    "cissoid," "lemniscate" all return zero hits). **All new fixtures,
    all with a clean/positive result** — worth stating plainly as new
    positive coverage, not filed as "findings" in the bug sense since
    nothing anomalous occurred.

---

## 5. Full reproduction of `errors_log.md`

`errors_log.md` contains exactly 7 top-level entries (`##` headers, one
per discrepancy or confirmation logged in real time; 2 of the 7 headers —
`horn_torus`'s two entries — each contribute a further sub-finding, so 9
distinct findings total). Every one is accounted for above; this section
maps each verbatim so nothing is lost in the higher-level synthesis:

1. **`squircle_quartic` — genuine smooth critical points misclassified
   Singular; edges built on off-curve Artificial placeholders instead.**
   Full detail: 2 correctly-located points at `(1,0)`/`(-1,0)`,
   `jacobian_rank=1`, `singular_values=[4.0,0.0]`; both edges wired to 4
   Artificial `:endpoint_fallback` vertices at `(±1,±1)` (off-curve,
   `|f(1,1)|=1`); render shows only 2 disconnected arc fragments, most of
   the curve missing. → covered in §2.2, §4 item 1.
2. **`folium_descartes` — run-to-run nondeterminism: identical
   fixture/config, 1 of 3 runs silently drops a curve segment.** Good run:
   5 vertices/5 edges, all referenced. Bad run: vertex id=4 (a genuine,
   on-curve `Boundary` vertex) referenced by zero edges; the unbounded
   branch's continuation segment is silently dropped, no exception. →
   covered in §4 item 3.
3. **`three_concurrent_lines_reducible` — vertical line component
   entirely absent from the edge graph; triple point misclassified
   Artificial not Singular.** Only the `x=y`/`x=-y` lines appear (4
   edges); the `x=0` line's 2 correctly-located boundary vertices are
   never referenced by any edge; the origin (triple point) is found and
   used to stitch the other two lines but typed `Artificial`, not
   `Singular`, 0 `Singular` vertices total. → covered in §2.1, §4 items
   4–5.
4. **`ellipsoid` — max residual ~170x above p90, small number of outlier
   mesh points.** Median `2.573e-08`/p90 `5.975e-08` in baseline range;
   p99/max both `~1.0003e-05`; `n_points=160`, roughly 1–2 outlier points.
   → covered in §3, §4 item 6.
5. **`whitney_umbrella` — CONFIRMATION of documented `OverflowError`
   prediction (not a new finding).** Reproduced exactly: `OverflowError:
   Cannot compute a start system.`, `intersect_bounding_object ->
   decompose_1d_curve -> slice_at_z -> _robust_slice_at_z ->
   decompose_3d_surface`, terminal frame `Solver.jl:873`. → covered in §4
   item 7.
6. **`horn_torus` — silent completely-empty decomposition, extending the
   documented "cone" prediction to a NEW fixture.** `(x^2+y^2+z^2)^2 -
   4(x^2+y^2)`, self-tangent pinch at the origin; `decompose_3d_surface`
   returns all-zero `(vertices=[], edges=[], faces=[], mesh)`, no
   exception, `wall_time=16.66s`. → covered in §2.2, §4 item 8.
7. **`horn_torus` — separate finding: `plot_surface_decomposition` throws
   on an empty mesh instead of handling it gracefully.**
   `ArgumentError: reducing over an empty collection is not allowed` from
   `_near_constant_colorrange` (`src/Visuals.jl:210`) via
   `plot_surface_decomposition` (`:282`). → covered in §4 item 10 (noted
   there as a 3x-occurring, single-root-cause issue across `horn_torus`,
   `cone`, `empty_surface`).
8. **`cone` — CONFIRMATION of documented "untested" silent-empty
   prediction (fixture 22, high priority).** Bare/default call returns
   all-zero decomposition, no exception, `wall_time=16.77s`; also hits the
   same `plot_surface_decomposition` crash as `horn_torus`. → covered in
   §2.2, §4 item 9.
9. **`quartic_superellipsoid` — surface-level recurrence of the
   squircle_quartic mechanism, worse: Artificial edge endpoints residual
   ~1.0 (nowhere near the surface).** Same misclassification mechanism as
   entry 1, plus a render-confirmed wedge-shaped notch/gouge near one
   pole (correction made in the log itself after visual inspection, kept
   as recorded). → covered in §4 item 2.

No entry in `errors_log.md` was dropped, merged away, or omitted from
this report.

---

## 6. What this report does not do

Per scope: no `src/` fixes are proposed anywhere above. Every item is
either a direct re-statement of measured data (JSON fields, render pixels,
`docs/DESIGN_NOTES.md` grep results) or an explicit confirms/new
determination against that data — diagnostic only.

---

## 7. Note on this file's provenance

The Phase 2 evaluator subagent that produced the content above had its
`Write` tool call blocked by a system-level restriction (subagents return
findings as text to the orchestrating session, not files) — so this file
did not exist on disk until the orchestrating Developer session wrote it
from the evaluator's returned text. Two `docs/DESIGN_NOTES.md` citations
central to §3/§4 (the taubin_heart "2.4e-6" fixed-axis figure, and the
"empty" grep claim) were independently re-verified by the Developer
session against the live file before this report was finalized, and both
held up exactly as stated.
