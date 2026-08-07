# Design notes

Investigation history, rejected alternatives, and measurement records for
this codebase's non-obvious design decisions — separated out (Audit 2,
2026-07) from the docstrings/comments where this material used to live
inline, so that source-level documentation stays focused on *current,
load-bearing* usage and maintenance rules, while the *story* of how each
rule was arrived at stays discoverable here instead of disappearing.

Each entry is dated and cites the source file/function it backs. Where a
source site still needs a rule stated for correct use, that rule stays in
the docstring; only the investigation narrative moved.

---

## Backlog: production robustness gaps

Real, live-confirmed gaps found during other investigations, logged here
rather than fixed inline because each needs its own separately-scoped
design decision. No existing consolidated backlog list was found in this
repo to add these next to (grepped `docs/`, `README.md`; nothing found),
so this section is new as of 2026-07. The Griffis-Duffy singular-curve
blueprint (a different, larger item, full detail kept externally — see
its own entry below) lives in this same section, not a separate one.

### TagBot-created tags never trigger a "stable" docs build — needs a PAT, not implemented here

Found 2026-07 during the v0.2.0 registry pre-flight. `.github/workflows/
Documentation.yml` is correctly configured for Documenter.jl's dev/stable
split (`on: push: tags: [v*]`, plus the default `deploydocs` versions
scheme), but checking the actual `gh-pages` branch content (not just the
workflow config) showed only a `dev/` folder exists — no `stable/`,
`v0.1.0/`, or `v0.1/` — and `versions.js` has `DOCUMENTER_STABLE = "dev"`,
i.e. stable is just aliased to dev because nothing else was ever built.
Root cause, confirmed by checking every historical workflow run: the only
existing tag, `v0.1.0`, was created by `.github/workflows/TagBot.yml`
using `secrets.GITHUB_TOKEN`. Pushes/tags made with the default
`GITHUB_TOKEN` do not retrigger other Actions workflows (GitHub's built-in
anti-recursion restriction) — every run in the list is a `push` event
with `headBranch: main`, none show a tag-ref trigger, confirming the
tag-triggered docs build has simply never fired. Low-risk partial fix
applied same day: added `workflow_dispatch:` to `Documentation.yml` so a
stable build can be kicked off manually (`gh workflow run
Documentation.yml`) whenever a real tag exists but didn't auto-trigger.
The full automatic fix — giving `TagBot.yml` a personal access token
(instead of `GITHUB_TOKEN`) so its own tag pushes *do* retrigger other
workflows, per TagBot's own documented recommendation for this exact
gap — was deliberately **not** implemented: creating a PAT is a
credential-issuing action that has to happen on Juan's own GitHub account
directly, not something to do on his behalf. Left as a follow-up decision
whenever the manual-dispatch workaround becomes enough of a recurring
annoyance to be worth it.

### Uncaught exception when `compute_critical_z_slices` finds zero critical z-values and the naive bbox midpoint is itself degenerate

Found 2026-07 while trying the Whitney umbrella (`x^2-y^2*z=0`) as a
candidate `decompose_3d_surface` validation fixture (see the Whitney
umbrella item below; this gap is independent of Whitney specifically).
`compute_critical_z_slices` returned `Float64[]` for this surface (no
critical z found at all), so `_slab_bounds` never split away from the
naive bbox midpoint — and that naive midpoint, `z=0`, is exactly where
`slice_at_z(F, 0.0, cfg)` itself throws `OverflowError: Cannot compute a
start system` (confirmed directly, not inferred). `_robust_slice_at_z`'s
existing retry mechanism (`docs/DESIGN_NOTES.md`, "Robust z-mid
selection") only catches the *softer* signal of a slice coming back with
bad vertices (`:endpoint_fallback` + `Singular` co-occurrence) — it has
no `try`/`catch` around the tracking calls themselves, so an outright
exception from a degenerate naive midpoint propagates uncaught and kills
the whole `decompose_3d_surface` call with no retry attempted. This is a
real production robustness gap independent of any specific fixture:
**any** surface whose only real critical z lands outside `bbox_z`, or
which has none at all, is one unlucky bbox choice away from this same
crash if its bbox midpoint happens to be numerically degenerate.

### HC.jl polyhedral-solve reliability on multiplicity≥2 / reducible critical-point systems

One recurring pattern, found independently four separate times during
2026-07 investigations — consolidated here into a single entry rather
than left scattered across separate bullets, so the connection isn't
missed. In every case, the critical-point-finding augmented system has a
solution of multiplicity ≥ 2 at exactly the point being searched for
(structurally guaranteed at cusps, reducible crossing points, and conical
apexes — that multiplicity *is* what makes them singular), and HC.jl's
default polyhedral `solve` unreliably finds it:

1. **Node curve** (`y^2-x^2=0`, two crossing lines): `compute_critical_points`
   on the y-derivative-augmented system `[y^2-x^2, 2y]` returns **zero**
   solutions — both tracked paths report `return_code=excess_solution`,
   neither converging to the genuine double root at the origin. See
   `paper_artifacts/VISUAL_ASSETS.md` Finding 1.
2. **Astroid's 4 cusps** (corrected equation, `(x^2+y^2-1)^3+27x^2y^2=0`):
   each cusp is a multiplicity-2 solution of its own augmented system
   (`solve` reports 12 raw paths converging to 4 physical points —
   several redundant paths per cusp). Across 5 independent runs, 2 came
   back with 1 or 2 cusps missing entirely (replaced by `Artificial`
   fallback vertices), not just coordinate jitter. See
   `paper_artifacts/VISUAL_ASSETS.md` Finding 5.
3. **Whitney umbrella's apex/singular z-axis** (`x^2-y^2*z=0`, tried as a
   `decompose_3d_surface` fixture): `compute_critical_z_slices` returned
   `Float64[]` — no critical z found at all — which is what let the naive
   bbox-midpoint slab land exactly on the degenerate `z=0` plane and crash
   (see the uncaught-`OverflowError` entry above for that downstream
   consequence; the empty critical-z result itself is this pattern).
4. **Cone's apex** (`x^2+y^2-z^2=0`, tried as a Whitney alternative): the
   apex is a genuine isolated critical point mathematically (`x=y=0`
   forced, then `z^2=0` — a double root in `z`), but
   `compute_critical_z_slices` returned `Float64[]` for it too,
   live-confirmed. Because the cone's `bbox_z` is symmetric, the
   resulting single whole-bbox slab's naive midpoint lands exactly on
   that undetected apex, and `decompose_3d_surface` silently returns a
   **completely empty** decomposition (0 vertices, 0 faces) — no crash
   this time, just nothing, which is arguably worse to notice.
   Independently reproduced 2026-08-06 by a live capability-survey run of
   the bare, default `decompose_3d_surface` call itself (not just
   `compute_critical_z_slices` in isolation): 0 vertices, 0 edges, 0
   faces, 0 mesh points, no exception — saved evidence at
   `dev/scratch/capability_survey/data/cone.json`.
   **2026-08-07 follow-up — neither tested workaround is safe to use**:
   - Asymmetric `bbox_z=(-4.0, 4.3)` (coarse survey density) produces a
     *non-empty* decomposition (2 vertices, 2 edges/faces, 98 mesh
     points) that looks superficially plausible but is geometrically
     corrupted: median residual is fine (`1.2e-7`), but `p99=12.2`,
     `max=16.0` — a visible seam/discontinuity around `z≈2.2-2.5` plus a
     spike artifact near the rim, confirmed by direct render inspection,
     not just the residual numbers. **This is more dangerous than the
     original empty-mesh failure, not a partial fix** — an empty mesh is
     an unambiguous failure signal; this produces a plausible-looking
     wrong answer that only reveals itself if residuals are actually
     checked. Do not present this as a usable workaround under any
     framing. `dev/scratch/capability_survey/data/cone_asymmetric_bbox.json`,
     `.../renders/cone_asymmetric_bbox.png`.
   - `projection=:random, rng=Xoshiro(42)` (same convention as the torus
     fixture) does not produce a working decomposition either — it fails
     loudly instead (`_robust_slice_at_z` gate failure after 8 perturbed
     retries), a preferable failure *mode* to the original silent
     emptiness, but still not a working result.
     `dev/scratch/capability_survey/data/cone_projection_random.json`.
   - **Net: no known safe workaround exists for the cone specifically**,
     unlike horn torus (item 5) and the original torus fixture.
5. **Horn torus's self-tangent pinch** (`(x^2+y^2+z^2)^2-4(x^2+y^2)=0`,
   genus-1 torus with `R=r` — the hole pinches to a point at the origin):
   the same mechanism as item 4, a different degenerate point. Bare,
   default `decompose_3d_surface` returns the same completely empty
   decomposition, no exception. Not previously tested in this file — a
   new fixture confirming the existing mechanism, not a new one.
   `dev/scratch/capability_survey/data/horn_torus.json`.
   **2026-08-07 follow-up — confirmed, clean workaround**:
   `projection=:random, rng=Xoshiro(42)` (same convention as the torus
   fixture) produces a genuinely working decomposition: 20 vertices (12
   `Critical`, 8 `Artificial`), 20 edges/faces, 1647 mesh points, all
   residuals small (median `2.3e-7`, p90 `1.2e-6`, max `8.5e-4`, no wild
   outliers). Recommended, at the same tier of confidence as the existing
   torus workaround entry above.
   `dev/scratch/capability_survey/data/horn_torus_projection_random.json`,
   `.../renders/horn_torus_projection_random.png`. Separately tested:
   asymmetric `bbox_z=(-4.0, 4.3)` does **not** help this fixture — still
   returns a completely empty decomposition, no exception.
   `dev/scratch/capability_survey/data/horn_torus_asymmetric_bbox.json`.

**Falsified prediction, flagged not resolved**: the geometric argument
that "the singular point sits at the origin, fixed under any rotation,
so `projection=:random` shouldn't help" correctly predicted the cone's
outcome above but not horn torus's success under the identical
treatment — the two fixtures diverge for a reason not investigated here
(likely something about how each system's specific degree/structure
interacts with HC.jl's polyhedral start-system construction under a
generic rotation, not the geometric fixed-point argument itself, which
only speaks to where the degenerate point is, not to whether the solver
can find it).

Working hypothesis, not yet verified against HC.jl internals: `solve`'s
default settings struggle specifically when the critical-point-finding
augmented system has a solution of multiplicity ≥ 2 at the point being
searched for — so this gap sits squarely in territory this project's own
fixtures keep landing on. Worth a dedicated investigation into whether a
different `solve` configuration (e.g. disabling early path truncation, or
a dedicated multiplicity-aware follow-up step) recovers these points
reliably.

### Singular-curve decomposition (Griffis-Duffy blueprint)

**Different problem class from the item directly above — do not
conflate the two despite sitting next to each other.** The multiplicity≥2
item above is a *detection-reliability bug*: HC.jl's `solve` sometimes
fails to find critical points that demonstrably exist. This item is a
*missing capability*: HomotopyGetsReal has no first-class representation
of a singular curve as its own decomposed object at all (today, a
singular point is only `Singular`-vertex metadata when `deflate=true` —
never a retained deflated system that gets curve-decomposed in its own
right), independent of whether any particular critical point was found
correctly.

Full phased blueprint (S0–S4), architectural detail, and the live
BertiniReal/Bertini1 investigation behind it live externally at
`/Users/juancagc/bertini_migration/GriffisDuffy_BertiniReal_run_and_HGR_blueprint.md`
(outside this repo, not committed here). Summary only:

- **Core architectural insight**: deflate → retain the deflated system →
  curve-decompose it embedded (in the ambient surface's own coordinates),
  fed by singular-locus samples from surface critical-point work. Smaller
  in scope than porting Bertini's full Numerical Irreducible Decomposition
  (NID) stage; larger than just flipping a metadata flag — HGR already
  generates the geometric samples BertiniReal would otherwise read from
  `witness_data` for smooth pieces, so what's actually missing is a
  witness-like object for the *singular* branch specifically (a general
  point + the deflated equations + a projection), not a whole new solver
  stage.
- **Phase S0 (documentation only, this entry) — oracle-clarity status**:
  an independent, read-only investigation against a live Bertini
  1.7.0 + BertiniReal 1.9.0 install (nothing in this repo touched) ran the
  packaged Griffis-Duffy fixture directly. The published README numbers
  (588 candidate minors, 8 needed, full deflation sequence `4, 1, 1`) do
  **not** currently reproduce: only the **leading term (4)** was confirmed
  live (`isosingular_summary = 4 0`, "Testing for a component of dimension
  4"); the rest failed downstream of a **real bug in stock BertiniReal's
  own Python deflation script** (subfunctions emitted as `Matrix([...])`
  then used inside scalar polynomials — `TypeError: unsupported operand
  type(s) for +: 'Add' and 'MutableDenseMatrix'`), not anything to do with
  HGR. A patched local reconstruction of the Python minors found **2347**
  nonzero minors (of 2352 combinatorial), not 588 — refuted for this
  path/tooling, not merely "not yet tried." **Do not cite 588/8 as
  independently re-verified against this install without a working Matlab
  (or repaired Python) deflation that matches the historical count** —
  this repo had no prior claim to correct (grepped clean, 2026-07: no hit
  for "588", "8 needed", or the full corank sequence anywhere in
  `docs/`, `test/`, `src/`, or `dev/scratch/`), so this status is recorded
  here for the first time, not amended.
- **S1–S4** (deflated-system retention, singular-locus sampling,
  first-class embedded singular-curve decomposition, Griffis CI target)
  are **not scoped for now** — real, separate future work, most likely
  after the Albatross talk given current deadline pressure. Not
  implemented or further investigated this session.

### External reference: independent architecture/math/BertiniReal-gap audit trio (2026-07-25)

**Future-work reference only — not active work.** An independent,
zero-context review (via Cursor) synthesized three separate audits
(architecture, math/performance, BertiniReal-gap analysis) into a
full 8-phase long-term roadmap toward making HGR "the definitive
pure-Julia BertiniReal-class library." Source documents live external
to this repo (`/tmp/HGR_MASTER_CONTEXT_AND_ROADMAP.md` and its two
source audits, on the machine Cursor ran on — ask Juan for the actual
file locations before referencing them precisely; not transferred
here).

No new active-bug findings — its value here is as independent
cross-validation: it substantially corroborates this project's own
existing backlog rather than surfacing anything new. Concrete matches
confirmed against current source before writing this entry: the
duplicate `compute_critical_points` solve under `incidence=true`
(`SurfaceDecomposition.jl`'s own docstring already calls the missing
dedup with `_slab_bounds`'s internal call "a deliberately deferred
optimization"); the linear `_cells_adjacent` scan (`filter(x -> x.id
== eid, cs.edges)`-style lookups, called once per consecutive
column-pair in `_check_continuity!`, not indexed); and the
positive-dimensional-critical-loci gap this project already
investigated directly via the torus (`compute_critical_points` cannot
represent a fold locus that's a curve rather than isolated points —
see the torus entry below; resolved for that specific fixture via
`projection=:random`, but the underlying representational gap is
general, not fixture-specific).

The roadmap's own Phase 1 (singular curves as first-class geometry)
overlaps almost exactly with the already-deferred Griffis-Duffy
singular-curve blueprint above — not a second, competing plan, the
same deferred scope observed independently twice.

### Not logged as backlog: torus, resolved via `projection=:random`

Tried a torus (`(x^2+y^2+z^2+3)^2-16(x^2+y^2)=0`, hole axis aligned with
the slicing z-axis) as a new-topology validation fixture, 2026-07. Found,
and confirmed by direct calculation, a different problem from the two
items above: at the fold `z=±1`, `∂f/∂x=∂f/∂y=0` **identically for every
point** on the circle `x^2+y^2=4` there — a genuine 1-dimensional
critical locus, not isolated points, which `compute_critical_points`
(built for isolated-solution homotopy continuation) cannot represent at
all. Not logged as a backlog item because it isn't a bug in the usual
sense — the pipeline was never designed to detect positive-dimensional
critical loci, and doing so would be a real new capability, not a fix.

**Resolved, 2026-07 (second investigation, same day): `projection=:random`
closes this — an existing capability, not a new one.** The fold's
positive-dimensional locus is a coordinate-ALIGNMENT artifact (the hole
axis exactly coincides with the slicing axis), not a real singularity —
the torus itself is smooth everywhere. A generic rotation breaks that
alignment. Confirmed live across 3 seeds (`1, 7, 42`, this project's own
seed precedent): `compute_critical_z_slices` on the generically-rotated
chart always finds an isolated, nonempty 4-value critical set (never
empty), and the full `decompose_3d_surface(...; projection=:random,
incidence=true)` run succeeds (seeds `42`/`7`: 8 vertices all `Critical`
— 0 `Singular`, correctly reflecting the torus has none — 8 edges/faces,
~78k-point/~160k-triangle meshes, `|f|` residuals up to ~2.2e-5 — looser
than sphere/ellipsoid's own but not garbage — 0 degenerate triangles).
Naked-edge count after incidence stitching is real but noisy run to run
even at a fixed seed (21–22 on one run, 144/24 on an earlier run of the
same two seeds) — the same cross-process HC.jl solver jitter already
documented elsewhere in this project (e.g. Taubin's own naked-edge
spread), just a larger swing here than typically seen; not investigated
further.

Both the z-aligned failure and the `projection=:random` fix are
independently reproducible from **one** script,
`dev/scratch/scratch_torus_validation.jl` (sections 1–4: the failure,
kept as-is rather than trimmed, ~1400s and a catastrophically wrong mesh;
sections 5–6: the fix). The torus is a usable `decompose_3d_surface`
validation fixture after all.

### `verify_isosingular_dimension`'s per-attempt `solve()` defaults to `compile=:mixed`, not `compile=:none` — unconfirmed Phase 2 lead (2026-08-04)

Found while investigating why `test_isosingular_deflation.jl` dominates
full-suite runtime (see "Isosingular deflation" section's own 2026-08-04
entry for the timing trace this came from). Diagnostic only — no `src/`
change made or proposed here.

`verify_isosingular_dimension` (`src/Solver.jl:419`) builds a fresh `System`
on every retry attempt inside `for attempt in 1:cfg.isosingular_verify_retries`
and calls `solve(Faug, [x0_c]; ...)` (`:448`) with no `compile=` kwarg, so it
falls through to `HomotopyContinuation.jl`'s global default,
`COMPILE_DEFAULT[] = :mixed`. This is structurally the same shape as
`build_face_tracker` (`src/FaceTracking.jl:263-272`) — a fresh small system
built per call, "called once per anchor, potentially hundreds of times per
surface" — which already has an explicit, commented `compile = :none`
override for exactly this pattern (see "Face tracking" section, "Patch
construction: `compile=:none`..." below: benchmarked directly, `compile =
:all` cost ~5,700x what `compile = :none` cost across 15 structurally
distinct systems). `:mixed` is not `:all`, so that 5,700x figure doesn't
transfer directly — HC.jl's `:mixed` heuristic may already avoid compiling
systems this small — but the underlying pattern (many small, single-use,
freshly-built systems) is the same one that motivated the FaceTracking fix,
and `verify_isosingular_dimension` has no comment anywhere addressing the
tradeoff. Not profiled; not confirmed as a real cost here.

**Before ever applying `compile=:none` to this call**: `verify_isosingular_dimension`'s
output (verdicts, corank sequences) must be re-validated against the existing
Hauenstein-Wampler ground-truth tests (Whitney umbrella, node, cusp —
`test/test_isosingular_deflation.jl`, `test/test_solver.jl`) to rule out
`compile` mode changing *numerical* tracking behavior, not just runtime.
Compile mode is not guaranteed numerically inert in general; this hasn't been
checked for this call site specifically.

### Capability survey (2026-08-06): new gaps found across a 24-fixture survey

*(Six distinct findings from `dev/scratch/capability_survey/` — a
curve/surface capability survey across 24 fixtures at coarse test-density
config. Full data/renders/logs live there, untracked; this entry logs
only what's genuinely new, not fixture-by-fixture confirmations of
already-documented mechanisms — see
`dev/scratch/capability_survey/summary_report.md` for the complete
picture including those.)*

1. **Smooth critical points correctly located but misclassified
   `Singular`, with edge/mesh construction bypassing them entirely** —
   `squircle_quartic` (`x^4+y^4-1`) and its surface analog
   `quartic_superellipsoid` (`x^4+y^4+z^4-1`). The misclassification
   itself confirms the "Stage 4c — `_deflation_applicable`" entry below
   (genuine smooth points, e.g. `jacobian_rank=1`, `singular_values=[4,0]`
   at `(±1,0)` for the curve case — exactly that entry's own `f=x-y^3`
   illustration, now with a live measurement). The NEW part: downstream
   edge/mesh construction doesn't use these correctly-found points at
   all — it wires to off-curve `Artificial` `:endpoint_fallback` vertices
   instead (residual up to ~1.0 on the surface case — nowhere near the
   actual surface). Visually confirmed in both renders: the squircle's
   curve render is missing most of the actual curve; the superellipsoid's
   mesh has a visible wedge-shaped notch cut out near one pole.
   `dev/scratch/capability_survey/data/squircle_quartic.json`,
   `.../quartic_superellipsoid.json`, `.../renders/squircle_quartic.png`,
   `.../renders/quartic_superellipsoid.png`.

2. **Cross-process nondeterminism silently drops a whole curve segment,
   zero error signal** — `folium_descartes` (`x^3+y^3-3xy`). The general
   jitter mechanism is already documented elsewhere in this file (torus
   naked-edge spread under "Watertightness measurements"; the full-suite
   536-vs-537 flake). This specific consequence is not: 1 of 3 identical
   `decompose_1d_curve` runs (same fixture, same config) silently
   orphaned a genuine on-curve `Boundary` vertex — referenced by zero
   edges — dropping the unbounded branch's continuation segment entirely,
   with no exception. `dev/scratch/capability_survey/errors_log.md`
   ("folium_descartes" entry).

3. **`three_concurrent_lines_reducible`** (`x(x-y)(x+y)`, triple point at
   the origin), two findings:
   (a) The vertical line component (`x=0`) is structurally invisible to
   x-parametrized slicing — never appears in the edge graph at all, even
   though its two `Boundary` vertices are correctly located. No
   precedent anywhere in this file.
   (b) The triple point is found but misclassified `Artificial`, not
   `Singular` — the first concrete measurement of a case the "HC.jl
   polyhedral-solve reliability" entry above already names in the
   abstract ("reducible crossing points") but had never tested. Not a
   new mechanism — the first data point for a previously-named,
   previously-untested sub-case.
   `dev/scratch/capability_survey/data/three_concurrent_lines_reducible.json`,
   `.../renders/three_concurrent_lines_reducible.png`.

4. **`ellipsoid` post-fix residual, first measurement** — max/p99
   residual jumps to ~1.0e-05, ~170x above its own p90 (1-2 outlier
   points out of 160 mesh points). Same mechanism as the "Adaptive
   re-anchoring"/ellipsoid-discovery entry below, but that entry's own
   validation numbers are all measured on the Taubin heart — never on
   the ellipsoid fixture that originally motivated it. This is the first
   documented post-fix residual figure for the ellipsoid itself.
   `dev/scratch/capability_survey/data/ellipsoid.json`.

5. **`plot_surface_decomposition` throws `ArgumentError` on any
   completely empty mesh** — uncaught, from `_near_constant_colorrange`'s
   `extrema` call over zero points (`src/Visuals.jl:210`). Hit 3x in the
   survey: `horn_torus`, `cone`, `empty_surface`. This is a different bug
   from the `_near_constant_colorrange` entry below (Float32 round-off on
   a *non-empty* near-constant range causing full-spectrum speckle) —
   same function, unrelated failure, do not conflate the two. **Good
   candidate for an actual small `src/` fix**: an early return with an
   informative message on an empty mesh, rather than a raw crash —
   low-risk, visualization-only, doesn't require re-validation against
   any ground-truth numerical test. Flagged, not implemented in this
   pass.

6. **Wholly empty real locus handled gracefully — but wholly
   undocumented** — `empty_curve` (`x^2+y^2+1`) and `empty_surface`
   (`x^2+y^2+z^2+1`) both decompose cleanly (0 counts, no exception).
   Verified via full-text grep: no entry anywhere in this file covers a
   wholly empty real variety (the existing "empty" mentions are all about
   an empty *critical-value set* mid-pipeline — the Whitney/cone
   discussions above). Logging the graceful-handling behavior itself as
   new, positive, previously-undocumented coverage. Also: plotting
   diverges 2D vs. 3D here — `empty_curve`'s render succeeds (a valid
   blank frame); `empty_surface`'s hits the same crash as finding 5
   above.

---

## Isosingular deflation

### Stage 1 — `estimate_corank` / `_corank_plateau_hint` (`src/Solver.jl`)

**Why `estimate_corank`'s `expected_rank` has no default (2026-07).** A
prior version of this function defaulted `expected_rank` to the row-rank
convention (`length(F.expressions)`). That default was silently wrong for
column-rank callers relying on it — discovered via a verification run
whose own print statement used the default instead of matching
`deflate_once`'s internal choice. The function serves two genuinely
different, co-equal conventions with no single sensible default between
them: `intersect_bounding_object`/`_classify_vertex_type`'s "full ROW rank
== smooth point" convention (`expected_rank = length(F.expressions)`,
meaningful on a bare, possibly-underdetermined curve/surface equation)
versus `deflate_once`'s "full COLUMN rank == isolated point in ambient
space" convention (`expected_rank = length(F.variables)`, confirmed
correct against the Hauenstein-Wampler `D_det` construction — see Stage 2
below). These coincide only when `F` is square (as `Faug` always is) and
diverge for every bare, single-equation curve/surface — exactly what
deflation is called on. Every call site now states its convention
explicitly.

**The `deflation_stabilized` → `_corank_plateau_hint` rename and the
`[1,1,1]` correction (2026-07).** A prior version of this function (named
`deflation_stabilized`) returned `false` for `[1,1,1]`, treating "ends in
0" as the only success signal. That was wrong on two counts — retracted
after finding it contradicted both Hauenstein-Wampler's Definition 5.18
directly and this project's own Whitney-umbrella-handle data (a real
`deflate_once` trace holding at corank 1 for 4 consecutive rounds,
`[3,1,1,1,1]`). `[1,1,1]` now correctly returns `true`. The function was
also renamed at the same time, from `deflation_stabilized` to
`_corank_plateau_hint`, to make explicit that it is a *necessary-but-not-
sufficient* pre-filter (Hauenstein-Wampler Section 6's own explicit
warning: "a necessary condition for stabilization is that two consecutive
terms in the deflation sequence must be equal, but this is not
sufficient" — their own `f_{k,l}=[x^k,y^l]` counterexample family proves a
sequence can plateau for multiple rounds and still decrease further), not
an authoritative stabilization test — that role belongs to
`verify_isosingular_dimension` (Stage 3).

### Stage 2 — `deflate_once` (`src/Solver.jl`)

**The witness-slice-construction proposal, investigation, and retraction
(2026-07).** An earlier version of `deflate_once`'s docstring assumed `F`
must already include generic slicing hyperplanes before `deflate_once` is
usable — a "construct a witness slice first" pipeline stage was proposed
on that basis and then directly investigated. It was confirmed
unnecessary: `deflate_once` called directly on the BARE, unsliced
curve/surface equation (Whitney umbrella `x^2-y^2*z`, both the handle
`(0,0,1)` and tip `(0,0,0)`) reproduces the Hauenstein-Wampler `D_det`
operator's own published deflation sequences exactly (`{3,1,1,...}` and
`{3,2,0,...}`), using the function's own default `expected_rank =
length(F.variables)` with no slice, no witness-point construction, and no
other change. The `minorSize=2`-exceeds-1-row failure this note used to
cite as evidence a slice was required is real (that guard still fires,
correctly, at a SMOOTH point of a bare curve, `rank(J)=1`), but it does
not mean a slice is needed in general — it means `deflate_once` was being
asked a question ("is this an isolated point?") that a smooth point of a
positive-dimensional variety can never answer yes to, independent of
slicing.

### Stage 3 — `verify_isosingular_dimension` (`src/Solver.jl`)

**Rejected constructions (2026-07).** Two alternatives to "append `d`
generic real hyperplanes to the FULL, UNMODIFIED `F_current`" were tried
and rejected:
- Discarding equations before tracking/certifying: produces false
  positives, demonstrated directly on the Whitney umbrella handle.
- `R(f)=A*f` randomization: produces a demonstrated **4-of-6
  false-positive rate** on the cusp's known-non-terminal round 1
  (residuals 0.08–1.79 against the original system, on tracks that
  reported `return_code == :success`). An implementation that checked
  only the tracker's own success flag (without the separate residual
  check against `F_original`) would have wrongly verified the cusp's
  round 1 as terminal on 4 of 6 random attempts. This same measurement is
  what motivated `IsosingularVerdict` being a 3-way enum rather than a
  `Bool` — collapsing "confirmed non-terminal" and "no clean answer either
  way" into a single `false` would hide exactly this distinction.

**The attempts-count retraction (2026-07-23).** An earlier version of
`verify_isosingular_dimension`'s docstring (and, duplicated, of
`Config.jl`'s `isosingular_verify_retries` field comment) claimed every
real case resolved in exactly 1 attempt, with zero instances of a second
attempt ever being needed. This was retracted after this project's own
Ambiguous-forcing investigation measured attempts of `1,2,1,4,1,1,4,1,1,1`
across 10 fresh trials on a comparable case, and a live run separately
observed `attempts=3` (see `test_solver.jl`'s HANDLE ground-truth test and
commit `dc93320`, "Fix over-specified attempts assertion in Stage 3
verify_isosingular_dimension test"). The real, guaranteed invariant is `1
<= attempts <= cfg.isosingular_verify_retries`, not a fixed count.

### Stage 4a — `resolve_isosingular_dimension` (`src/Solver.jl`)

**The one-extra-round cost of a real (nonzero) plateau, confirmed
directly.** `_corank_plateau_hint` needs to see a repeated value to
recognize a plateau, so a genuine nonzero-limit case always costs one
deflation round beyond where the corank first reaches its true limit.
Concretely, on the Whitney umbrella handle: `corank_sequence == [3, 1,
1]`, `rounds == 2` — corank already reached its terminal value `1` after
round 1 (`F1`, 4 equations), but the loop cannot know that yet (`[3,1]`
has no repeat) and must deflate once more (`F2`, 8 equations) before
`_corank_plateau_hint([3,1,1])` finally sees the repeat and triggers
verification. This cost does not appear for a `d==0` resolution
(node/cusp/tip below) — corank hitting `0` always triggers verification on
the SAME round it's first observed, no repeat needed.

**Round counts, measured directly against every ground-truth case
available in this project, not estimated:**

| case | corank_sequence | rounds | verdict |
|---|---|---|---|
| node (bare curve) | `[2, 0]` | 1 | Resolved, dim 0 |
| cusp (bare curve) | `[2, 1, 0]` | 2 | Resolved, dim 0 |
| Whitney umbrella tip | `[3, 2, 0]` | 2 | Resolved, dim 0 |
| Whitney umbrella handle | `[3, 1, 1]` | 2 | Resolved, dim 1 |

Maximum observed: 2 rounds. `cfg.max_deflations = 10` is, on the identical
footing as `cfg.isosingular_verify_retries`, an explicit safety margin for
cases this project hasn't exercised yet (deeper singularities, real
surfaces under Stage 4c) — not a value any ground-truth case has ever
needed more than a fifth of.

**Later validation (Stage 4c, Taubin heart fixture, 2026-07-23).**
Individually inspecting all 22 real `deflate=true` firings during a
`decompose_3d_surface(...; deflate=true)` run (not just the ones surviving
to final output) found 17 of 22 at 1 round, but **5 genuine outliers at 4
rounds** (`corank_seq=[2,1,1,1,0]`), all at the slice-level critical points
`x=(±1,0)` — the heart's left/right-most z-slice extremes. All 22 still
resolved cleanly (`Resolved`, zero `Ambiguous`/`Exhausted`/`attempts>=15`),
so 4 rounds at a genuine local degeneracy is confirmed normal behavior, not
a red flag by itself — but it's a real, non-spurious deviation from the
`[1,2]` range seen everywhere else, worth knowing about before assuming a
similar count elsewhere is a bug. Separately, every one of those 22
firings resolved to `isosingular_dimension=0`, and every `d=0` resolution
has `attempts=0` by construction (`verify_isosingular_dimension` needs
zero hyperplanes when `d=0`, so its retry mechanism never runs) — the
soft-flag attempts range of `[1,4]` used during this validation was only
ever scoped to `d>=1` cases that actually exercise that retry loop;
`attempts=0` at `d=0` is not itself evidence of anything and should not be
treated as a flag in future runs.

### Stage 4c — `_deflation_applicable` (`src/Solver.jl`)

**Why this exists, and why `corank > 0` alone (the originally proposed
gate) is wrong.** Found via a real crash on the Taubin heart fixture's own
crit-slices (`ArgumentError: deflate_once: minorSize=2 exceeds available
rows`), not invented speculatively. `compute_critical_points`'s `deflate =
true` path triggers on `v_type == Singular`, which is classified against
`Faug` (the caller's pre-augmented curve system, or the auto-built
`[f,f_x,f_y]` surface system) — NOT against `F_original`'s own bare
Jacobian. A point can be `Faug`-singular (e.g. a fold w.r.t. the
x-projection, where `Faug`'s own square Jacobian is rank-deficient) while
`F_original` itself is perfectly regular there (nonzero gradient, just a
higher-order-degenerate critical point in one projection direction —
concretely reproduced by `f = x - y^3` at the origin: `f_y = -3y^2 = 0`
there, so `Faug = [f,f_y]` has a singular 2x2 Jacobian, but `f_x = 1 !=
0`, so `F_original`'s own 1x2 Jacobian has full rank 1). For a bare,
single-equation `F_original` the maximum achievable rank is 1, so under
the ambient (`expected_rank = nv`) convention `corank` can NEVER reach `0`
for such a system regardless of whether the point is genuinely singular —
confirmed directly: `corank == 1` at BOTH the fold point above and at any
ordinary smooth point of a bare curve. `corank > 0` therefore excludes
nothing; the actual crash-preventing condition is on `minor_size`, not on
`corank`'s bare sign.

**Cost, measured, not assumed.** On a case that genuinely exercises
tracking (the Whitney umbrella handle, not a case that resolves via the
free `corank == 0` shortcut), this check costs ~0.008ms against a real
`resolve_isosingular_dimension` resolution's ~1486ms — roughly **200,000x
cheaper**. Safe to call unconditionally for every `Singular`-classified
vertex; no cheaper pre-check is needed.

### The `ResolveVerdict`/`IsosingularVerdict` enum-naming-collision near-miss

*(Canonical entry — this consolidates what were previously two
near-duplicate copies of the same story, one in `HomotopyGetsReal.jl`'s
module header and one in `Solver.jl`'s `ResolveVerdict` docstring. Neither
source file owns this account; both point here.)*

`Solver.jl` is a single flat namespace — every file under `src/` is
`include`d directly into one module, with no submodule separation, so two
`@enum` blocks anywhere in `src/` can collide on a member name with no
compile error. Julia's `==` across two different enum types just silently
returns `false` rather than erroring, so a collision doesn't announce
itself; it quietly breaks whichever comparison assumed the wrong one was
in scope.

This nearly happened between `IsosingularVerdict.Inconclusive` (Stage 3)
and a first draft of `ResolveVerdict` (Stage 4a), which originally reused
`Inconclusive` verbatim from `IsosingularVerdict`. That silently shadowed
the earlier binding at module scope: every bare `Inconclusive` reference
in the file, including inside the already-shipped
`verify_isosingular_dimension`, started resolving to the *new* enum's
value instead, and `resolve_isosingular_dimension`'s own orchestration
loop's `vr.verdict == Inconclusive` comparison would have silently always
been `false` (comparing across two different enum types, which Julia
allows and just returns `false` for — no error, no warning). Caught before
ever being tested, not after — `ResolveVerdict`'s middle value was renamed
to `Ambiguous` instead of qualified, so the possibility can't recur for
this specific pair.

**Convention established going forward** (not applied retroactively to
`IsosingularVerdict`/`ResolveVerdict`'s own already-shipped, exported
members): before adding any new `@enum` to this module, grep the whole
`src/` tree for every candidate member name first, and prefix multi-value
"verdict"-shaped enums with a short tag tied to their own owning function
(e.g. `ResolveResolved`/`ResolveAmbiguous`/`ResolveExhausted`, not bare
`Resolved`/`Ambiguous`/`Exhausted`) rather than relying on a word simply
sounding unlikely to collide — that's exactly what `Inconclusive` sounded
like the first time too.

### Full-suite test timing and the 536-vs-537 assertion-count variance (2026-08-04)

Found while independently verifying `docs/ORCHESTRATOR_BRIEFING.md`'s
"537/537" test-count claim for a `CLAUDE.md` draft (not reasoned from the
doc, confirmed live: three full `HOMOTOPYGETSREAL_RUN_SLOW_TESTS=1` runs this
session).

**536 vs 537 is expected variance, not a bug.** Independently summing the 12
per-file `Test Summary` rows from a live run gives 536, one less than
`Pkg.test()`'s own reported aggregate of 537. Cause, confirmed by reading the
source, not inferred: `test/test_isosingular_deflation.jl`'s "Historical
curves (N=2)" testset (`:99-139`) loops over 4 curve/config combos —
including the astroid fixture (`f1` at `:109`, `(x²+y²-1)³+27x²y²`) — calling
`decompose_1d_curve(F, cfg; deflate=true)` live for each, then firing one
`@test verdict == Resolved` per resulting `Singular` vertex that carries
`:isosingular_verdict` metadata (`:127-131`). That count depends on how many
vertices a given live run happens to flag; it is not fixed. Matches the
paper's own §5.3 documentation of this mechanism exactly.

**Per-file timing breakdown**, live run, all 12 files,
`HOMOTOPYGETSREAL_RUN_SLOW_TESTS=1`:

| File / testset | Pass | Total | Time |
|---|---|---|---|
| Isosingular deflation (Stage 5) | 19 | 19 | 18m15.2s |
| Projection (Phase 8) | 44 | 44 | 2m01.3s |
| Solver (Phase 2) | 85 | 85 | 1m22.9s |
| Taubin heart (Phase 5 slow) | 37 | 37 | 1m08.3s |
| Visuals (Phase 6) | 24 | 24 | 41.8s |
| SurfaceDecomposition (Phase 5) | 60 | 60 | 14.0s |
| Incidence (Phase 9a) | 36 | 36 | 4.1s |
| Topology (Phase 3) | 36 | 36 | 2.0s |
| Docstring rendering | 113 | 113 | 1.9s |
| VertexRegistry (Phase 10) | 46 | 46 | 0.5s |
| Types and Config (Phase 1) | 18 | 18 | 0.3s |
| PathTracking (Phase 4) | 18 | 18 | 0.2s |

`test_isosingular_deflation.jl` alone is ~76% of full-suite runtime — more
than 3x the sum of all 11 other files combined (18m15s vs 5m42s).

**Retracted comparison, logged so it isn't re-derived**: a first pass at this
investigation compared that 18m15s whole-file figure against the file's own
`:103-104` comment (commit `0c73491`, 2026-07-24), which benchmarks its
"Historical curves" testset **alone** at 6m24.8s, and called the gap "~3x the
documented baseline." That comparison is invalid — the file has two more
slow-gated testsets not covered by the original benchmark ("Taubin heart:
..." at `:157` and `:164`) — and was retracted before reaching `CLAUDE.md`. A
fair single-testset comparison would need its own timing pass, not done here.
Separately confirmed via `git blame`: `isosingular_verify_retries` and
`max_deflations` (`HomotopyConfig`) both predate that 6m24.8s benchmark by
one day (added commit `6b4aaf7`/`4e4dbd2`, 2026-07-23; benchmark written
2026-07-24) — the config knobs were already in place when it was measured.

**The flake.** Three live full-suite attempts this session: one errored
before reaching a `Test Summary`, two passed clean. Root cause of the error
is **not confirmed** — that run's output was piped through `tail -15` before
capture, which discarded the actual error message and testset location along
with everything else but generic `Pkg.API` stack-unwind frames, and it
hasn't reproduced since (2 clean runs after it, including a dedicated
attempt to recapture it). `docs/ORCHESTRATOR_BRIEFING.md` item 6 documents a
different, already-diagnosed nondeterminism
(`MixedSubdivisions.uniform_lifting_sampler` reading Julia's global RNG,
affecting torus mesh counts) — plausible analogy, not evidence, for this
specific error. No further live reruns planned this session; next time the
full suite runs for a real reason (not a diagnostic one), capture its output
via direct redirection to a file, never through a pipe that can silently
truncate on failure — the pattern this session eventually settled on for the
two runs that did land cleanly.

See "Backlog: production robustness gaps" for a `compile=:none` lead found
during this same investigation.

**Follow-up, 2026-08-06**: a full/fast suite run can also crash with a
native segfault (signal 11) inside GLMakie's window-creation path
(`test_visuals.jl`), when the machine's display has gone to sleep/locked
while the system itself stays awake — confirmed via `pmset -g assertions`
showing `UserIsActive: 0`, and conclusively isolated from any
source-code cause via a pristine-file control test (the completely
unmodified `src/Visuals.jl`, zero code changes, reproduced the identical
segfault in the same environment). This is a **plausible, not
confirmed**, explanation for at least some of the previously-unconfirmed
full-suite crash instances documented above — the original failed run's
trace was lost to the `tail`-piping mistake described above and was
never actually compared against this specific mechanism, so this is a
new lead, not a retroactive proof. Practical implication: an unattended
full-suite run (e.g. a scheduled Routine or an overnight run) should
keep the display awake, or expect this as a known failure mode unrelated
to code correctness.

---

## Face tracking (Phase 5) — `src/FaceTracking.jl`

### Patch construction: `compile=:none` over parametric `ParameterHomotopy`

Patch equations are built with literal `x0,y0,a,b` coefficients baked in
per anchor, via `build_tracker(...; compile = :none)`, rather than a
single parametric patch with `x0,y0,a,b` threaded through as extra
`ParameterHomotopy` parameters. This was benchmarked directly (Phase 5):
across 15 structurally distinct systems, `compile = :all` cost ~5,700x
what `compile = :none` cost, dominated by Julia re-specializing `track`'s
call graph against each new `CompiledSystem` type; `:none`'s per-step
tracking throughput was not worse for the small (1-2 equation) systems
this pipeline builds. See `build_tracker`'s own docstring
(`src/PathTracking.jl`) for the general `compile` tradeoff this choice
draws on.

### Adaptive re-anchoring

**The ellipsoid discovery (2026-07).** A fixed literal patch line is only
guaranteed to keep intersecting the level curve as `z` sweeps away from
the anchor when the surface's gradient is radially symmetric about the
axis the curve shrinks toward — true for a sphere (an algebraic
coincidence of radial symmetry, not a general guarantee), false in
general. This was found empirically while validating this file against
an asymmetric ellipsoid (`x^2+4y^2+9z^2=1`, `scratch_phase5_check.jl`
section 7): for a general surface the fixed line can lose transversal
intersection with the curve entirely before reaching `z_bottom`/`z_top`,
which no amount of `_track_path_segment!`-style bisection can fix, since
the system genuinely has no nearby solution beyond that point rather than
just a hard-to-resolve one.

**Why both gates, confirmed empirically neither alone was sufficient.**
Tracking anchor `(-0.722,0.346)` toward the ellipsoid's pole at `z=-1/3`:
the transversality-drift (cosine) check alone re-anchors proactively, but
a coarse `z_targets` grid (e.g. `midslice_sample_density = 8`) can still
take one large single hop that crosses the entire remaining
transversality margin in one step — measured with 8 hops toward `z=-1/3`,
`cos_angle` was `0.72` after hop 6 (above even a generous re-anchor
threshold), yet hop 7 already landed with `residual ≈ 0.076`. The
residual-based bisection gate alone is also insufficient on its own,
since without proactive re-anchoring a fixed patch line can lose
transversal intersection entirely (see the ellipsoid discovery above),
a failure the residual check has nothing to bisect against — there is no
nearby solution to converge toward. Both mechanisms are needed together;
see `_sweep_hop!`'s own docstring for the current two-gate mechanism.

**The `_gradient_at` sign-error near-miss.** `_sweep_hop!` deliberately
recomputes the anchor's raw gradient `(fx0,fy0)` via a direct
`_gradient_at` call rather than reconstructing it from
`patch_direction`'s `(a,b)=(fy0,-fx0)` return. This was found the hard
way during validation: an earlier version reconstructed `(fx0,fy0)` by
inverting `patch_direction`'s swap-and-negate, got the inversion wrong,
and the bug was silent — it just negated the computed `cos_angle`,
causing the transversality check to fire on literally every hop instead
of only when genuinely needed, with no exception or test failure to flag
it. The current direct-call rule (stated in `_sweep_hop!`'s own
docstring) exists specifically to make that mistake structurally
impossible to reintroduce by accident.

---

## Surface decomposition, Phases 8-9 (`src/SurfaceDecomposition.jl`)

### Robust z-mid selection: rejected alternatives and gate calibration — `_robust_slice_at_z`

`_robust_slice_at_z` chooses which literal `z` to hand to `slice_at_z`
for one slab, defending against a naive midpoint that coincidentally
lands on a NON-REDUCED (repeated-factor) plane curve. The failure was
discovered while validating this file against the Taubin heart surface
(`(x^2+(1.2y)^2+z^2-1)^3 - x^2z^3 - 0.1(1.2y)^2z^3`,
`scratch_phase5_taubin_check.jl` section 6); on the Taubin heart's
degenerate slab, retry attempt 1 already succeeds (see the retry report
in `scratch_phase5_taubin_check.jl`).

**A direct residual-magnitude tiebreaker was rejected.** Flagging `z_mid`
suspect by `maximum(|_residual_at(...)|)` over every edge sample point,
instead of inspecting vertex types, was empirically tested: raw
`sample_edge` output is linearly interpolated between homotopy-tracked
points, so its residual has a substantial baseline even on healthy
slices — a fully clean control slice at `z=0.05` already measures max
raw residual `0.105`; the fully clean `[1.0648,1.2367]` slab measures
`0.144`; the narrow-but-fine `[1.0,1.0648]` slab measures `0.268`; the
genuinely-degenerate `z_mid=0` slab measures `1.000`. There is no
scale-free gap to split a threshold into — the "fine" cases already span
nearly 3x among themselves, and the "bad" case is only ~3.7x above the
worst "fine" case. The `Singular`-typed-vertex co-occurrence check the
function actually uses has no such scale problem (it is boolean, not a
magnitude) and correctly separates all four cases.

**Why the gradient gate exists at all.** Confirmed by re-running the
Taubin heart regression once the `Singular` co-occurrence refinement was
in place: the `[-1,1]` slab's retry landed on a topologically CLEAN
`z_mid=0.02` (2 edges, 0 `Artificial`/`Singular` vertices) — yet the
downstream sweep was STILL catastrophic (max `|f|` via `track_face`:
`1.43`). Root cause: near `z=0`, every first partial derivative of `f`
scales like `O(z^2)` together, so no vertex ever gets classified
`Singular` even though the patch system is nearly singular there.

**A same-point gradient ratio was rejected in favor of the cross-z
reference ratio the function actually uses.** `hypot(patch_direction(...))
/ |f_z|` at the candidate's own anchor does not correlate with sweep
quality: the two already-healthy reference slabs measured `0.60` and
`1.16`, while already-confirmed-bad candidates near `z=0` measured
`2.08`-`4.29` — HIGHER than the healthy baseline, because `f_z` is
suppressed by the same `O(z^2)` factor as `f_x`/`f_y` near `z=0`, so the
ratio between them stays `O(1)` right through the degenerate
neighborhood. The adopted cross-z reference ratio was measured instead:
the two already-healthy reference slabs scored `0.82`/`0.98`; the
`[-1,1]` slab's bad candidates scored `0.0014`-`0.0057`; its
eventually-accepted good candidate (`z=0.06`, independently confirmed via
`track_face` to give max `|f|=2.4e-6`) scored `0.013`. The chosen
default, `cfg.z_mid_gradient_ratio_tol = 0.01`, sits with almost two
orders of magnitude of margin on both sides of this gap.

**Why the reference scale is computed eagerly, not lazily.** Measured
2026-07 (see
`dev/scratch/scratch_robust_slice_eagerness_check.jl`): with
`bbox_z = (-0.96, 1.3)` (an ordinary asymmetric bounding-box crop), the
bottom slab's naive midpoint `z=0.02` is topologically clean (no
`:endpoint_fallback`, no `Singular`) yet its downstream sweep measures
max `|f| ≈ 1.6`, versus `2.4e-7` at the gradient-gate-chosen
`z ≈ 0.059`. A retry-armed LAZY variant (gradient gate active only from
the first retry onward) was evaluated and REJECTED for exactly this
reason: it silently accepts that slab. The measured eager cost is a ~3x
multiplier on a healthy slab's slicing time (sphere/ellipsoid: 3.0x),
accepted as the price of the guarantee.

**The min-vs-max false positive and its Phase 8 resolution.** An earlier
version of the gradient gate compared the candidate's MIN anchor
gradient against the reference's MAX, which false-fires on curves with
multiple, legitimately very-differently-conditioned branches. On the
fixed-axis `[1.0, 1.0648]` slab this cost one harmless avoidable retry
(candidate min `0.035` vs reference max `3.94`, ratio `0.0089`, just
under the `0.01` threshold). The Phase 8 five-seed rotated-Taubin
regression then showed the FATAL form of the same artifact: narrow slabs
between two close genuine critical values (seed 3: `[-0.864, -0.858]`,
seed 4: `[0.884, 0.900]`, nowhere near the singular ellipse) carry a
just-born tiny branch whose Critical anchor gradient (`4e-4`/`3.5e-3`) is
structurally weak across the ENTIRE slab, so every retry candidate failed
identically and the ladder exhausted. The structural-heterogeneity skip
(the function's current gradient-gate refinement) resolves both:
measured `ref_min/ref_max` is `0.00306` on the `[1.0, 1.0648]` slab (gate
skipped, false positive gone) and `0.372`/`0.0723` on the `[-1, 1]`/
`[1.0648, 1.2367]` slabs (gate active, all true positives preserved).

**Phase 8 transversal-singular-curve measurements.** On the rotated
Taubin heart (seed 1), the naive midpoint of slab `[-0.956, -0.24]` had 2
`Singular` vertices (the ellipse crossings), 2 benign
`:endpoint_fallback` vertices, and per-edge first-sample anchor gradients
`{4.01, 4.01, 1.44, 8.6e-15, 1.8e-13, 0.28}` against reference scale
`3.98` — BOTH gates false-fired on every retry candidate before the
topology-gate and gradient-gate-exclusion refinements existed, so 5 of 9
slabs threw. For the gradient-gate exclusion refinement, a mid-edge
anchor was evaluated as an alternative and REJECTED with data:
`sample_edge` chords sit measurably off-curve precisely in the degenerate
cases (the z=0.02 catastrophe candidate reads a healthy `0.353` at its
chord midpoint while its downstream sweep is provably catastrophic),
which would blind the gate exactly where it is most needed.

**Validation of the refined gates** (measured 2026-07): fixed-axis
Taubin — `[-1,1] -> 0.06` retried exactly as before (reference
heterogeneity `0.372`), `[1.0, 1.0648] ->` naive midpoint (heterogeneity
`0.00306`, resolving the documented false positive), `[1.0648, 1.2367]
->` naive midpoint (heterogeneity `0.0723`); full fixed-axis decompose
max `|f| = 2.4e-6`. Rotated Taubin, seeds 1-5: zero retries, zero throws;
seed 3 (previously fatal) full decompose median world-`|f|` `2.1e-8` with
21/5677 points above `1e-4`, confined to the singular-curve band.

### Slab boundary merging — `_slab_bounds`

`_slab_bounds` merges z-slab boundaries closer than `cfg.min_slab_width`
because path endpoints landing on a point singularity carry
~accuracy^(1/multiplicity) scatter in their z-estimate: ~2e-4 observed on
the rotated Taubin heart, beyond `vertex_match_tol`'s own reach, which
would otherwise mint sliver slabs centered ON a singular point that no
`z_mid` choice can slice. All existing fixed-axis fixtures have
critical-z gaps >= 0.065 — 65x the default `min_slab_width` floor of
1e-3 — so none of their bounds are affected by the merge.

### Phase 9a/9b: crit-slice incidence and the `critical_point` coincidence fix — `SurfaceIncidence`, `_cells_adjacent`, `_decompose_crit_slice`

*(Canonical entry — consolidates what were previously two near-duplicate
copies of the same investigation, one in `SurfaceIncidence`'s docstring
and one in `_cells_adjacent`'s. Neither source site owns this account;
both point here.)*

`SurfaceIncidence.continuity_ok` flags, per face, whether every pair of
CONFIDENT consecutive boundary columns lands on combinatorially adjacent
crit-slice cells (`_cells_adjacent`). A direct diagnostic measurement
(2026-07) — ground-truth residual/geometry checks, resolution-sensitivity
re-runs at higher `edge_sample_density`/`midslice_sample_density`, and
direct coordinate comparison of the flagged cells themselves, see
`dev/scratch/scratch_continuity_ok_diagnosis.jl` — found that the
ORIGINAL measured firing rate (fixed-axis 10/14 faces `continuity_ok`/16
flagged pairs, rotated seed 1 15/22 faces/17 flagged pairs) was dominated
by a FALSE POSITIVE, not genuine branch-jumping or resolution-limited
ambiguity: 12 of the 16 fixed-axis pairs (and a comparable share of the
rotated ones) paired a `:critical_point` landing against a
`:crit_slice_vertex` landing that are, by construction, the SAME physical
fold-tip location — confirmed coincident to machine precision (~1e-16) by
direct coordinate comparison — which `_cells_adjacent`'s original rule
had no case for (only exact `(kind, id)` matches or `:edge`-involving
pairs were ever considered adjacent).

The tell that distinguished this from genuine ambiguity: these 12 pairs
did NOT shrink at higher sampling density the way the OTHER (genuine) 4
fixed-axis violations at z=1.0648 did (which vanish completely by
`edge_sample_density = 16`) — they persisted and grew MORE confidently
separated, the opposite of what resolution-limited ambiguity should do.

Fixed by adding `:critical_point` coincidence cases to `_cells_adjacent`
(coincidence checked against `cfg.vertex_match_tol`, the same tolerance
`weld_mesh`'s own clustering uses): a `:critical_point` and a
`:crit_slice_vertex`/edge-endpoint that are the SAME physical location
within that tolerance now count as adjacent. Re-measured after the fix:
fixed-axis drops to 4/14 faces flagged with exactly the 4 genuine
z=1.0648 pairs remaining (0 at the fold-tip boundaries, confirmed still
resolving to 0 at `edge_sample_density = 16`, unchanged); rotated seed 1
drops to roughly 5-7/22 faces (~10-11 pairs, run-to-run jitter from the
same HC.jl cross-process nondeterminism documented elsewhere in this
codebase). The residual flags on both fixtures now concentrate at the
SAME multi-face edge-type boundaries the watertightness measurements
below document as an open coverage-gap issue — this is exactly the
evidence-gathering `continuity_ok` exists for, not a surprise finding.

Separately, `_decompose_crit_slice` was measured across every fixture
(fixed-axis and rotated Taubin included): nodal/saddle-type crit-slices
decompose like the nodal-cubic fixture, fold-type extremes come back
empty or as isolated `Singular` vertices, and nothing throws.

`_inward_row_points` exists because the boundary row's own inter-column
spacing self-referentially defeats a ratio check near a converging fold
point. Measured on the fixed-axis Taubin cusp: boundary-row spacing
`0.00355` vs a landing distance of `0.0037` (fails a `1.5x` ratio against
itself) while the row-inward spacing (the SAME columns, one ring earlier)
is `0.0982` (comfortably passes).

### Phase 9b: monotone snap-target assignment — `_monotone_snap_targets`

Motivated directly by measurement: plain nearest-target assignment (no
monotonicity constraint) was checked on the fixed-axis Taubin `z=1.0`
boundary and visits targets non-monotonically on more than half of
same-edge column runs (4/8 sides monotone without this constraint, 16/56
segments spanning >= 2 targets) — which would zigzag the welded boundary.
The two-pass monotone DP eliminates that by construction.

### Phase 9c: edge chaining, loft triangulation, and the Option A investigation — `_chained_edge_polylines`, `weld_mesh`

`_chained_edge_polylines` chains adjacent crit-slice edges before
`_split_t_junctions` runs, because a T-junction at a point shared between
two DIFFERENT crit-slice edges is otherwise invisible to it
(`_split_t_junctions` only fan-splits within a single polyline). Confirmed
empirically (2026-07) to be a real, common case: 10 of 16 residual naked
edges at the fixed-axis Taubin z=1.0/z=1.0648 boundaries (measured with
per-edge-only polylines) had endpoints on DIFFERENT edge_ids; chaining
closed roughly half of those.

`_append_loft_triangles!`'s seam-capping triangle is confirmed empirically
(2026-07) to close most but not all such seams; the residual is part of
the watertightness measurement below, not silently treated as fully
solved.

`weld_mesh`'s three-mechanism `incidence` path (snap, loft, chain+split)
was compared directly (2026-07) against a HIGHER-fidelity but
higher-risk alternative: BertiniReal-style targeted homotopy tracking
toward known crit-slice vertices, replacing free-sweep `track_face` hops
("Option A"). Measured precision showed no mechanism for Option A to beat
direct reuse of already-Newton-polished crit-slice coordinates (parity at
best, given it necessarily adds tracker predictor-corrector error on top
of converging to the SAME target), so it was not pursued. The residual
naked edges after all three mechanisms (see the watertightness
measurement below) break down as roughly a third genuine
cross-edge-junction cases the chaining did not resolve, a third seam
artifacts at run boundaries the capping did not fully reach, and a third
UNDIAGNOSED after three separate investigation attempts — confirmed NOT a
resolution artifact (coarser/finer `edge_sample_density` does not trend
this to zero). Full closure is DEFERRED as a future, separately-scoped
phase.

**Checked directly, 2026-07 (fourth investigation attempt): the
`sample_edge` straight-chord bug fixed below is NOT the explanation for
the undiagnosed third.** Before the fix, cross-referencing residual naked
edges against their nearest crit-slice edge found a plausible partial
link at `z=1.0` (naked edges landing 0.08-0.17 away from severely
under-sampled reference points — see `sample_edge`'s own docstring for
the fix). After the fix (`sample_edge` now Newton-projects every
interpolated point onto the true curve, residuals down from up to 0.28 to
~5e-7), the naked-edge count was re-measured on the identical fixed-axis
Taubin fixture and found **unchanged**: 188 unstitched / 32-34 with full
`incidence=true` (previously 31-35) — within the same cross-process
jitter range, not a reduction. The plausible link the cross-reference
suggested did not materialize into a measurable improvement once the
underlying point accuracy was actually fixed; recording this so a future
reader doesn't re-investigate the same already-ruled-out hypothesis. The
undiagnosed third remains undiagnosed.

### Watertightness measurements

*(Canonical entry — consolidates what were previously three
near-duplicate tellings of the same progression, in `weld_mesh`,
`_naked_mesh_edges`, and `decompose_3d_surface`'s `incidence` keyword.
None of the three source sites owns this account; all point here.)*

`_naked_mesh_edges` counts undirected mesh edges belonging to exactly one
triangle — a surface closed within the bounding box should have none.
Measured on the fixed-axis Taubin fixture (verified deterministic across
sessions): the plain tolerance-only `weld_mesh` (`incidence = nothing`)
leaves 188 naked edges (10 at z=-1, 70 at z=1.0, 84 at z=1.0648, 24 at
z=1.2367). Phase 9b's incidence-based snap-stitching
(`weld_mesh(...; incidence = ...)`, step 1 only) closes this to 58
(0/26/32/0) — the fold/point-type boundaries (cusp, tips) close
COMPLETELY. Phase 9c's coordinated loft (steps 2-3, layered on the same
`incidence` path) reduces it further to 31-35 across repeated decomposes
(the spread is cross-process HC.jl solver jitter, not nondeterminism in
the loft logic itself; see `test/test_taubin.jl`'s Phase 9c section) — a
further ~40% reduction over 9b, all of it at the two multi-edge,
multi-face EDGE-type boundaries (the singular notch, the saddle pair).
The 188 and 58 baselines are asserted exactly in `test/test_taubin.jl`
(fully deterministic); the Phase 9c result uses loose bounds there for
the documented jitter reason, not a weaker guarantee about the mechanism.

Measured cost of `decompose_3d_surface(...; incidence = true)` on the
fixed-axis Taubin fixture: ~+9 s on a ~16 s decompose (4 crit-slices plus
one extra critical-point solve for fold anchors; the snapping/splitting
mechanisms themselves are comparatively negligible, no new solves).

### Density-dependence and cross-fixture extension of the Watertightness gap (2026-08-06 capability survey)

Follow-up to "Watertightness measurements" above, from
`dev/scratch/capability_survey/`'s 24-fixture survey and its targeted
render-gap follow-ups (`investigate_render_gaps.jl`,
`investigate_taubin_seam.jl`, `investigate_taubin_seam_matched_zoom.jl`,
`investigate_taubin_density_choice.jl` — full logs/renders there, not
duplicated here).

**Residual anomaly that motivated this investigation**: `taubin_heart`'s
max residual at this survey's coarse density (`2.85e-05`) is ~12x the
existing fixed-axis production-density baseline stated elsewhere in this
file (`2.4e-6`, see "Validation of the refined gates" above). This gap
remains **unconfirmed** — nothing below resolves it. A plausible but
explicitly unconfirmed contributing factor: at coarse density,
`taubin_heart`'s bare mesh has substantially more unstitched connectivity
than production-density validation typically shows (see the naked-edge
data immediately below) — logged as a plausible link, not a settled
explanation.

Naked-edge counts (`_naked_mesh_edges`, `incidence=true` in every case):

| Config | torus naked edges | taubin_heart naked edges |
|---|---|---|
| Coarse survey (`edge_sample_density=6`, `midslice=8`) | 8 (from 80 bare) | 29 (from 132 bare) |
| Paper-figure config (`edge=8`, `midslice=8`; taubin only) | — | ~30-31 |
| Full production defaults (`edge=50`, `midslice=100`) | — | ~56-63 |

(a) `incidence=true` substantially reduces naked edges for both fixtures
at every density tested, but does not fully close them at any of
them — matches this file's own "reduces... without fully closing" caveat
above exactly, now confirmed across multiple density regimes rather than
only the one production-density measurement the original entry used.

(b) **New, not covered by the entry above**: sampling density is a
second, independent lever on this limitation's *visual* severity,
distinct from the documented lever (pipeline stage: bare → `incidence=true`
→ Phase 9c loft). Confirmed via matched-axis-limit rendering (controlling
for plot zoom/bbox as a confound — the first comparison attempt here
didn't control for this and had to be redone): the same underlying gap
looks markedly smaller at higher density even though the raw naked-edge
count rises (more, individually shorter perimeter edges sampling a
fixed-or-shrinking gap area at finer resolution — not a contradiction).

(c) **New, and not merely an internal robustness note**: the actual
archived, already-published paper figure
(`paper_artifacts/figures/taubin_singular_structure.pdf`, the §5
Validation incidence-overlay figure) was generated at
`edge_sample_density=8` — essentially the same density as this survey's
own coarse config — and exhibited this same defect, confirmed by direct
visual inspection of the figure's own preview PNG (partially obscured by
its chosen camera angle and vertex-marker overlay, but genuinely
present). **Resolved, 2026-08-06**: the intermediate density
(`edge_sample_density=20`, `midslice_sample_density=25` — reviewed
against the current `edge=8` figure and against full production
defaults, `edge=50`/`midslice=100` at ~138MB, impractically heavy) was
accepted. `paper_artifacts/figures/taubin_singular_structure.pdf` and
`taubin_talk_hook.pdf` were both regenerated in place at this density
(originals preserved as `*_PRE_DENSITY_FIX_2026-08-01.*`, not deleted),
and `paper_artifacts/taubin_singular_structure_example.jl`'s own density
defaults were updated to match — this is now the reproducible default
for both figures, not a one-off scratch render. As of this writing, the
regenerated files still need to be manually uploaded into the paper's
and talk's Cowork sessions by Juan — this repo-local change doesn't
propagate there automatically.

(d) **Methodological note, useful beyond this survey**: raw
`_naked_mesh_edges` counts are not comparable across fixture types
without care. A bbox-clipped rim on an *unbounded* surface (e.g.
`hyperboloid_one_sheet`: 196 naked edges at default density, purely from
its clipped boundary) is naked *by construction* — there is no surface
beyond the box to stitch to — unlike genuine interior stitching cracks on
a closed/near-closed surface (torus, taubin_heart). Separately confirmed:
`hyperboloid_one_sheet`'s own visual jaggedness at coarse density is a
sampling-density artifact (smooths out almost completely at default
density), not the naked-edge/stitching mechanism above — a different
fixture, a different cause, do not conflate the two.

### Known limitation: generic projections over singular curves — `decompose_3d_surface`

Measured (rotated Taubin heart, seed 1, 2026-07): median world-`|f|`
residual `6e-9` over 2617 mesh points, but ~2.5% of points exceed `1e-4`
(max `0.26`), ALL confined to `|z_world| <= 0.14` around the singular
plane (the surface spans `|z| <= 1.24`). Away from the singular locus the
decomposition is unaffected.


---

## Generic projection support (Phase 8) — `src/Projection.jl`

*(See also "Known limitation: generic projections over singular curves"
under "Surface decomposition, Phases 8-9" — that entry is a DIFFERENT
investigation, mesh-quality degradation near a singular locus AFTER a
projection is applied, not this section's construction/validation of the
projection itself.)*

### `random_orthogonal_matrix`'s Haar-uniformity verification

Plain LAPACK `qr()` of a Gaussian matrix is not itself Haar-distributed;
the sign-correction step (multiplying columns by the signs of `diag(R)`)
is what restores Haar measure on O(n). Verified empirically (20k
samples): det always +1, `||mean(column)||` consistent with the
1/sqrt(N) null, `E[(col)_z^2] = 0.3323` vs 1/3 for uniform.

### `projection_orthonormality_tol` — the Audit 1 Item 4/2a fix

*(Canonical entry — consolidates what were previously two near-duplicate
copies of the same fix history, one in `Config.jl`'s
`projection_orthonormality_tol` field comment and one in
`_resolve_projection`'s own docstring. Neither source site owns this
account; both point here.)*

Before 2026-07-23 (Audit 1, Item 4/2a), the orthonormality acceptance
threshold inside `_resolve_projection` (`norm(Q'Q - I) <= tol`) was a
bare `1e-8` Float64 literal, with no `HomotopyConfig` field to source it
from at all — a gap against this project's own "every numerical
tolerance threaded through `HomotopyConfig`" commitment (see
`Config.jl`'s own header). Fixed by adding
`cfg.projection_orthonormality_tol`, preserving the exact
already-working `1e-8` default while making it genuinely configurable.
Deliberately its own field rather than a reuse of `jacobian_rank_tol`
(same default magnitude at the time, but a different physical quantity —
a matrix-identity defect, not a Jacobian singular-value cutoff).

### `_verify_projection_ok`'s degeneracy-discrimination measurements

Measured: `z - x^2` with `Q = I` flagged (probe max exactly 0), a 1e-9
near-degenerate rotation flagged (ratio 6e-10), a 1e-7 rotation passes,
and the sphere / Taubin heart / cylinder all pass.

## Plotting (Phase 6) — `src/Visuals.jl`

GLMakie was confirmed (empirically, not assumed) to render headlessly in
this project's dev sandbox before this file was written — no CairoMakie
fallback needed, `Project.toml`'s existing `GLMakie` dependency
(already listed, never previously `using`d) was sufficient.

### `_near_constant_colorrange`'s Phase 6 "02b" discovery

Found via the Phase 6 "02b" investigation: a mathematically-exact-1.0
quantity's Float32 round-off, ~6e-8 relative, auto-scaled a colorbar
into full-spectrum speckle — confirmed via `radial_fn`'s printed
min/max on the unit sphere: literally `0.99999994` to `1.0`.
