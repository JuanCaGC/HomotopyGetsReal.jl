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

### Not logged as backlog (informational only)

Tried a torus (`(x^2+y^2+z^2+3)^2-16(x^2+y^2)=0`, hole axis aligned with
the slicing z-axis) as a new-topology validation fixture, 2026-07. Found,
and confirmed by direct calculation, a different problem from the two
above: at the fold `z=±1`, `∂f/∂x=∂f/∂y=0` **identically for every point**
on the circle `x^2+y^2=4` there — a genuine 1-dimensional critical
locus, not isolated points, which `compute_critical_points` (built for
isolated-solution homotopy continuation) cannot represent at all. Not
logged as a backlog item because it isn't a bug in the usual sense — the
pipeline was never designed to detect positive-dimensional critical
loci, and doing so would be a real new capability, not a fix. Reorienting
the torus so its hole axis is *not* aligned with the slicing axis
(matching `prototipo_viejo_julia/HomotopyGetsReal.jl`'s own original
orientation) avoids this specific problem — confirmed live that
`compute_critical_z_slices` then finds 4 clean, isolated critical
z-values (`[-3,-1,1,3]`) — but the full `decompose_3d_surface` validation
run against that corrected orientation was not completed this session
(the first, wrong orientation alone cost 1419s and produced a
catastrophically wrong mesh before this was diagnosed).

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
