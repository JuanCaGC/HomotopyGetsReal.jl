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
