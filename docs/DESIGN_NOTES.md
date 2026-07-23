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
