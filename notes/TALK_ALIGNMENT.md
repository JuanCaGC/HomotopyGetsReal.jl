# Talk deck vs manuscript alignment

**Sources read, not edited:** `paper_artifacts/figures/src/talk.tex`
(current 25-slide Albatross deck) and
`paper/HomotopyGetsReal_paper_current.md`. Supporting context from
`notes/REPO_TECHNICAL_SUMMARY.md` and `notes/TALK_PREP.md` only where
needed to decide whether a slide claim is paper-backed.

The lean 23-slide file `talk_lean_23slides.tex` is **not** the current
deck; it is behind on test count (537) and version (0.2.0). This note
is about `talk.tex`.

---

## Narrative overlay

| | Paper | Talk |
|---|---|---|
| Spine | Architecture → examples → residuals → Bertini_real | Hook → gap → pipeline → **astroid headline** → Taubin walkthrough → projection → tests → comparison → Whitney |
| Headline result | Two extensions (generic projection + diagnostic deflation), residual table | Astroid 4-vs-6, then Taubin 14-vs-20 |
| Curves | Buried in §3.6 and §6 | First empirical block (slides 9–12) |
| Certification | Not named | Slide 8 names `certify` and says it is unused |
| Exact methods | CAD in §1 only | Dedicated msolve / AlgebraicSolving.jl slide (audience-specific) |

This is a defensible talk structure for Albatross (exact-methods room,
Rémi Prébet named on slide 6). It is **not** a slide-for-slide of the
paper. Core theorems and definitions that the paper actually has *are*
present, but in a different order and at different depth.

---

## What aligns

- **Problem statement** (slide 4) matches the paper’s opening: cell
  decomposition of real points via homotopy, not CAD.
- **Gap** (slide 5): Bertini_real needs a witness set / NID; HGR takes
  the bare equation. Matches §2.4 and §6.
- **Pipeline diagram** (slide 7): curve `{f, ∂_y f}` and surface
  `{f, ∂_x f, ∂_y f}` match the **code** and the technical summary.
  They do **not** match paper §2.3’s written example (see
  `notes/MANUSCRIPT_REVIEW.md`). On this point the talk is the one that
  is correct.
- **Astroid 4 vs 6** (slides 9–10): same finding as paper §6, including
  Amethyst–Hauenstein–Wampler as third confirmation, and “finer
  partition, not missing topology.”
- **Taubin plain vs incidence** (slides 17–18): 14 vertices / 0
  Singular vs 20 / 2 Singular. Matches paper §4.2 and
  `results.json` overlay counts.
- **Generic projection / torus** (slide 19): matches paper §4.3 / §6
  (fixed-axis fails completely; `projection=:random` is the fix).
- **Whitney sequences** (slide 22): `[3,2,0]` tip, `[3,1,1]` handle,
  round-5 minor explosion. Matches paper §3.5 and the tests.
- **Deflation diagnostic-only** (slides 22–23) matches paper §3.5
  Scope and §6.
- **Griffis–Duffy** (slide 23): leading term only, Bertini_real tooling
  bug. Matches paper §6 citation caution.
- **Availability** (slide 24): v0.2.1, MIT, `Pkg.add`. Matches §7.

---

## Discrepancies (paper-backed)

Severity: **must fix before the talk** / **should fix** / **framing**.

### Must fix

1. **Test count: talk 538/538 vs paper 537/537.**
   Slide 20 prints `538 / 538`. Paper §5.3 prints `537 / 537` and
   documents ±1 jitter from the astroid deflation loop. The technical
   summary (`REPO_TECHNICAL_SUMMARY.md:675`) already lists 536, 537,
   and 538 in different files. Putting **538** on a full-bleed slide
   picks the high end of a known unstable count. Either use `537 ± 1`
   as the paper does, or pin a live run the morning of the talk and
   put *that* number on the slide.

2. **“One of the three validation surfaces” (slide 19).**
   The same slide lists four surfaces (sphere, ellipsoid, Taubin,
   torus) then says the torus is “the only one of the **three**
   validation surfaces” for which the default fails completely. Paper
   §6 says “the only one of the **four**.” Use four.

3. **Naked-edge drop 188 → 34 (slide 16) vs paper 31–35 range.**
   `results.json` `pipeline_steps_surface.stitching.n_naked_after` is
   34, so the slide matches **that run**. Paper §3.9 reports 31–35
   across repeated decomposes and refuses a single count. A single
   number on a slide is fine if spoken as “this run”; it is not fine
   as an invariant. Add “this run” or “~30–35.”

### Should fix (accuracy vs paper)

4. **Talk never shows the residual table (paper Table 3 / Figure 5).**
   Slide 8 defines residuals; no numbers. The paper’s actual
   quantitative claim (sphere mean 2.85×10⁻⁸, ellipsoid max 1.12×10⁻⁴
   above the test gate, torus looser and characteristic) is absent.
   For an exact-methods audience this is the slide most likely to be
   asked. Consider one compact Table-3 slide, or a single residual
   histogram (`residual_histograms.pdf` exists and is unused).

5. **`sample_edge` regression (paper §5.1, Figure 6) is not on the
   current deck.** It is on `talk_lean_23slides.tex` slide 18
   (`08_sample_edge_fix_before_after.pdf`). The paper treats this as
   the concrete illustration of what residual testing is *for*. The
   current talk dropped it. Restoring one slide would back slide 8
   with a number (0.4998 → 7.3×10⁻⁷) that exists in `results.json`.

6. **Ellipsoid is named, never shown.** Slide 19 lists it as a
   genus-zero control. The paper’s reason it exists (asymmetric, cannot
   hide a swapped axis; production max residual **above** the 10⁻⁴ test
   gate) never appears. One sentence on slide 19, or a residual-table
   row, would be enough.

7. **Sphere is absent.** Paper §4.1 is the smoothness control
   (critical z = ±1, 19,504 vertices). The talk jumps from architecture
   to astroid. Fine for time; not aligned with the paper’s example
   order.

8. **Future-work wording vs paper/code.** Slide 23: deflation “does not
   yet **retain** the deflated system.” Paper §3.5 and
   `Solver.jl` (as summarized in the technical notes) **do retain**
   `F_final`; they do not **consume** it. Say “does not use,” not
   “does not retain.”

9. **“Curves, surface, and singular structure: same pipeline”
   (slide 3).** True at the level of shared tracking / config. Easy to
   hear as “singular curves are first-class output.” Paper §4.2 / §6
   spends a paragraph saying the opposite for the default API. Slide
   17 exists to correct this; keep the spoken “default” on 17, and do
   not let slide 3 overclaim.

### Framing (not errors)

10. **msolve / AlgebraicSolving.jl (slide 6)** has no paper section.
    Appropriate for this workshop; do not add it to the paper solely
    to match the talk. If a proceedings version is derived from the
    talk, it needs a short related-work paragraph.

11. **Astroid as “headline result”** vs paper’s two extensions
    (projection + deflation). The astroid is a *comparison*, not a new
    theorem. The talk is allowed to lead with it. Spoken line should
    stay “same variety, do the counts agree?” not “we beat
    Bertini_real.”

12. **Zoom diagram (`\ZoomDiagram`)** is defined in `talk.tex` and
    never used. Dead code, not a content error.

13. **Pipeline step numbering on slides 11–16** (curve steps 1–6,
    surface steps 1–7) does not match paper §3.3’s six stages. It
    *does* match the figure-generation scripts. Fine for the talk;
    do not cite “Step 4” from a slide as if it were paper §3.4.

---

## Suggested slide updates (do not apply here)

| Slide | Change |
|---|---|
| 19 | “four validation surfaces”; torus is the one that fails completely |
| 20 | `537 ± 1` or a pinned live count, not a bare 538 |
| 16 | “188 → 34 on this run (paper: 31–35)” |
| 8 or new | One residual row (sphere / ellipsoid max / torus\*) or Figure 5 |
| 18→19 | Optional: restore `sample_edge` before/after from the lean deck |
| 23 | “does not use the deflated system,” not “does not retain” |
| 7 | Keep the TikZ `{f, ∂_y f}` — it is correct; fix the paper instead |

If time forces a cut, do **not** cut slides 17–18 (plain vs incidence)
or 22 (Whitney). Those are the paper’s two extensions made visible.
Cut sphere/ellipsoid visuals first; they are already missing.

---

## Mapping: paper sections → slides

| Paper | Talk coverage |
|---|---|
| §1 CAD / motivation | Slides 4–5 (compressed) |
| §2.1 homotopy | Slide 4, one italic line |
| §2.3 critical fibers | Implied by pipeline; formula not shown (talk’s TikZ is the better statement) |
| §2.4 no NID | Slide 5 |
| §3 types / config | Slide 7 caption only; Table 2 never shown |
| §3.5 deflation | Slides 22–23 (late) |
| §3.6–3.8 pipeline | Slides 11–16 (worked examples, not definitions) |
| §3.9 watertightness | Slide 16, one number |
| §4.1 sphere | Missing |
| §4.2 Taubin | Slides 13–18 |
| §4.3 torus | Slide 19 |
| §5 residuals | Slide 8, no numbers |
| §5.2 degenerate slice | Slide 15 caption (“perturbed away from a degenerate naive midpoint”) |
| §5.3 tests | Slide 20, different count |
| §6 astroid | Slides 9–10 (promoted) |
| §6 comparison table | Slide 21 (compressed; paper Table 5 is richer) |
| §7 availability | Slide 24 |
