# Patch summary (2026-08-17)

Factual corrections only. Narrative, architecture sections, and timings were not rewritten. Source of the findings: `notes/MANUSCRIPT_REVIEW.md`, `notes/TALK_ALIGNMENT.md`, `notes/DATA_CONSISTENCY.md`.

## `paper/HomotopyGetsReal_paper_current.md`

| Location | What | Was | Now |
|---|---|---|---|
| §2.3 (line 50) | Plane-curve critical system (example for projection onto \(x\)) | `{f, ∂ₓf, ∂ᵧf}` | `{f, ∂ᵧf}` |
| §3 opening (line 60) | Library size | approximately 3700 lines of Julia | approximately 5800+ lines of Julia (live `src/*.jl` = 5867) |
| §4.3 torus mesh (line 212) | Triangle count | 160,086 | 160,381 (matches `paper_artifacts/data/results.json`) |
| Table 3, Taubin heart row (line 236) | Median residual | \(2.31\times 10^{-8}\) | \(3.44\times 10^{-8}\) (JSON `3.437…e-8`) |

**Left unchanged on purpose**

- §3.4 surface augmentation `{f, ∂ₓf, ∂ᵧf}` (three variables) — correct for surfaces.
- Table 3 sphere-coarse *mean* \(2.31\times 10^{-8}\) — different cell, matches JSON.
- Wall-clock timings (sphere 3.94 s / 23.93 s; Taubin 18.52 s / 24.86 s).
- Test-suite count 537/537.
- Torus vertex count 78,340.

## `paper_artifacts/figures/src/talk.tex`

| Slide | What | Was | Now |
|---|---|---|---|
| 19 (generic projection, line 369) | Validation-set size | “three validation surfaces” | “four validation surfaces” |

**Left unchanged on purpose**

- Slide 19 “the other three surfaces (sphere, ellipsoid, Taubin heart)” — correct (torus vs the other three).
- Slide 20 test count 538/538 (paper still 537; ±1 jitter; not patched).
- Slide 16 naked-edge drop 188 → 34 (in JSON range 31–35; not patched).
- `talk_lean_23slides.tex` still says “three validation surfaces”.

## Compile

From `paper_artifacts/figures/src/`: `./compile_talk.sh` exited 0.

- `paper_artifacts/figures/exports/talk.pdf` — 25 pages
- `paper_artifacts/figures/exports/talk_lean_23slides.pdf` — 23 pages

Warnings only (overfull boxes, PDF 1.7 vs 1.5). No TeX errors.

## Still open (not patched this pass)

See the three review notes. Highest remaining items for the lead editor: talk 538 vs paper 537; torus triangle count vs other observed runs if Figure 4 is retied; stale `results.json` `tests` block (162/178); paper timings vs JSON; talk future-work “does not retain” vs `F_final` retained unused.
