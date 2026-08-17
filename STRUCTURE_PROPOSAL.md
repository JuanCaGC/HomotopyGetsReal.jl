# Structure proposal (2026-08-17)

Audit of everything outside `src/`, `test/`, and the already-organized
`paper_artifacts/figures/` tree. Status: **executed the same day** — the
layout in `REPO_ORGANIZATION.md` is this proposal after the moves.

This is a **registered Julia package** (General registry, Documenter.jl CI,
`Pkg.test()`). It is **not** a Lean project and **not** a LaTeX-first
monorepo. The proposal follows computational-math hygiene (separate
manuscript, talks, notes, generated artifacts, compiled exports) *without*
breaking Julia package conventions.

## What is actually here (verified)

| Kind | Finding |
|---|---|
| Lean | **None.** Zero `.lean` files, no `lakefile`, no `lean-toolchain`. |
| Paper `.tex` | **None.** The manuscript is markdown: `HomotopyGetsReal_paper_current.md`. `paper_artifacts/figures/exports/main.pdf` is a compiled snapshot of unknown provenance, not produced by a repo `.tex`. |
| Talk `.tex` | Two Beamer decks, already in `paper_artifacts/figures/src/` after the figures reorg. |
| Package docs | Documenter.jl at `docs/` (`docs/make.jl`, `docs/src/`, CI in `.github/workflows/Documentation.yml`). |
| Generated paper/talk figures | `paper_artifacts/figures/{paper,talk,exports,archive}/` (do not reshuffle). |
| Generation scripts + numbers | `paper_artifacts/*.jl`, `results.json`, `residuals/`, mixed with logs at the `paper_artifacts/` root. |
| Working notes | `docs/DESIGN_NOTES.md` (tracked, load-bearing); `dev/scratch/` (tracked evidentiary scripts + untracked talk-prep notes). |
| Reference implementations | `prototipo_viejo_julia/` (tracked); `codigo_cplusplus/` (gitignored). |
| Academic PDFs | `referencias_academicas/` (one tracked author PDF; two gitignored manuals). |

## Current tree (in-scope, before this reorg)

```
HomotopyGetsReal/                  Julia package root — keep
├── src/  test/  docs/  examples/  UNTOUCHED (Pkg + Documenter + CI)
├── paper_artifacts/               gitignored; messy at its own root
│   ├── *.jl                       generation scripts mixed with logs
│   ├── *.log
│   ├── results.json, residuals/
│   ├── astroid_bertini_real_sampled.{pdf,png}   duplicate of figures/paper/curves/
│   ├── VISUAL_ASSETS.md
│   └── figures/                   already organized; leave
├── dev/scratch/                   evidentiary record (cited by path)
│   ├── scratch_*.jl, capability_survey/, oscar_investigation/
│   ├── TALK_PREP.md, REPO_TECHNICAL_SUMMARY.md   talk-prep, not scratch scripts
├── prototipo_viejo_julia/         UNTOUCHED
├── referencias_academicas/        UNTOUCHED
├── codigo_cplusplus/              gitignored; UNTOUCHED
└── figures_organization.md        keep at root (figures-specific map)
```

## Rejected layouts (and why)

### Root-level `lean/`

There is nothing to put in it. An empty Lean tree would misrepresent the
project to coworkers.

### Root-level `talks/` that relocates `talk.tex`

The figures reorg just set `\figpaper` / `\figtalk` relative to
`paper_artifacts/figures/src/`. Moving the decks again would re-break
those paths for no gain. Talk sources stay where they compile.

### Root-level `artifacts/` rename of `paper_artifacts/`

Dozens of citations used `paper_artifacts/results.json` (now
`paper_artifacts/data/results.json`), plus `paper_artifacts/VISUAL_ASSETS.md`
and many `docs/DESIGN_NOTES.md` paths under that prefix. Renaming the
directory is a citation bomb, not a cleanup.

### Moving `docs/`, `examples/`, `dev/scratch/{capability_survey,oscar_investigation,scratch_*.jl}`

- `docs/` is Documenter's required location (`docs/make.jl`, Documentation.yml).
- `examples/oscar_integration.jl` is the public example; it cites
  `dev/scratch/oscar_investigation/` by path.
- Scratch scripts are the evidentiary record `test/test_isosingular_deflation.jl`
  and `docs/DESIGN_NOTES.md` cite by filename. Moving them rewrites history
  of a lab notebook.

### Flattening `prototipo_viejo_julia/` into `src/`

Forbidden by `.cursorrules`. Leave it as the archived prototype.

## Proposed layout (what we will do)

Match the spirit of `paper/` / `notes/` / `exports/` / `artifacts/`
*inside the constraints above*:

```
HomotopyGetsReal/
├── paper/                         NEW — manuscript only
│   └── HomotopyGetsReal_paper_current.md
├── notes/                         NEW — talk/exam prep, not scratch validation
│   ├── TALK_PREP.md
│   └── REPO_TECHNICAL_SUMMARY.md
├── paper_artifacts/               stays the artifacts root (gitignored)
│   ├── Project.toml, Manifest.toml, setup.jl, common.jl, generate_all.jl
│   ├── scripts/                   NEW — all example/regen .jl except common/setup/generate_all
│   ├── data/                      NEW — results.json, residuals/, figure_regen_report.json
│   ├── logs/                      NEW — *.log
│   ├── figures/                   unchanged
│   └── VISUAL_ASSETS.md
├── STRUCTURE_PROPOSAL.md          this file
├── REPO_ORGANIZATION.md           final map + compile commands
└── figures_organization.md        figures-specific map (already exists)
```

No `lean/`. No second `exports/` at repo root (`paper_artifacts/figures/exports/`
already holds compiled decks). No `talks/` directory; talk sources remain
`paper_artifacts/figures/src/`.

## Rationale per folder

| Folder | Role |
|---|---|
| `paper/` | Manuscript source, separated from figure binaries and from Documenter (`docs/`). Previously buried in `figures/src/` next to Beamer aux concerns. |
| `notes/` | Human working notes that are *not* the Documenter site and *not* the scratch-phase validation harness. |
| `paper_artifacts/scripts/` | Regenerators. A Julia mini-project should not mix 7 run logs with the scripts that produce figures. |
| `paper_artifacts/data/` | Measured numbers (`results.json`) and residual arrays. Cited as the paper's residual source. |
| `paper_artifacts/logs/` | Regenerable stdout captures. Not inputs. |
| `paper_artifacts/figures/` | Visual assets; already split into `paper/`, `talk/`, `src/`, `exports/`, `archive/`. |
| `docs/` | Package documentation (Documenter) + `DESIGN_NOTES.md`. |
| `dev/scratch/` | Dated investigation record. Stay put because other files cite those paths. |

## Path updates required when executing

- Every `paper_artifacts/*_example.jl` `include` of `common.jl` → `joinpath(@__DIR__, "..", "common.jl")`.
- `common.jl`: `RESULTS_PATH` / `RESIDUALS_DIR` under `data/`.
- `generate_all.jl`: include from `scripts/`; print `data/results.json`.
- `results.json`: rewrite stored absolute paths (`residuals/` → `data/residuals/`; stale `figures/*.pdf` → `figures/paper/...` and talk-hook → `figures/talk/...`).
- `CLAUDE.md`, `VISUAL_ASSETS.md`, `figures_organization.md`, `TALK_PREP.md`, `README.md` Layout table.
- Delete byte-identical duplicate `paper_artifacts/astroid_bertini_real_sampled.{pdf,png}` (canonical copy is `figures/paper/curves/`).

## Verify

- Recompile `paper_artifacts/figures/src/talk.tex` and `talk_lean_23slides.tex` (Tectonic). There is no paper `.tex` to compile; do not invent one.
- Confirm `results.json` residual paths exist on disk after the `data/` move.
- Do **not** run `Pkg.test()` or `docs/make.jl` for this change: no `src/` or Documenter source moved.
