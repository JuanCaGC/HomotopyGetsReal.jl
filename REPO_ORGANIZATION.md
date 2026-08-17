# Repository organization

Canonical map for coworkers. Figures detail lives in `figures_organization.md`.
The decision record for this layout is `STRUCTURE_PROPOSAL.md`.

This is a **Julia package** (HomotopyGetsReal.jl, v0.2.1). Package-required
paths (`src/`, `test/`, `docs/`, `Project.toml`) are the Documenter/Pkg
layout and are not to be moved. There is **no Lean code** and **no paper
`.tex`**.

## Layout

```
HomotopyGetsReal/
├── src/                           package source — do not reorganize
├── test/                          Test.jl suite — do not reorganize
├── docs/                          Documenter.jl site + DESIGN_NOTES.md
├── examples/                      public examples (oscar_integration.jl)
├── paper/                         manuscript (markdown)
│   └── HomotopyGetsReal_paper_current.md
├── notes/                         talk/exam prep (not scratch validation)
│   ├── TALK_PREP.md
│   └── REPO_TECHNICAL_SUMMARY.md
├── paper_artifacts/               gitignored generators + figures + numbers
│   ├── scripts/                   Julia figure regenerators
│   ├── data/                      results.json, residuals/, figure_regen_report.json
│   ├── logs/                      regenerable stdout captures
│   ├── figures/                   visual assets (see figures_organization.md)
│   │   ├── paper/                 current paper renders
│   │   ├── talk/                  talk-only renders
│   │   ├── src/                   Beamer sources + talk scripts
│   │   ├── exports/               compiled talk.pdf, main.pdf
│   │   └── archive/               dated superseded renders
│   ├── common.jl, generate_all.jl, setup.jl
│   └── VISUAL_ASSETS.md           what each render shows
├── dev/scratch/                   evidentiary investigation record (cited by path)
├── prototipo_viejo_julia/         archived Julia prototype
├── referencias_academicas/        manuals (two gitignored; one tracked)
├── codigo_cplusplus/              gitignored Bertini_real tree
├── STRUCTURE_PROPOSAL.md
├── REPO_ORGANIZATION.md           this file
└── figures_organization.md        figures-only map
```

### What each top-level extra folder is for

| Path | Put here | Do not put here |
|---|---|---|
| `paper/` | The JSAG/manuscript markdown | Figure binaries, Beamer, Documenter pages |
| `notes/` | Spoken-prep and exam-style notes | `docs/DESIGN_NOTES.md` (that stays in `docs/`); scratch `scratch_*.jl` |
| `paper_artifacts/` | Regenerable paper/talk pipeline | Package source |
| `dev/scratch/` | Dated validation scripts and survey data cited from tests/DESIGN_NOTES | New talk prep (that goes in `notes/`) |
| `docs/` | Documenter + the design-notes lab notebook | Talk decks |

## Compile commands

### Talk decks (Beamer)

There are two `.tex` files, both under `paper_artifacts/figures/src/`.
Paths `\figpaper` / `\figtalk` are relative to that directory.

```bash
cd paper_artifacts/figures/src
./compile_talk.sh
```

Writes:

- `paper_artifacts/figures/exports/talk.pdf` — current 25-slide deck
- `paper_artifacts/figures/exports/talk_lean_23slides.pdf` — older 23-slide cut

Requires `tectonic` or `pdflatex` with Beamer + metropolis. From `src/`, Tectonic:

```bash
cd paper_artifacts/figures/src
tectonic -p --keep-logs --keep-intermediates -o ../exports talk.tex
tectonic -p --keep-logs --keep-intermediates -o ../exports talk_lean_23slides.tex
```

### Paper

There is **no** `paper/*.tex` in this repository. The manuscript is
`paper/HomotopyGetsReal_paper_current.md`. A PDF snapshot exists at
`paper_artifacts/figures/exports/main.pdf`; it is not produced by a
checked-in build script.

### Package docs (Documenter.jl)

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Needs a desktop session so GLMakie can load. Output: `docs/build/index.html`.
On CI, wrap with `xvfb-run` (see `.github/workflows/Documentation.yml`).

### Tests

```bash
julia --project -e 'using Pkg; Pkg.test()'                          # fast
HOMOTOPYGETSREAL_RUN_SLOW_TESTS=1 julia --project -e 'using Pkg; Pkg.test()'
```

### Regenerating paper/talk figures

```bash
julia --project=paper_artifacts paper_artifacts/setup.jl   # first time only
julia --project=paper_artifacts paper_artifacts/generate_all.jl
# or individual scripts:
julia --project=paper_artifacts paper_artifacts/scripts/curve_examples.jl
julia --project=paper_artifacts paper_artifacts/scripts/sphere_example.jl
julia --project=paper_artifacts paper_artifacts/scripts/ellipsoid_example.jl
julia --project=paper_artifacts paper_artifacts/scripts/taubin_example.jl
julia --project=paper_artifacts paper_artifacts/scripts/taubin_singular_structure_example.jl
julia --project=paper_artifacts paper_artifacts/scripts/torus_example.jl
julia --project=paper_artifacts paper_artifacts/scripts/pipeline_steps_curve_example.jl
julia --project=paper_artifacts paper_artifacts/scripts/pipeline_steps_surface_example.jl
julia --project=paper_artifacts paper_artifacts/scripts/residual_histograms.jl
```

Numbers land in `paper_artifacts/data/results.json`. Renders land under
`paper_artifacts/figures/paper/` and `paper_artifacts/figures/talk/`.

## Gitignore vs tracked

`paper_artifacts/` is gitignored (regenerable). `paper/` and `notes/` are
*not* gitignored: they hold authoring sources that used to be buried
inside the gitignored tree. Add them when you want them in git.

`docs/DESIGN_NOTES.md` is tracked. `docs/BERTINIREAL_AUDIT.md` and
`docs/ORCHESTRATOR_BRIEFING.md` are gitignored local files; leave them
at those paths (`CLAUDE.md` says so).

`dev/scratch/capability_survey/` and `dev/scratch/oscar_investigation/`
are mostly tracked evidentiary record. Do not relocate; tests and
`examples/oscar_integration.jl` cite them by path.

## Known stale pointer (not silently "fixed")

`dev/scratch/README.md` refers to `test/ASSERTION_AUDIT.md`. That file
is `dev/scratch/ASSERTION_AUDIT.md`. Left as written; the scratch README
is itself historical.
