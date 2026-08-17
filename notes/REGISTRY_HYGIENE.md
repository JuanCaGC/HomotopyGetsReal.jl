# Registry hygiene audit (2026-08-17)

Empirical pass against the live working tree, `git ls-files`, and
`.gitignore`. Nothing in `src/`, `test/`, or `docs/` was deleted or moved.
The numbered sections below are the pre-cleanup audit. The cleanup itself
was executed the same day; see **Cleanup Executed**.

Package is already in General as HomotopyGetsReal v0.2.1
(`Project.toml` uuid `4051f461-4a5e-492e-9237-10b01d865883`).
`Pkg.add` fetches the **tagged tree**, not this machine’s gitignored
directories. `git clone` of GitHub fetches **history**, including blobs
that are now ignored.

---

## Cleanup Executed (2026-08-17)

`git rm --cached` only — local files were not deleted. `.gitattributes`
`export-ignore` rules were added. Core package paths were not edited.

### Index metrics

| | Files | Size |
|---|---|---|
| Tracked before | 249 | 15.93 MB |
| Tracked after (index) | **111** | **2.84 MB** |
| `git archive` with `export-ignore` | — | **700 KB** uncompressed, **~211 KB** gzip |

### Untracked from the index (still on disk)

- 78 PNGs under `dev/renders/` and `dev/scratch/`
- 44 logs under `dev/scratch/`
- 28 JSON dumps under `dev/scratch/`
- root `Manifest.toml` (libraries pin `Project.toml` only)

`dev/` still tracks **41 `.jl` scripts and 6 markdown notes**. Scratch
renders/dumps stay on disk and are now gitignored (`dev/**/*.png`,
`dev/**/*.json`, `*.log`).

### Core code

`src/`, `test/`, and `docs/` were unaffected: no staged edits, no
deletes, no moves. Unrelated unstaged leftovers from earlier work
(`docs/DESIGN_NOTES.md`, `CLAUDE.md`, `README.md`,
`examples/oscar_integration.jl`) were left alone.

---

## 1. Footprint

### On-disk working tree (this checkout)

| Path | Disk | Git status |
|---|---|---|
| Whole tree | **573 MB** | mix |
| `paper_artifacts/` | **228 MB** | ignored (`paper_artifacts/`) |
| `paper_artifacts/figures/paper/` | 167 MB | ignored (torus.pdf 81 MB, ellipsoid/sphere ~20 MB each) |
| `paper_artifacts/figures/archive/` | 38 MB | ignored |
| `paper_artifacts/figures/exports/` | 17 MB | ignored (`talk.pdf` 7.0 MB, lean deck 4.6 MB, `main.pdf` 2.4 MB) |
| `paper_artifacts/data/` | 2.5 MB | ignored (`results.json`, residual dumps; torus residuals 1.6 MB) |
| `paper_artifacts/logs/` | 192 KB | ignored |
| `dev/` | **167 MB** | mostly tracked; two PDFs ignored (14 MB + 138 MB) |
| `.git/` | **88 MB** | history |
| `codigo_cplusplus/` | 62 MB | ignored |
| `referencias_academicas/` | 11 MB | 1 PDF tracked; two manuals ignored |
| `docs/` | 2.7 MB | source tracked; `build/`, `BERTINIREAL_AUDIT.md`, `ORCHESTRATOR_BRIEFING.md` ignored |
| `src/` | 296 KB | tracked |
| `test/` | 220 KB | tracked (`test/output/` ignored) |
| `notes/` | 180 KB | **untracked** |
| `paper/` | 68 KB | **untracked** |
| `prototipo_viejo_julia/` | 36 KB | tracked |
| `examples/` | 8 KB | tracked |

`paper_artifacts/` by extension (all untracked): 43 PDF, 38 PNG, 9 log,
8 JSON, 2 tex, plus Beamer `*.aux` / `*.nav` / `*.toc` from
`compile_talk.sh`. Those intermediates are why the `.gitignore` update
names LaTeX suffixes even though the parent directory was already ignored.

### What Git actually tracks (clone checkout)

| | |
|---|---|
| Tracked files | **249** |
| Tracked bytes | **15.93 MB** |
| Of which `dev/` | **13.96 MB / 197 files** |
| Of which `referencias_academicas/Homotopy gets real.pdf` | **1.62 MB** |
| Everything else | **~0.35 MB** |

Tracked by extension: **78 PNG, 73 jl, 44 log, 28 json**, 15 md, 4 toml,
3 yml, **1 pdf**. Every PNG, log, and JSON is under `dev/`.

Largest tracked files (all already in the index):

- `referencias_academicas/Homotopy gets real.pdf` — 1.6 MB
- `dev/scratch/renders_pipeline_steps/*.png` — up to 762 KB
- `dev/scratch/renders_paper_final/*.png` — 200–486 KB
- `dev/scratch/oscar_investigation/setup_and_time.log` — 102 KB

No `paper_artifacts/` file is tracked. That part of isolation is already
correct.

### Clone size vs package-only size

| What a user gets | Size | Notes |
|---|---|---|
| `git clone` (`.git` + tracked checkout) | **~104 MB** | `.git` 88 MB + 16 MB files |
| Tagged tree `Pkg.add` would archive today | **~16 MB** | full 249-file HEAD, including `dev/` |
| Package-only (src, test, docs, examples, Project.toml, LICENSE, README, CITATION, CHANGELOG, .github, CLAUDE.md, prototipo) | **~0.71 MB / 45 files** | `git ls-files` minus `dev/` and the academic PDF |
| Strict library core (src, test, Project.toml, LICENSE, README) | **436 KB** | what General actually needs to run tests |

`.git` is 88 MB because history still holds blobs that `.gitignore` no
longer admits, notably:

- `codigo_cplusplus/python/docs/tutorials/anaglypy_pictures/schneeflocke_raw_smooth.gif` — **22 MB**
- other anaglyph GIFs / collages — 1–11 MB each
- `referencias_academicas/bertini_real_manual.pdf` — 8.1 MB (now ignored, still in history)
- `referencias_academicas/BertiniUsersManual.pdf` — 1.2 MB (same)

Those do **not** ship in a new `Pkg.add` of v0.2.1, but they do inflate
every full clone and every CI `actions/checkout` of this repo.

Untracked-and-not-ignored markdown (`paper/`, `notes/`,
`REPO_ORGANIZATION.md`, `STRUCTURE_PROPOSAL.md`,
`figures_organization.md`) is **~250 KB**. Safe to add later; they do not
move clone size.

---

## 2. Julia General Registry readiness

### Package core — clean enough

| File | Status |
|---|---|
| `Project.toml` | name, uuid, MIT, v0.2.1, authors, `[compat]` including `julia = "1.12"`, `[extras]`/`[targets]` for Test |
| `src/HomotopyGetsReal.jl` | module name matches package name; 11 `src/*.jl` files |
| `test/runtests.jl` | standard; slow Taubin gated on `HOMOTOPYGETSREAL_RUN_SLOW_TESTS` |
| `docs/` | Documenter (`docs/Project.toml`, `docs/src/`, `docs/make.jl`); `docs/build/` ignored |
| `LICENSE` | MIT, 2026, Juan Camilo Gonzalez |
| `README.md` | install via `Pkg.add`, Julia 1.12+, test instructions |
| `CITATION.cff` | version matches Project.toml 0.2.1 |
| `.github/workflows` | CI (fast on PR, slow on cron), Documentation, TagBot |

`src/`, `test/`, and `docs/` source are independent of `paper_artifacts/`
and of `codigo_cplusplus/`. Tests write under `test/output/` (ignored).
Documenter CI develops the package from `pwd()`; it does not need the
paper tree.

### Issues that are not blockers for a already-registered v0.2.1, but are not best practice

1. **Root `Manifest.toml` is tracked (75 KB).** Pkg’s rule for
   *libraries* is: commit `Project.toml`, do **not** commit
   `Manifest.toml`. Apps/environments commit both. HomotopyGetsReal is a
   library. `docs/Manifest.toml` is on disk and **not** tracked (the
   comment in `.gitignore` that says it is tracked is stale). Satellite
   `dev/scratch/*/Manifest.toml` are already ignored.

2. **`dev/` is inside the registered tree.** 78 PNGs + 44 logs + 28 JSON
   files are what a General tarball of HEAD would include. They are not
   imported by `src/`. They are an evidentiary lab notebook, cited from
   `docs/DESIGN_NOTES.md` and `examples/oscar_integration.jl`. Fine in
   *this git repo*; they should not ride along in the next registered
   tag.

3. **No `.gitattributes` `export-ignore`.** GitHub’s tag tarball (what
   Registrator/Pkg uses) honors `export-ignore`. That is the least
   disruptive way to keep one git remote and a small package archive.

4. **GLMakie is a hard `[deps]` entry.** Every `using HomotopyGetsReal`
   pulls an OpenGL stack (CI already needs xvfb). Registry-legal;
   optional weakdep / extension would be a later design change, not a
   hygiene fix.

5. **`paper/` and `notes/` are untracked.** They do not affect the
   current registered tree. Adding them as markdown would not inflate
   Pkg.add in any meaningful way.

### Isolation of non-code directories

| Directory | Isolated from Pkg.add today? | Verdict |
|---|---|---|
| `paper_artifacts/` | yes (gitignored, 228 MB) | keep ignored; do not track |
| `codigo_cplusplus/` | yes (gitignored); **no** for clone history | keep ignored; history rewrite is optional and costly |
| `paper/`, `notes/` | N/A (untracked, tiny) | safe to track as markdown; do not put PDFs here |
| `dev/` | **no** — 14 MB in the tagged tree | untrack binaries or `export-ignore` before the next tag |
| `prototipo_viejo_julia/` | in the tree, 24 KB | harmless; still not package code |
| `referencias_academicas/Homotopy gets real.pdf` | in the tree, 1.6 MB | author-owned; keep or move to the paper repo |

---

## 3. `.gitignore` update (this pass)

Added, without removing existing rules:

- `paper_artifacts/data/`, `paper_artifacts/logs/`,
  `paper_artifacts/figures/exports/` (explicit, redundant with
  `paper_artifacts/`, documents the heavy paths)
- `*.pdf` with `!referencias_academicas/Homotopy gets real.pdf`
- LaTeX/Tectonic/Beamer suffixes: `*.aux *.bbl *.blg *.fdb_latexmk
  *.fls *.nav *.out *.snm *.synctex.gz *.toc *.vrb *.xdv *.dvi`
- `*.log` (the `!dev/scratch/**/*.log` exception was removed once those
  logs were untracked; see Cleanup Executed)

**Done in Cleanup Executed:** `git rm --cached` of the 78 PNGs / 44 logs
/ 28 JSON and root `Manifest.toml`. `.gitattributes` `export-ignore` is
in place. Section 4 below is the original recommendation list, kept as
the audit record.

---

## 4. Recommended for removal from Git *tracking*

Do not delete the files on disk. Suggested `git rm --cached` groups,
largest first:

| Group | Tracked now | Why |
|---|---|---|
| `dev/scratch/**/*.png`, `dev/renders/**/*.png` | 78 files, ~10+ MB | regenerable renders; DESIGN_NOTES already cites the findings |
| `dev/scratch/**/*.log` | 44 files | solver transcripts, not package API |
| `dev/scratch/**/*.json` | 28 files | capability-survey dumps; paper numbers live in ignored `paper_artifacts/data/` |
| root `Manifest.toml` | 75 KB | library lockfile should not be pinned for General users |

Keep tracked: `src/`, `test/`, `docs/src/`, `docs/make.jl`,
`docs/Project.toml`, `docs/DESIGN_NOTES.md`, `docs/README.md`,
`examples/`, `Project.toml`, `LICENSE`, `README.md`, `CHANGELOG.md`,
`CITATION.cff`, `.github/`, `.claude/agents/`.

Optional keep: `prototipo_viejo_julia/` (tiny),
`referencias_academicas/Homotopy gets real.pdf` (author-owned; 1.6 MB),
`dev/scratch/**/*.jl` scripts (the actual evidentiary *code*, not the
renders).

History shrinkage (`git filter-repo` of `codigo_cplusplus/` GIFs) would
cut clone size from ~104 MB toward ~20 MB. That rewrites SHAs, fights
TagBot/CI, and is **not** required for the next registry bump. Do it only
on an explicit decision.

---

## 5. Unified repo vs split

**Recommendation: keep one git remote (`HomotopyGetsReal.jl`), do not
split the registered package out.**

Reasons:

1. v0.2.1 is already registered against this GitHub URL. A split forces
   a repo-move dance in General for no user-facing gain: `Pkg.add` never
   saw `paper_artifacts/` and would not see `paper/` / `notes/` either
   unless they are committed.
2. The 228 MB paper/talk binary tree is **already** isolated by
   `.gitignore`. Clone size is a *history* problem (`codigo_cplusplus`
   GIFs) plus a *tagged-tree* problem (`dev/` PNGs), not a “paper lives
   in the same folder” problem.
3. `paper/` and `notes/` are markdown (currently 68 KB + 180 KB). They
   belong next to `docs/DESIGN_NOTES.md` for the author; they are not a
   second package.

**Do this instead of a split, before the next tag:**

1. `git rm --cached` the `dev/` PNG/log/JSON groups (keep the `.jl`
   scripts).
2. Add `.gitattributes` `export-ignore` for `dev/`,
   `prototipo_viejo_julia/`, and `referencias_academicas/` so even a
   forgotten binary cannot enter the Registrator tarball.
3. Stop committing root `Manifest.toml`.
4. Leave `paper_artifacts/` gitignored. If JSAG/Zenodo need the 81 MB
   torus PDF and `results.json`, publish them as a **release asset or
   Zenodo record**, not as a second Julia package.

**When a split *would* be worth it:** if you want a public
`homotopygetsreal-paper` that *tracks* the compiled figures and
`results.json` so a reviewer can clone and recompile the talk without
this repo’s Julia package layout. That would be a reproducibility
archive (maybe with `paper_artifacts/` and `paper/`), **in addition to**
the registered package, not a replacement. Do not register that archive
with General.

Until then: unified repo, ignored artifacts, smaller tracked tree.
