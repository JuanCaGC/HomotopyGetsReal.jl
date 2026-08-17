# Overleaf Git sync

Local working copies of the Albatross talk and the JSAG paper live in two
Overleaf projects. This clone talks to them through two Git remotes. The
authentication token is stored only in this machine's `.git/config` (untracked)
and is never written to Markdown, TeX, or Julia sources.

## Remotes

| Remote | Overleaf project | Project ID | Branch |
|---|---|---|---|
| `overleaf-talk` | Talk deck | `6a6fc9f551a6522b216b3962` | `main` |
| `overleaf-paper` | Paper | `6a678fcd283cb58b0cb332c2` | `main` |

Clone URL form (token omitted on purpose):

```
https://git.overleaf.com/<project-id>
```

Username is `git`. The password is the Overleaf Git token, already saved in
`.git/config` for this clone. To rotate the token later:

```bash
git remote set-url overleaf-talk  https://git:<NEW_TOKEN>@git.overleaf.com/6a6fc9f551a6522b216b3962
git remote set-url overleaf-paper https://git:<NEW_TOKEN>@git.overleaf.com/6a678fcd283cb58b0cb332c2
```

Do not put the token in this file, in the sync script, or in any commit.

The two Overleaf histories have **no commits in common** with `origin/main`.
Never merge `overleaf-talk/main` or `overleaf-paper/main` into `main`. The
helper script uses dedicated worktrees under `paper_artifacts/overleaf/`
(gitignored) so the Julia package tree stays untouched.

## Helper script

```bash
paper_artifacts/scripts/sync_overleaf.sh fetch  [talk|paper|all]
paper_artifacts/scripts/sync_overleaf.sh pull   [talk|paper|all]
paper_artifacts/scripts/sync_overleaf.sh status [talk|paper|all]
paper_artifacts/scripts/sync_overleaf.sh sync   [talk|paper|all] [--commit]
```

| Command | What it does |
|---|---|
| `fetch` | `git fetch` the named remote(s). Auth check. |
| `pull` | Fetch, then fast-forward the matching worktree. Aborts on divergence. |
| `status` | Remote tip plus worktree dirty state. |
| `sync` | Pull, copy local assets into the worktree, stop unless `--commit`. |
| `sync … --commit` | Same, then commit in the worktree and push to Overleaf `main`. |

`push` is an alias for `sync`. Without `--commit` the copy is left uncommitted
so it can be reviewed in `paper_artifacts/overleaf/{talk,paper}`.

The script never prints remote URLs (they contain the token).

## What gets copied

**Talk** (`overleaf-talk`):

- `paper_artifacts/figures/src/talk.tex` → Overleaf `main.tex`, with
  `\figpaper` / `\figtalk` rewritten from `../paper` / `../talk` to
  `./figures` so the deck compiles at the Overleaf project root.
- `paper_artifacts/figures/src/talk_lean_23slides.tex` → `figures/talk_lean_23slides.tex` (same rewrite).
- `paper_artifacts/figures/paper/**` → `figures/**` (keeps `curves/`, `surfaces/`, `pipeline_steps_*`).
- `paper_artifacts/figures/talk/**` → `figures/**`.

**Paper** (`overleaf-paper`):

- Overleaf already holds the split TeX (`main.tex`, `section*.tex`, bibliography).
- Local `paper/HomotopyGetsReal_paper_current.md` is **not** copied over those files.
- `paper_artifacts/figures/paper/**` → Overleaf `figures/**`.

**Excluded:** any file larger than 50 MB (Overleaf's per-file limit). That
currently skips `paper_artifacts/figures/paper/surfaces/torus.pdf` (~81 MB).
The talk uses `figures/talk/torus_slide_raster.png` instead.

## Conflict rule

`pull` and `sync` use `--ff-only`. If Overleaf and the worktree have diverged,
the script stops and leaves the worktree dirty. Resolve inside
`paper_artifacts/overleaf/<talk|paper>`, then re-run. Do not force-push unless
you intend to overwrite the Overleaf project.

## Connection check (2026-08-17)

Both remotes accepted the token. `git fetch overleaf-talk` and
`git fetch overleaf-paper` succeeded. Each remote tracks `main`
(`overleaf-talk/main`, `overleaf-paper/main`). Histories are unrelated to
this package's `main`, as expected.
