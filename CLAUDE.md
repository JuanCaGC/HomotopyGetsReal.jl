# CLAUDE.md

This file supersedes `docs/ORCHESTRATOR_BRIEFING.md`. That file is
gitignored and stays on disk, but its header should be marked historical —
this file is now the canonical entry point for any agent/session working in
this repo.

## What this is

HomotopyGetsReal.jl is a Julia reimplementation of "Homotopy gets real" —
numerical algebraic geometry for decomposing real algebraic curves and
surfaces into vertices/edges (1D) or vertices/edges/faces/mesh (2D surfaces
in 3-space), built on
[HomotopyContinuation.jl](https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl).

- **Version**: `0.2.1`, registered in Julia's General registry (confirmed
  live against the registry as of this writing — `Versions.toml` for this
  package in `JuliaRegistries/General` lists `0.1.0`, `0.2.0`, and `0.2.1`).
- **License**: MIT, copyright Juan Camilo Gonzalez, 2026.
- **Author**: Juan Camilo Gonzalez, sole author. `git shortlog -sne --all`
  shows only `JuanCaGC` (two personal emails) plus the `Documenter.jl` CI
  bot — no analytical co-author on this codebase.

## Correctness standards

Four mechanisms are actually in use to establish correctness. This is not
a proposal — all four have live, checkable evidence in this repo (or, for
the third, in the adjacent repo it deliberately lives in).

1. **Pointwise residual checks against the defining polynomial.** Every
   decomposition claim can be checked by plugging mesh/vertex coordinates
   back into the original `f` and asserting the residual is near zero.
   Example: `test/test_surfacedecomposition.jl:263` computes
   `|x²+4y²+9z²-1|` at every welded ellipsoid mesh vertex and asserts
   `all(<=(1e-4), ell_residuals)`. The same pattern (`mean`/`p90`/`p99`
   residual stats per fixture) is computed for the paper figures in
   `paper_artifacts/data/results.json` (e.g. torus: mean `2.2e-6`, p99 `1.0e-5`).

2. **Cross-validation against Hauenstein–Wampler's published isosingular
   deflation sequences.** `src/Solver.jl`'s deflation code
   (`estimate_corank`, `deflate_once`, `verify_isosingular_dimension`) is
   built directly against Hauenstein–Wampler's `D_det` construction,
   Definition 5.18, and Algorithm 6.3. `test/test_isosingular_deflation.jl`
   reproduces their published Whitney-umbrella deflation sequences exactly
   on the bare, unsliced equation `x²-y²z`: `[3,2,0]` at the tip `(0,0,0)`
   (isolated point, resolves in two rounds) and `[3,1,1,...]` at the handle
   `(0,0,1)` (genuine plateau on the 1-dim singular locus) — raw asserted
   values at `test/test_isosingular_deflation.jl:82` and `:92`. Matches the
   paper's §3.5 pairing exactly; see `docs/DESIGN_NOTES.md` "Stage 2" for
   the write-up.

3. **Coordinate-level comparison against Bertini_real's actual output (the
   astroid case).** This comparison itself lives *outside* this repo, in
   the Correctness Auditor's external Cursor workspace (which has the real
   Bertini_real/Bertini1 install; gitignored/local by design, not part of
   this repo) — see `Astroid_edge_endpoint_crosscheck.md` there. It's
   genuine coordinate data, not a description: on `(x²+y²-1)³+27x²y²=0`, Bertini_real's real
   6-edge decomposition and HGR's 4-edge decomposition were reconciled
   exactly — the 2 extra BertiniReal edges come from its random-projection
   subdividing 2 of HGR's 4 cusp-to-cusp arcs at two smooth projection-critical
   points `S±`, confirmed on-curve (`|f(S±)| ~ 3.2e-10`, distance to the
   parametric astroid `~1.2e-5`). This repo should keep linking to that
   comparison, not re-deriving or duplicating it.

4. **The `Test.jl` suite as a topology-invariant encoding.** Assertions
   encode expected topology, not just "it ran": vertex-type counts and
   coordinates (`test/test_topology.jl:33-36`, e.g. a node curve asserted
   to have exactly 1 `Critical`, 1 `Singular`, 2 `Boundary` vertices at
   specific coordinates), critical-value tolerances (`test/test_surfacedecomposition.jl:65`,
   sphere critical z-slices asserted `isapprox(sort(z_crits), [-1.0, 1.0];
   atol=1e-6)`), and isosingular resolution across curve fixtures including
   the astroid (`test/test_isosingular_deflation.jl:127-131`, every flagged
   `Singular` vertex's deflation `verdict == Resolved`).

   **Test suite: 538/538** (fast subset: 478/478 in ~3.5 min via default
   `Pkg.test()`; full suite via `HOMOTOPYGETSREAL_RUN_SLOW_TESTS=1`:
   538/538 in ~30-34 min — `README.md`'s Testing section is kept in sync
   with this figure). The baseline moved from 537 to 538 permanently at
   commit `a53b6c4` (`plot_surface_decomposition` empty-mesh fix) — a
   fresh worktree-isolated rerun of the fast subset gave 477/477 (×3) at
   `a53b6c4~1` and 478/478 (×3) at `a53b6c4` itself, so the shift is
   reproducible and commit-aligned, not per-process RNG jitter. The exact
   code path from that fix to the count is not identified (no test file
   in `test/` exercises the empty-mesh branch it added); the full-suite
   538 figure carries this same +1 by inference and was not independently
   rerun this pass (~30-34 min). Separately, on top of that baseline, a ±1
   count variance (537 vs 538) is still expected, not a bug: the astroid
   fixture's isosingular-deflation test loop fires a variable number of
   assertions depending on live solver
   output (see paper §5.3, and `docs/DESIGN_NOTES.md`'s 2026-08-04 entry
   for the full trace). The full suite has occasionally errored outright
   (1 of 3 live attempts this session); root cause unconfirmed — check the
   log in full if it recurs rather than just rerunning. Isosingular
   deflation testing dominates full-suite runtime; a concrete but
   unconfirmed optimization lead is logged in `docs/DESIGN_NOTES.md`'s
   backlog.

## Ground rules

- **Direct empirical verification is ground truth.** A `grep`, a live test
  run, a `results.json` read — that's what counts. A claim written in a
  status document (`ORCHESTRATOR_BRIEFING.md`, a paper draft, a Slide deck,
  even this file) is not itself ground truth until independently
  re-checked against the live repo.
- **Flag genuine inconsistencies directly** — don't assume they're fine,
  don't paper over them, don't silently pick the version that sounds more
  recent or more authoritative.
- **Never resolve a flagged inconsistency unilaterally.** Report it and let
  Juan or the next role in the pipeline (Auditor, Arbitrator) decide.
- **Push back on underspecified asks** instead of quietly picking a
  default and running with it.

## Roles

- **Developer** (this environment): full read/write/execute access to this
  repo. Git discipline: never `git add -A`; always stage specific files by
  name. Propose before anything hard-to-reverse — registry-trigger tag
  pushes, force-pushes, anything that publishes externally.
- **Self-review pass**: a worktree-isolated pass Developer spawns before
  reporting any change as done. Isolation: worktree. Read-only plus test
  execution — runs `Pkg.test()` and `docs/make.jl` in that isolated
  worktree before Developer reports a task complete. Operational note:
  subagents cannot use their own `Write` tool to persist files
  ("Subagents should return findings as text, not files") — a subagent
  must return findings as text in its final message; the
  parent/orchestrating session is responsible for writing anything to
  disk. Also: a worktree-isolated subagent starts from the last commit,
  not the current working tree — it does not see uncommitted changes,
  and its own git operations are restricted to its own worktree (even
  read-only ones against the shared checkout are refused). To review an
  uncommitted change, copy the modified file(s) into the subagent's
  worktree directly via the filesystem (not `git add`/commit) before
  invoking it.
- **Correctness Auditor**: lives in Cursor, external to this repo. Has the
  actual installed Bertini_real/Bertini1 source and cross-checks new HGR
  claims against its real output (see mechanism 3 above). This repo links
  to that work (an external, gitignored/local-by-design workspace); it
  does not duplicate it — investigation and BertiniReal-install artifacts
  stay out of this repo's git history by design.
- **Arbitrator**: a separate Claude.ai Project chat with persistent memory,
  for Juan. Synthesizes across Developer/Auditor output, verifies directly
  rather than trusting descriptions, never edits code or paper content
  itself.

## Superseded / external docs — don't duplicate

- `docs/ORCHESTRATOR_BRIEFING.md` — superseded by this file. Gitignored;
  keep it on disk but mark its header historical.
- `bertini_migration/HGR_ORCHESTRATOR_BRIEF.md` (external — the
  Correctness Auditor's own local Cursor workspace, gitignored, not part
  of this repo) — stays external by design. Link to it; don't pull its
  content or the BertiniReal-install artifacts it documents into this
  repo.
