# Manuscript review

**Sources read, not edited:** `paper/HomotopyGetsReal_paper_current.md`,
`notes/REPO_TECHNICAL_SUMMARY.md`, plus live `src/` line counts and
`test/runtests.jl` only as a check on claims the two documents disagree
about.

**Rule used:** a claim in the paper is not treated as ground truth if the
technical summary, or a live file, contradicts it. Disagreements are
flagged, not silently picked.

---

## Overall

The manuscript is a competent JSAG-shaped software paper: honest about
scope (no NID, deflation diagnostic-only), has a real Bertini_real
comparison, and the residual table plus Whitney-umbrella sequences are
the right kind of evidence. It is weaker as a *mathematical* paper than
as an *implementation* paper. The load-bearing theorem in §2.3 is
sketched rather than stated, and one formula in that section is wrong
relative to the code. Against `REPO_TECHNICAL_SUMMARY.md`, the paper is
more polished and less complete: several caveats the notes treat as
talk-critical never appear.

Footer of the markdown (`paper/…:361`) says it is a transcription of
`main.tex`. There is no `main.tex` in this repository
(`REPO_ORGANIZATION.md`). Treat this file as the manuscript, not as a
shadow of a missing LaTeX tree.

---

## Mathematical rigor

### What is solid

- MidSlice-First (§3.6) is stated clearly and matches
  `src/Topology.jl:13-19` as described in the technical summary: never
  start tracking from a critical/singular vertex.
- Isosingular deflation (§3.5) is the strongest mathematical section.
  Whitney tip `[3,2,0]` and handle `[3,1,1,…]`, node `[2,0]`, cusp
  `[2,1,0]`, and the rejection of naive plateau / randomized-reduction
  tests line up with `REPO_TECHNICAL_SUMMARY.md` §1.3 and with
  `test/test_isosingular_deflation.jl` as cited there.
- The two-gate singularity test (§3.4) matches the code
  (`_classify_vertex_type`: rank deficiency **or** small trailing
  singular value).
- Scope limits in §2.4 and §6 (no NID, AABB vs sphere, diagnostic
  deflation) are stated as limitations, not as virtues.

### Load-bearing error: plane-curve critical system

§2.3 (`paper/…:50`) says the critical values of π|_V(F) come from
`{F(x), det J(x)} = 0`, then gives as the running example

> `{f, ∂ₓf, ∂ᵧf}` for a plane curve projected onto x.

That is **three equations in two unknowns**. The technical summary
already flags this (`notes/REPO_TECHNICAL_SUMMARY.md:665`, citing
`notes/TALK_PREP.md:718-737`): the code is `{f, f_y}` at
`Topology.jl:369`. Surfaces correctly use `{f, ∂ₓf, ∂ᵧf}`
(`Solver.jl` auto-augment). The talk TikZ diagram also has the code’s
version (curve lane `{f, ∂_y f}`, surface lane `{f, ∂_x f, ∂_y f}`).

This is not a wording slip. A referee who tries to reproduce the
zero-dimensional solve from §2.3 will write down an overdetermined
system. **Fix the example in §2.3 to `{f, ∂_y f}` for a plane curve
projected onto x**, and keep `{f, ∂ₓf, ∂ᵧf}` only for a surface
projected onto z.

The surrounding paragraph’s topology claim (Ehresmann / constant fiber
between consecutive critical values) is the right theorem, but it is
not named, not cited to a precise statement in Lu et al. / Besana et
al., and not given hypotheses (properness, smoothness of the regular
part, compactness via the bounding box). For JSAG this may be
acceptable; it is not a substitute for a theorem environment.

### Other mathematical gaps (paper vs notes)

| Topic | Paper | Technical summary |
|---|---|---|
| Completeness of the critical-point solve | Implicit in §2.3 | Conditional: a missed critical value silently merges slabs; cone returns empty mesh, no exception (`REPO_TECHNICAL_SUMMARY.md:615`) |
| Multiplicity ≥ 2 at the search point | Not discussed | Named as the largest correctness risk (node, astroid cusps, Whitney apex, cone) |
| `certify` / α-theory | Absent | HC.jl ships it; unused; would not certify decompositions anyway |
| `_verify_projection_ok` | Described in §6 as “vanishing of augmenting partials” | Crash guard, **not** Bertini_real’s transversality determinant |
| Path tracking always Float64 | §3.1 implies BigFloat is a full-precision mode | Tracking is Float64-only; BigFloat is post-hoc Newton polish |
| Table 2 vs live `HomotopyConfig` | 10 fields | Live struct has additional knobs: `projection_orthonormality_tol`, `z_mid_*`, `min_slab_width`, `incidence_snap_tol_ratio`, `isosingular_verify_retries`, `max_deflations`, `bbox_*` |

None of these need to become the paper’s headline. The completeness
conditional and the Float64 tracking boundary are the two that, if
omitted, oversell §3.1’s “same code path in BigFloat.”

### Line-count claim

Paper §3: “approximately 3700 lines.” Live `wc -l src/*.jl` is **5867**.
The technical summary recorded 5767. The paper figure is stale by ~2000
lines. Drop it or re-measure; do not keep 3700.

---

## Logical flow and structure

The promised map in the introduction is:

§2 math → §3 architecture → §4 examples → §5 validation → §6
Bertini_real → §7 availability.

That is the right JSAG skeleton. Execution problems:

### 1. Ellipsoid is in the evidence, not in the story

§4 opens with **three** surfaces (sphere, Taubin, torus). Table 3 and
the abstract use **four** (plus ellipsoid). §4.3 then refers to “the
sphere, ellipsoid, and Taubin heart” as if ellipsoid had been
introduced. There is no §4.x for the ellipsoid, no listing, and no
figure placeholder, even though `ellipsoid.pdf` exists and is the
asymmetric control the test suite was built around
(`REPO_TECHNICAL_SUMMARY.md` §3.2).

Either add a short §4.1b (why the ellipsoid exists: x/y/z swap cannot
hide) or stop citing it as a peer of the three narrated examples.

### 2. Curves are half the pipeline and almost none of §4

§3.6 is a full 1D section. §6’s astroid comparison is the only
coordinate-verified Bertini_real result. §4 contains **zero curve
examples**. A reader finishes “Usage examples” believing this is a
surface package. The talk inverted this (astroid first). The paper
should at least show the astroid (or the nodal cubic, which the test
suite actually asserts) before §6.

### 3. Two different “six steps”

§3.3 numbers six stages as: (1) solver core, (2) 1D decomposer, (3)
shared tracker, (4) surface sweep, (5) weld, (6) visualization.

Bertini_real / `.cursorrules` / the talk diagram number: (1) critical
points, (2) bounding object, (3) midslice, (4) connect-the-dots /
shared engine, (5) merge, (6) sample.

The paper even says “the same six steps as Bertini_real” and then
defines a different six. Visualization is not a Bertini_real pipeline
step. This is the weakest transition in the architecture section.
Pick one numbering and keep it through Figure 1, §3.4–§3.10, and §6.

### 4. Generic projection is promised in §4 and defined in §6

§4.3 uses `projection = :random` and points at “§6” for the mode.
The reader has not yet been told what the mode is. A short subsection
in §3 (or a forward-reference that actually defines it) is missing.
§6 then carries both the comparison *and* the only real statement of
the algorithm.

### 5. Validation split across §3.9, §4, and §5

Naked-edge counts live in §3.9; timings in §4; residuals in §5.1;
retry ladder in §5.2; tests in §5.3. A referee looking for “is the
Taubin heart correct?” has to assemble four sections. Consider a
single per-fixture validation block, or a pointer table at the start
of §5.

### 6. Missing sections relative to the notes’ “honest inventory”

Not required for JSAG, but currently absent and load-bearing if a
referee asks:

- No related-work paragraph on exact/real-root tools (the talk’s
  msolve slide has no paper analogue).
- No statement that the **torus is not in `test/`**
  (`REPO_TECHNICAL_SUMMARY.md` §3.3).
- No cone / empty-mesh / silent-failure example, which is the actual
  completeness caveat behind §2.3.
- Acknowledgments mention AI assistance; good. The email field is still
  a TODO (`paper/…:5`).

---

## Weak transitions (specific)

| From → to | Problem |
|---|---|
| §2.3 → §3.4 | Example system in §2.3 disagrees with the augmentation actually used in §3.4 / the code |
| §3.3 “same six steps” → §3.4–§3.10 | Different six |
| §3.1 BigFloat → §3.7 | Tracking engine is Float64; the precision story is unfinished |
| §4 open (“three surfaces”) → Table 3 (four) | Ellipsoid appears from nowhere |
| §4.3 → §6 | Generic projection used before it is defined |
| §5.1 ellipsoid max 1.12×10⁻⁴ vs test gate 10⁻⁴ | Explained in prose; still easy to misread as a failing production mesh |
| §6 astroid | Strongest empirical comparison arrives after the examples, with no figure placeholder |

---

## Suggested structural edits (do not apply here)

1. Correct §2.3’s plane-curve system to `{f, ∂_y f}`.
2. Reconcile the “six steps” with Figure 1 and with Bertini_real.
3. Add a short curve example (astroid or nodal cubic) to §4, or move the
   astroid comparison out of §6 into §4/§5.
4. Either give the ellipsoid a subsection or demote it to a table-only
   control and stop calling it one of “the three examples.”
5. Define generic projection in §3, not only in the comparison.
6. Re-measure or drop the 3700-line claim.
7. Align timings and the Taubin median residual with
   `paper_artifacts/data/results.json` (see `notes/DATA_CONSISTENCY.md`).
8. Add one sentence that torus numbers are measurements, not test
   invariants.

Priority for a referee-facing pass: **(1), (2), (7)**. The formula
error and the stale numbers are the ones that get caught.
