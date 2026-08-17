# Data and figure consistency

**Sources read, not edited:** `paper/HomotopyGetsReal_paper_current.md`
against `paper_artifacts/data/results.json`. Figure files checked for
existence under `paper_artifacts/figures/` (layout in
`figures_organization.md`). Live `src/` / `test/` used only to confirm
claims that are not in `results.json`.

`results.json` mixes **current figure measurements** (sphere, ellipsoid,
Taubin, torus, curves, pipeline steps) with a **stale test-count
block** (`tests.fast` / `tests.full`, dated by its own note to the
artifact-generation session that added `[extras] Test`). Do not quote
the `tests` object as the suite size.

Rounding: scientific notation in the paper is compared to JSON at the
digits the paper actually prints. “Matches” means the printed digits
are a correct rounding of the JSON value.

---

## 1. Residual table (paper Table 3 vs JSON)

### Sphere production (`results.json` `sphere.default`)

| Stat | Paper | JSON | Verdict |
|---|---|---|---|
| n | 19,504 | 19504 | match |
| mean | 2.85×10⁻⁸ | 2.850793…e-8 | match |
| median | 2.62×10⁻⁸ | 2.622716…e-8 | match |
| p₉₀ | 5.32×10⁻⁸ | 5.318477…e-8 | match |
| p₉₉ | 7.05×10⁻⁸ | 7.049153…e-8 | match |
| max | 7.58×10⁻⁷ | 7.584689…e-7 | match |
| min | (caption: some exact 0s) | 0.0 | match |

### Sphere coarse (`sphere.test`)

| Stat | Paper | JSON | Verdict |
|---|---|---|---|
| n | 152 | 152 | match |
| mean | 2.31×10⁻⁸ | 2.313642…e-8 | match |
| median | 2.14×10⁻⁸ | 2.136178…e-8 | match |
| p₉₀ | 4.92×10⁻⁸ | 4.917831…e-8 | match |
| p₉₉ | 5.61×10⁻⁸ | 5.610274…e-8 | match |
| max | 5.61×10⁻⁸ | 5.610274…e-8 | match |

### Ellipsoid production (`ellipsoid.default`)

| Stat | Paper | JSON | Verdict |
|---|---|---|---|
| n | 19,500 | 19500 | match |
| mean | 2.70×10⁻⁷ | 2.699127…e-7 | match |
| median | 2.79×10⁻⁸ | 2.790703…e-8 | match |
| p₉₀ | 5.85×10⁻⁸ | 5.851432…e-8 | match |
| p₉₉ | 8.99×10⁻⁸ | 8.992520…e-8 | match |
| max | 1.12×10⁻⁴ | 1.116566…e-4 | match |

Paper §5.1 correctly notes this max is **above** the test-suite gate
`all(<=(1e-4))` used on the coarse ellipsoid mesh
(`test/test_surfacedecomposition.jl:263` as cited in the technical
summary). That gate is **not** applied to this production mesh. Do not
smooth the 1.12×10⁻⁴.

### Ellipsoid coarse (`ellipsoid.test`)

| Stat | Paper | JSON | Verdict |
|---|---|---|---|
| n | 160 | 160 | match |
| mean | 5.26×10⁻⁷ | 5.258816…e-7 | match |
| median | 2.57×10⁻⁸ | 2.573162…e-8 | match |
| p₉₀ | 5.96×10⁻⁸ | 5.960512…e-8 | match |
| p₉₉ | 1.000×10⁻⁵ | 1.000264…e-5 | match (printed as 1.000×10⁻⁵) |
| max | 1.000×10⁻⁵ | 1.000264…e-5 | match |

### Taubin heart (`taubin.paper`) — **one real error**

| Stat | Paper | JSON | Verdict |
|---|---|---|---|
| n | 1,638 | 1638 | match |
| mean | 1.034×10⁻⁷ | 1.033625…e-7 | match |
| **median** | **2.31×10⁻⁸** | **3.437564…e-8** | **MISMATCH** |
| p₉₀ | 1.847×10⁻⁷ | 1.847077…e-7 | match |
| p₉₉ | 2.000×10⁻⁶ | 2.000157…e-6 | match |
| max | 2.42×10⁻⁶ | 2.420909…e-6 | match |

The printed median is not a rounding of 3.44×10⁻⁸. It is closer to the
sphere’s median band. **Replace Table 3 Taubin median with 3.44×10⁻⁸**
(or 3.437×10⁻⁸). The §5.1 sentence “median residual clusters near
2×10⁻⁸” is still roughly true after the correction (3.4×10⁻⁸), but
do not leave 2.31 in the table.

### Torus (`torus`) — characteristic, mostly consistent

Paper marks the row with `*` and says median ~1.4×10⁻⁶, max ~2.0×10⁻⁵.

| Stat | Paper | JSON | Verdict |
|---|---|---|---|
| n | 78,340 | 78340 | match (this run) |
| mean | 2.211×10⁻⁶ | 2.196845…e-6 | slight: paper is 2.211, JSON 2.197 |
| median | ~1.4×10⁻⁶ | 1.367701…e-6 | OK as ~ |
| p₉₀ | 5.491×10⁻⁶ | 5.469955…e-6 | slight (paper 5.491 vs 5.470) |
| p₉₉ | 1.008×10⁻⁵ | 1.001707…e-5 | slight |
| max | ~2.0×10⁻⁵ | 1.782241…e-5 | OK as ~; JSON is 1.78×10⁻⁵ |

Paper §4.3 triangle count **160,086**; JSON `n_mesh_triangles` is
**160381**. Same section also lists other observed pairs
(78,360/168,962 and 78,339/160,346). 160,086 is **not** the JSON run
tied to Figure 4’s current `torus.pdf` (the JSON `figure_pdf` points at
the file next to these 160381 triangles). Either reprint 160,381 as
“the run behind Figure 4” or stop claiming Figure 4 and Table 3 are the
same decompose.

Naked edges: JSON 18, with an explicit non-determinism note (21 vs 144
across runs). Paper §4.3 does not quote a naked-edge number (good).

---

## 2. Combinatorial counts and timings

### Sphere §4.1

| Claim | Paper | JSON | Verdict |
|---|---|---|---|
| Vertices / edges / faces | 2 / 2 / 2, all Critical | 2 Critical, 2, 2 | match |
| Mesh | 19,504 verts, 39,004 tris | 19504 / 39004 | match |
| Coarse mesh | 152 / 300 | 152 / 300 | match |
| Warm time | **3.94 s** | **0.957 s** (`default`) | **MISMATCH** |
| First-call | 23.93 s | 26.65 s is the *coarse* `test` run | different experiment; do not equate |

The 3.94 s figure is not in `results.json`. Do not present 3.94 s and
the Table 3 mesh as one measurement without a matching log.

### Taubin §4.2

| Claim | Paper | JSON | Verdict |
|---|---|---|---|
| z_crit ≈ {−1, 1, 1.0648, 1.2367} | yes | −1.0000000357, 0.9999999344, 1.0647678179, 1.2366591700 | match at 10⁻⁷ as claimed |
| Plain vertices | 14 = 10 Critical + 4 Artificial, 0 Singular | same | match |
| Edges / faces | 14 / 14 | 14 / 14 | match |
| Mesh | 1,638 / 3,124 | 1638 / 3124 | match |
| Overlay | 20 = 14 C + 4 A + 2 S | `taubin_singular_structure.paper_figure` overlay same | match |
| `compute_critical_z_slices` | **18.52 s** (called “first-call”) | **0.159 s** | **MISMATCH** (JIT vs warm, or different session) |
| `decompose_3d_surface` | **24.86 s** first-call | **34.24 s** | **MISMATCH** |
| Incidence overlay wall | (not in §4.2) | 47.96 s | paper silent; fine |

Retry ladder §5.2: naive z_mid = 0 rejected; accepted 0.06 after 5
retries; 4 Artificial + 2 Singular on the naive slice. JSON
`robust_slice_slab_m1_p1` matches all of that.

### Torus §4.3

8 Critical, 0 Singular, 8 edges, 8 faces: JSON match. Mesh vertices
78,340 match; triangles do not (above).

### Curves (`curve_examples`) — paper §6 vs JSON

Paper does not put these in Table 3; §6 astroid: 4 cusps, 4 edges.
JSON: astroid 4V/4E all Singular; cusp 3V/2E; nodal cubic 4V/4E
(1 Critical, 1 Singular, 2 Boundary). Matches the technical summary
and the test-encoded nodal cubic. Node figure is **illustrative**
(JSON `note`); paper never claims an automated node decomposition
(good — the technical summary records that `decompose_1d_curve` finds
0 critical vertices there).

### `sample_edge` regression (paper §5.1)

Paper: max residual 0.4998 → 7.3×10⁻⁷. JSON
`pipeline_steps_curve.step8_before_after`: **0.499781… → 7.33857…e-7**.
Match.

### Stitching (paper §3.9 vs JSON)

Paper: 188 naked before incidence; 31–35 after, run-to-run jitter.
JSON `pipeline_steps_surface.stitching`: **188 → 34**. Compatible
(34 is inside 31–35). Talk slide 16 prints 34 as if exact; see
`notes/TALK_ALIGNMENT.md`.

### Test suite size

| Source | Number |
|---|---|
| Paper §5.3 | 537/537, ±1 jitter |
| Talk slide 20 | **538/538** |
| `results.json` `tests.fast` | **162/162** (1m40.4s) |
| `results.json` `tests.full` | **178/178** (3m15.1s) |
| `CLAUDE.md` / `README.md` | 537 ±1; fast 477/477 |

The JSON `tests` block is **stale by an order of magnitude** and still
unmarked except for a Project.toml note. The technical summary already
flags this (`REPO_TECHNICAL_SUMMARY.md:673`). **Do not quote
`results.json` tests.** Paper 537 vs talk 538 is the live ±1 issue,
not a JSON issue.

Paper “thirteen test files”: `test/` has 12 `test_*.jl` plus
`runtests.jl`. Imprecise, not a numerical contradiction.

---

## 3. Figures: placeholders vs files

Paper placeholders (binary figures not inlined in the markdown):

| Placeholder | Intended file | On disk? | Used in talk? |
|---|---|---|---|
| Figure 1 pipeline TikZ | in `talk.tex` `\PipelineDiagram`, not a PDF | n/a (talk has it; paper does not) | yes, slide 7 |
| Figure 2 sphere | `figures/paper/surfaces/sphere.pdf` | yes | no |
| Figure 3 Taubin | `taubin.pdf` (mesh-only) | yes | slides 13, 17 |
| — (text: 20-vertex overlay) | `taubin_singular_structure.pdf` | yes | slide 18 |
| Figure 4 torus | `torus.pdf` (~81 MB) | yes | talk uses `talk/torus_slide_raster.png` instead |
| Figure 5 histograms | `residual_histograms.pdf` | yes | no |
| Figure 6 sample_edge | `pipeline_steps_curve/08_sample_edge_fix_before_after.pdf` | yes | **not** in current `talk.tex`; yes in lean deck |

**Missing from the paper as placeholders, present on disk and in the
talk or table:**

- `ellipsoid.pdf` — Table 3 row, no Figure.
- `curves/astroid.pdf` and `astroid_bertini_real_sampled.pdf` — §6
  comparison, no Figure. Talk slides 9–10 use both.
- Pipeline-step sequences (curve 01–06, surface 01–07) — talk slides
  11–16; paper never claims them.

**Stale path inside JSON (not a paper bug):**
`pipeline_steps_curve.output_dir` still says
`…/figures/pipeline_steps_curve` (pre-reorg). Files actually live in
`figures/paper/pipeline_steps_curve/`. `figure_pdf` fields for the
current surfaces were rewritten to the new layout.

**JSON `figures.count`: 4** and a note “Sections 1–4 request 4 PDFs
(not 5)” — that is the original generation request, not the paper’s
Figure 1–6 plan. Do not treat `figures.count` as the paper figure list.

---

## 4. Unverified or paper-only numbers (not in JSON)

These cannot be checked against `results.json`. They are not therefore
false; they are **unverified by this artifact file**.

| Claim | Where | Check elsewhere? |
|---|---|---|
| HomotopyGetsReal ≈ 3700 lines | §3 | Live `wc -l src/*.jl` = **5867**. Unverified in JSON; **false vs live tree** |
| 22 Taubin deflation firings, 17×1 round, 5×4 rounds, sequence `[2,1,1,1,0]` | §3.5 | Not in JSON. Technical summary / DESIGN_NOTES may hold it; re-verify before submission |
| False-positive rate 4/6 on randomized reduction | §3.5 | Not in JSON |
| Sphere warm 3.94 s / first-call 23.93 s | §4.1 | Not in JSON (JSON 0.96 s) |
| Taubin 18.52 s + 24.86 s | §4.2 | Not in JSON (JSON 0.16 s + 34.24 s) |
| Fixed-axis torus residuals ~10⁸ | §4.3, §6 | Not in JSON (only the `:random` run is stored). DESIGN_NOTES records the failure; do not cite 10⁸ from JSON |
| Bertini_real astroid 32 vertices / 20 edges, filtered 4 cusps / 6 edges | §6 | Coordinate comparison lives in the **external** Auditor workspace (`CLAUDE.md` mechanism 3), not in this repo’s JSON |
| Bertini_real Griffis–Duffy 4,1,1 / 588 minors | §6 | Paper already says do not cite as independently verified |
| 537/537 live full suite | §5.3 | JSON `tests.full` is 178; use a fresh `Pkg.test()` log, not JSON |

---

## 5. Priority list (do not apply here)

1. **Correct Table 3 Taubin median** 2.31×10⁻⁸ → 3.44×10⁻⁸.
2. **Do not quote `results.json` `tests` (162/178).**
3. **Reconcile torus triangle count** 160,086 vs JSON 160,381, or stop
   tying Figure 4 to Table 3.
4. **Reconcile or drop §4 wall-clock times** (sphere 3.94 s vs 0.96 s;
   Taubin 24.86 s vs 34.24 s; crit-z 18.52 s vs 0.16 s). If they are
   first-call vs warm, say so *and* print the JSON warm numbers too.
5. **Add figure placeholders** for astroid side-by-side and ellipsoid,
   or stop calling them paper figures.
6. **Re-measure src line count** or delete “3700 lines.”
7. Pin talk 538 vs paper 537 with a live run (`notes/TALK_ALIGNMENT.md`).

Items 1–4 are numerical. Item 1 is the only Table-3 cell that is
simply wrong.
