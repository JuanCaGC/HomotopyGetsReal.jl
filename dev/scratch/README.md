# Archived scratch validation scripts

These `phase*_check.jl` files were the standalone validation harness used during the Phases 1–6 rebuild. They are **superseded by** the formal suite in `test/`.

Kept for historical reference — assertion-level coverage was migrated per `test/ASSERTION_AUDIT.md`.

`scratch_robust_slice_eagerness_check.jl` (2026-07, Phase 7.5 hardening) is the evidence record for why `_robust_slice_at_z`'s gradient gate runs eagerly on the naive-midpoint attempt: a cropped-`bbox_z` Taubin heart whose naive midpoint is topology-clean but gradient-degenerate (sweep max `|f| ≈ 1.6` if accepted, `2.4e-7` with the eager gate).

`scratch_stage4c_historical_curves.jl`, `scratch_stage4c_taubin_resolve_trace.jl`, `scratch_stage4c_taubin_validation.jl`, `scratch_whitney_umbrella_tip.jl`, and `scratch_whitney_umbrella_validation.jl` (2026-07, isosingular deflation Stage 4c) are the investigation records `test/test_isosingular_deflation.jl` cites directly by filename as the basis for its own historical-curves and Whitney-umbrella testsets — not independently re-derived there.

`scratch_torus_validation.jl` and `scratch_cone_validation.jl` (2026-07) are the `decompose_3d_surface` validation trials behind `docs/DESIGN_NOTES.md`'s backlog entries on the torus's positive-dimensional fold locus and the cone's undetected-apex failure — both surfaces failed distinctly; neither is a usable fixture yet.
