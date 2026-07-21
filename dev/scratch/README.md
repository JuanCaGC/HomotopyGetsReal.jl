# Archived scratch validation scripts

These `phase*_check.jl` files were the standalone validation harness used during the Phases 1–6 rebuild. They are **superseded by** the formal suite in `test/`.

Kept for historical reference — assertion-level coverage was migrated per `test/ASSERTION_AUDIT.md`.

`scratch_robust_slice_eagerness_check.jl` (2026-07, Phase 7.5 hardening) is the evidence record for why `_robust_slice_at_z`'s gradient gate runs eagerly on the naive-midpoint attempt: a cropped-`bbox_z` Taubin heart whose naive midpoint is topology-clean but gradient-degenerate (sweep max `|f| ≈ 1.6` if accepted, `2.4e-7` with the eager gate).
