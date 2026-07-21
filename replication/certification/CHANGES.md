# CHANGES — package realignment to the Theorem-16 chain

## What changed vs the raw pilot

**New package estimator + inference (all in new `R/` files, promoted verbatim up to
renaming from the validated replication code; `sgscatm()` pilot left intact):**

- `R/refine.R` — exact-objective engine + damped per-document Gauss-Newton LS z-step
  (`.sg_z_step`, trust-region cap `dz_cap`), read-out (`.sg_readout_gn`), and the k-step
  refinement loop with B-stationarity stop (`.sg_refine`). Promoted from
  `replication/basin_check/01_functions.R` (derivative-checked there) and
  `replication/feasibility/02_constrained_refine.R`.
- `R/anchors.R` — anchor-word recovery (`.sg_anchor_pipeline`) + pooled multinomial-EM
  polish (`.sg_phi_em_polish`, the single deviance step). From `feasibility/01_anchors.R`.
- `R/orient.R` — general-linear alignment (`.sg_gl_align`), never Procrustes.
- `R/chain.R` — `sgscatm_chain()` orchestrator (Stage 1 pilot → Stage 2 anchored GL
  orientation → Stage 3 k-step LS refinement) + `print`/`coef` methods.
- `R/sandwich.R` — Lemma-17 sandwich `vcov.sgscatm_chain()` (from `fs_sandwich`),
  `chain_se()`, and permutation+sign metric helpers `perm_sign_align()`,
  `perm_sign_coverage()` (replacing the old Procrustes alignment for reported RMSE/coverage).
- `tests/testthat/test-{chain,sandwich,anchors}.R`.

`sgscatm()` (the raw pilot) is unchanged and retained as a documented fast/secondary
estimator; `sgscatm_vcov()` and `ilr_se()` retained as reference SEs for G4.

## Key deviation from the brief, with evidence (logged per protocol)

**The default refinement is `refine = "joint"` (the V4 variant), not the brief's
`frozen_phi` baseline.** Direct measurement (`replication/certification/diag_refine.R`,
M=1000, N=500, K=5, b_max=0.5, L=400):

| refinement | RMSE(Bz) sweep 1 | 3 | 5 | 10 | 20 |
|-----------|:---:|:---:|:---:|:---:|:---:|
| **frozen_phi** | 0.250 | 0.300 | 0.322 | 0.353 | 0.384 |
| **joint (V4)** | 0.409 | 0.396 | 0.365 | 0.307 | 0.239 |

`frozen_phi` is **best after ~1 sweep and then monotonically degrades** — a direct empirical
confirmation of **Proposition 21** (the fixed-Φ criterion optimum is inconsistent at fixed L).
The B-stationarity rule (`|ΔB|/rms(B) < 1e-3`) does not fire before this degradation, so a
frozen-Φ chain run to the rule/cap returns a drifted B. `joint` re-estimates Φ each sweep,
is stable, improves monotonically, and is the estimator the manuscript's reported feasible
numbers actually come from (`replication/feasibility/04_run.R` chooses V4). It corresponds to
the Prop-19 alternating mode.

**Consequence for the manuscript:** Section 3.7's "frozen Φ̂, no Φ update" prose does not match
the implemented (and better) estimator. Either reconcile the text to the joint refinement, or
present frozen-Φ with the required aggressive early stop (k≈1) and apply the Lemma-17 sandwich
at that early iterate. See `PAPER_EDITS.md`.

`refine = "frozen_phi"` remains available and its Prop-21 behaviour is unit-tested
(`test-chain.R`).

## Open item flagged during certification

- **G0 principal-angle sub-metric did not reproduce machine-precision recovery** (measured
  ~1 rad at the base design L=200, vs the manuscript's <1e-5 in Table 1). The other two G0
  sub-metrics DID reproduce: the raw-pilot norm ratio collapses as M^{-1/2}
  (0.060→0.042 from M=1000→2000), and the chain RMSE descends toward the oracle floor. The
  angle discrepancy is under investigation (candidate cause: at λ=1 the pilot subspace mixes
  the covariate block `√λ Q_C`, so the document-score subspace comparison needs either the
  λ=0 pilot or the word-side signal subspace as the target). Reported as a finding, not tuned.

## Reproduce

```
Rscript replication/certification/cert_sim.R          # G0, G7, G2c (reduced reps)
devtools::test()                                       # chain / sandwich / anchor unit tests
```
