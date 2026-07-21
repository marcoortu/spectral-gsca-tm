# Pre-registration — Full re-certification of strada a) (frozen-Φ chain + calibrated SE)

**Written before running.** Two phases; Phase 0 is STOP-IF-FAIL and Phase 1 runs only if it passes.
Nothing tuned; a FAIL is a finding. Package `sgscatm`, `devtools::load_all(".")`.

## Starting point (settled in DIAG_REPORT.md)

- Deliverable estimator = `sgscatm_chain(refine="frozen_phi")` **run to convergence** (B-stationarity);
  in-regime the frozen bias-vs-k curve is flat → no early stop. `joint` = cautionary Prop-19 drift.
- Plain Lemma-17 sandwich is calibrated **at clean orientation** (A3: oracle-frozen@conv + plain,
  L=1e4, coverage 0.934/0.950/0.946, SE/SD≈1) — the analytic SE valid as δ_Φ→0.
- Feasible chain SE is ~6× too small (SE/SD 0.13–0.17); split-document jackknife does **not** repair it
  (halves share frozen Φ̂ → misses anchor re-estimation variance). The deficit is variance, not bias.
  **Calibrated feasible SE = full-chain document bootstrap.**
- Two-gate regime (Corollary 18): calibrated where √M/L→0 (large L) AND √M·δ_Φ→0 (clean anchors).

## Phase 0 — STOP-IF-FAIL: `chain_boot_se()` + in-regime coverage

`chain_boot_se(fit, W, C, B=200, ...)`: resample M documents with replacement → re-run
`sgscatm_chain(refine="frozen_phi")` on the resample (re-running anchors + orientation + refinement,
no cached Φ̂) → perm+sign align each B̂z* to the **point estimate** B̂z → SE = per-entry SD across B.
Unit-tested (`test-boot.R`: SE O(M^{-1/2}); alignment invariance).

- **Gate 0a (quick):** N=200, L=1e4, M=2000; R=30 outer × B=100 inner. Report SE_boot/empirical-SD and
  true coverage (perm+sign to Bz,0). **PASS** = SE/SD∈[0.85,1.25] AND coverage∈[0.88,0.97] → proceed.
  **FAIL** = SE≪SD → STOP, write finding (check the bootstrap re-runs the anchor stage).
- **Gate 0b (full, if 0a passes):** N=200, L=1e4, M∈{2000,5000}; R=50 × B=200. PASS = coverage∈[.88,.97]
  and SE/SD∈[.9,1.2] at M≥2000; over-coverage OK, only under-coverage fails. Report M=1000 with the
  bad-anchor-bias caveat (orientation-bias-limited, out of regime, not a bootstrap failure).

Output `output/gate0_bootstrap.md`. Phase 1 runs only if Gate 0a passes.

## Phase 1 — full re-certification (only if Gate 0a passes)

Package: `refine="frozen_phi"` documented default; `joint` cautionary; primary SE = `chain_boot_se`
(feasible) with `vcov.sgscatm_chain` (plain sandwich) alongside; pilot + `sgscatm_vcov` retained as
secondaries. Remove dead duplicate `sgscatm()` in `egscatm_fit.R` (logged). R CMD check clean.

- **A Simulation** (frozen-Φ@conv, perm+sign): T1 subspace/collapse (correct Table-1 to per-direction
  cosines, not <1e-5 angle); T2 coverage vs √M/L @L=200 (incidental gate; bootstrap calibrates
  variance not the 1/L bias); T2-inreg = Gate 0b; T3 crossover vs STM (frozen-Φ, in-regime); b(z)/G0.
- **B CRC** (swap pilot→frozen-Φ chain): T4 delocalization; T5 Wald vs PERMANOVA; **T6 speed reported
  honestly** — point (fast) and point+SE (bootstrap path comparable to competitors; sandwich path
  fast); Fig2 forest (primary-SE CI); **Fig3 G4 fork = sandwich vs bootstrap ratio on CRC**; Fig4
  known genera.
- **C BES** (swap → frozen-Φ chain): T7/Fig5 recovery strong/weak vs STM + runtime; out of regime, no
  coverage claim.

## G4 fork (decides primary SE)

CRC plain-sandwich/bootstrap SE ratio: ∈[0.8,1.25] → **primary SE = analytic sandwich** (fast; speedup
holds via point+sandwich), bootstrap = robustness. Diverges → **primary SE = bootstrap**; speed claim
scoped to the point estimate. Record which branch the data select.

## Protocol

Minimal edits; reuse `ilr_se` loop, `sim_dgp.R`, `perm_sign_*`, phase2/BES scripts. Fixed seeds
base+regime·1000+rep, PSOCK. Monotone/stability checks. Compute knobs if bound: inner B 200→100,
outer R 50→30 (disclosed).
