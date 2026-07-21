# Pre-registration — Phase 1 full re-certification (frozen-Φ chain + analytic sandwich)

**Written before running.** Nothing tuned; a FAIL is a finding. Paper edits deferred to a later step.
Package `sgscatm`, `devtools::load_all(".")`. Metric: perm+sign align B̂z to Bz,0 (never Procrustes).

## Deliverable (fixed by DIAG_REPORT/gate4/robust_bootstrap)

Estimator = `sgscatm_chain(refine="frozen_phi")` to B-stationarity. Primary SE = analytic Lemma-17
sandwich (`vcov.sgscatm_chain`), calibrated in the clean-anchor regime δ_Φ→0 (A3; CRC δ_Φ proxy 1.0).
Fallback = full-chain bootstrap (`chain_boot_se`), conservative (1.4–1.9×). Two-gate regime (Cor. 18):
Wald calibrated when √M/L→0 AND √M·δ_Φ→0.

## Package fixes (minimal, logged)

- Small-M B-stationarity: add an F-convergence backstop to `.sg_refine` so the stop fires cleanly at
  small M (Gate-4 caveat: M=136 ran to cap). Leave core estimator behavior intact.
- Remove the dead duplicate `sgscatm()` in `egscatm_fit.R` (canonical def is `sgscatm_fit.R`); repoint
  the two replication scripts that source it. document(); tests.

## Gates / tables

- **T1** subspace recovery + raw-ratio M^{-1/2} collapse + chain RMSE→oracle. Report per-direction
  principal **cosines** (correcting the false <1e-5 angle). M∈{1000,2000,5000}.
- **T2** coverage vs √M/L, fixed L=200, M∈{1000,2000,4000,5000}, R≥50, **clean anchors** (exclusive
  topics); report **feasible AND oracle** coverage side by side. SE=sandwich. Expect nominal at small
  √M/L, degrading as it grows (incidental gate; sandwich variance does not remove the 1/L bias).
- **T2-inreg-feasible (KEY):** feasible frozen-Φ + sandwich, **L=1e4**, genuinely clean anchors
  (near-one-hot, N=100, M∈{1000,5000}, exclusive topics), R≥50, vs known Bz,0. **PASS** = coverage
  ≥0.90 (over-coverage OK) AND SE/SD∈[0.9,1.6]; **under-coverage (SE/SD<0.9, coverage<0.90) = FAIL**
  (sandwich misses variance even with clean anchors → bootstrap required in-regime). Bootstrap reported
  alongside as the conservative check.
- **T3** crossover vs STM, frozen-Φ, in-regime, M∈{1000,5000}, weak/strong: feasible/oracle/STM.
- **b(z)/G0 (G8):** round-5 harness; MC+Richardson vs closed-form b(z); CRN finite-diff vs G0.
- **CRC (G1/G4/G5/G6):** chain frozen-Φ, sandwich primary. T4 delocalization (genus+species); T5 Wald
  vs PERMANOVA (BMI SE-sensitive; study_condition robust); T6 speed three paths (point / point+sandwich
  fast / point+bootstrap slow); Fig2 forest (sandwich CI); Fig3 sandwich-vs-bootstrap (Cor-18 gate-2
  illustration); Fig4 known genera. G4 from gate4_crc.
- **BES (T7/Fig5):** frozen-Φ chain vs STM, recovery strong/weak + runtime; out of regime (√M/L≈86),
  no coverage claim.

## Certification record

Update MANIFEST.md + GATE_REPORT.md: every Table/Figure → script → gate → chain-certified vs draft
number, flag what moved; one line per gate PASS/FAIL with deciding number. Short pointer to sections
whose numbers moved (full PAPER_EDITS deferred).

## Protocol / compute

Minimal edits; reuse sim_dgp, perm_sign_*, phase2/BES scripts, round-5 harness. Fixed seeds
base+regime·1000+rep, PSOCK. Coverage tables closed-form (cheap); drivers = STM (T3/BES) + b(z)/G0 MC;
bootstrap only for CRC Fig3 (have it) + one T2-inreg cell. Disclose any reduced-rep knob. No tuning.
