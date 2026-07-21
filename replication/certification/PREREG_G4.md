# Pre-registration — G4 CRC fork + robust-bootstrap plan

**Written before running.** Phase A (CRC fork) runs now; Phase B (robust bootstrap) is frozen and
runs only after A. Nothing tuned; a FAIL is a finding. Package `sgscatm`, `devtools::load_all(".")`.

## Context (gate0_bootstrap.md / DIAG_REPORT.md)

Full-chain bootstrap captures anchor variance but is conservative (SE/SD 1.61, cov 0.979 in sim);
plain Lemma-17 sandwich is calibrated at clean orientation (A3, SE/SD≈1). Two-gate regime (Cor. 18):
sandwich valid as δ_Φ→0, bootstrap otherwise. The G4 CRC fork decides which is the paper's **primary**
SE on real data, where M≈136 makes the bootstrap cheap and taxa are exclusive (anchors plausibly clean).

## Phase A — G4 CRC fork (run now)

Estimator swapped to `sgscatm_chain(refine="frozen_phi")` on the CRC genus path. Data: Zeller 2014
via SIAMCAT (`output/phase2_data.rds`, genus prev≥10%, renormalized). Counts reconstructed as
`round(rel_abund × depth)`, depth=1e6 (chain anchors need counts; large depth ⇒ in-regime). K=5.
Covariates: study_condition, age(std), BMI(std), gender.

Steps: fit chain → B̂z; SE_sandwich = `vcov.sgscatm_chain`; SE_boot = `chain_boot_se(B=200)`
(perm+sign to point est). Report per-entry ratio SE_sandwich/SE_boot: **median, per-covariate median,
% within ±25% of parity**; both SE magnitudes and B̂z. δ_Φ proxy = anchor exclusivity/TV of Φ̂.
Significant covariates (joint Wald) under **both** SEs (sandwich full vcov; bootstrap empirical vcov).

**Frozen fork rule:** median(SE_sandwich/SE_boot) ∈ [0.8,1.25] AND ≥60% entries within ±25%
→ **primary SE = analytic Lemma-17 sandwich** (speed/closed-form survive; bootstrap = robustness).
Else → **primary SE = bootstrap**; sandwich = asymptotic (δ_Φ→0); speed scoped to point estimate.
Emit `output/gate4_crc.md`.

## Phase B — robust bootstrap (frozen; run only if bootstrap is primary)

Pre-registered a priori (not swept): (primary) **MAD-scale** SE = 1.4826×MAD of bootstrap B̂z*
vs SD-scale; (alternative) **fixed 10% anchor-TV trim** then SD. Acceptance on the Gate-0 cell
(N=200,L=1e4,M=2000): SE/SD ∈ [0.9,1.2] AND coverage ∈ [0.90,0.97]; over-coverage OK. Report
SD-scale, MAD-scale, empirical SD side by side. If it passes → re-run Gate 0 (M∈{2000,5000}); if that
passes the tight band → clear full Phase-1 re-certification. **Guard:** MAD constant and 10% trim are
fixed a priori; do not sweep to hit the band. Emit `output/robust_bootstrap.md`.
