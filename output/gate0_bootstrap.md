# Gate 0 — full-chain bootstrap SE + in-regime coverage (STOP-IF-FAIL)

Pre-registered in `replication/certification/PREREG_RECERT.md`. Estimator = `sgscatm_chain(refine="frozen_phi")`
run to B-stationarity; SE = `chain_boot_se()` (resample documents → re-run the **whole** chain, incl.
anchors + orientation + frozen refinement, per replicate → perm+sign align each B̂z* to the point
estimate → per-entry SD). Cell: N=200, L=1e4, M=2000, R=30 outer × B=100 inner. Metric: perm+sign to Bz,0.

## Deciding numbers (Gate 0a)

| metric | value | pre-registered PASS band | |
|--------|-------|--------------------------|---|
| coverage of Bz,0 | **0.979** | [0.88, 0.97] | above (conservative) |
| SE_boot / empirical-SD | **1.605** | [0.85, 1.25] | above (over-wide) |
| bias norm | 0.106 | — | (covered by the wide SE) |
| wall-clock | 5546 s (~92 min) | — | 30×101 chain fits, 8 cores |

## Verdict: numeric band **not met**, but the failure mode is the **opposite** of the stop trigger

The pre-registered STOP-IF-FAIL was **"SE ≪ SD (bootstrap does not capture the anchor variance)"**.
The result is **SE ≫ SD**: the full-chain bootstrap **over-captures** the anchor variance
(SE/SD = 1.605) and **over-covers** (0.979 ≥ nominal). So, unlike the split-document jackknife
(SE/SD 0.13, coverage ~0.4 — genuine under-coverage), the bootstrap **does** capture the
anchor-orientation variance and yields **valid, conservative** intervals. It simply over-estimates by
~60%, most plausibly because a minority of bootstrap resamples select a poor anchor set (the anchor
stage is the unstable component), inflating the bootstrap SD above the across-dataset SD; the
point-estimate bias (0.106) is also larger than the clean-data diagnostic (~0.037), consistent with
resample anchor instability.

Because the deviation is **over-coverage** (A5: over-coverage is acceptable, only under-coverage is a
failure), this is best read as a **conservative pass on validity** that **misses the tight calibration
band** — a curable calibration/efficiency refinement, **not** a structural wall and **not** the
"bootstrap can't see the anchor variance" failure the gate was built to catch.

## Consequence for the plan (honoring the pre-registration)

Per the literal STOP-IF-FAIL rule (SE/SD ∉ [0.85,1.25]), the **full Phase 1 re-certification is not
auto-launched**. The decisive scientific question is nonetheless answered:

- **The full-chain bootstrap fixes the jackknife's under-coverage** — it captures the anchor variance
  and covers ≥ nominal. The deliverable feasible SE is therefore **usable but conservative** (~1.6×).
- **Two curable refinements** (either would tighten SE/SD toward 1): (i) an **outlier-robust**
  bootstrap (trim/winsorize resamples whose anchor recovery is poor, e.g. by anchor-TV or by a
  studentized B̂z* screen), or (ii) stabilize the anchor stage on resamples (more EM polish / larger
  min-docfreq). Both are calibration tuning, to be pre-registered separately before use.
- **Cost:** ~90 min per M=2000 cell (100 chain refits/SE at ~10 s/fit). This makes the bootstrap SE
  **expensive** and directly informs the G6 speed framing and the G4 fork: the point+bootstrap path is
  **not** fast (comparable to the competitors), so any speed claim must be scoped to the point estimate
  or the point+sandwich path.

## What still stands (from A3 / DIAG_REPORT, unchanged)

- **At clean (oracle) orientation, in-regime, the PLAIN Lemma-17 sandwich is calibrated** (coverage
  0.934/0.950/0.946, SE/SD≈1). The two-gate regime (Cor. 18) is intact: √M/L→0 AND √M·δ_Φ→0.
- On real CRC (large depth + exclusive taxa) both gates plausibly close, so the **G4 fork** — plain
  sandwich vs bootstrap ratio on CRC (M≈136, where the bootstrap is **cheap**) — remains the decisive,
  low-cost test of whether the analytic sandwich suffices as the paper's primary SE. That test does not
  depend on this gate and can be run independently.

## Bottom line

The full-chain document bootstrap **captures the anchor variance and over-covers** (valid, conservative;
SE/SD 1.61, coverage 0.979) — the opposite of the jackknife's under-coverage. It **misses the tight
pre-registered band** (over-wide, not under), so by the literal STOP-IF-FAIL the full Phase-1 rerun is
held; the honest next step is an outlier-robust bootstrap (pre-registered) and/or the cheap CRC G4 fork.
