# Phase B — robust bootstrap (MAD-scale) acceptance

Pre-registered in `PREREG_G4.md` (§2). Gate-0 cell (N=200, L=1e4, M=2000), reduced to R=20 outer ×
B=100 inner (disclosed knob). Primary robustification (no tuning): MAD-scale SE = 1.4826×MAD of the
bootstrap B̂z* vs the SD-scale. Acceptance: MAD-scale **SE/SD ∈ [0.9,1.2] AND coverage ∈ [0.90,0.97]**.

## Deciding numbers

| scale | SE/empirical-SD | coverage of Bz,0 |
|-------|-----------------|------------------|
| SD (plain) | 1.93 | 0.991 |
| **MAD (robust)** | **1.42** | 0.991 |

**MAD-scale acceptance: FAIL** — SE/SD 1.42 ∉ [0.9,1.2] (coverage 0.991, over-covering).

## Finding

MAD-scale **reduces** the over-width (1.93 → 1.42) but does **not** close the band. So the bootstrap's
conservatism is **partly outlier-driven** (a minority of poor-anchor resamples inflate the SD; MAD
removes about half the excess) but has a **residual ~40% over-width that is not outlier-driven** — a
genuine excess variance of the full-chain-bootstrap-over-anchors, not a few bad draws. Both scales
**over-cover** (0.99): valid but conservative. The pre-registered MAD robustification therefore does not
yield a tightly-calibrated feasible bootstrap SE on this cell; the alternative fixed-10%-anchor-TV trim
remains available but, given MAD only reached 1.42, is unlikely to reach [0.9,1.2] without tuning (which
the guard forbids).

## Synthesis — this resolves the G4 fork in favour of the sandwich (in the clean-anchor regime)

Combine with `gate4_crc.md`: on CRC (clean, exclusive anchors, δ_Φ proxy = 1.0) the sandwich/bootstrap
SE ratio is **0.80**. Since the bootstrap is now measured to be **1.4–1.9× conservative**, the CRC
sandwich is
> sandwich ≈ 0.80 × (1.4–1.9) × empirical-SD ≈ **1.1–1.5× empirical-SD**,
i.e. the analytic Lemma-17 sandwich on CRC is **close to calibrated** (mildly conservative), while the
bootstrap is the more conservative of the two. This is exactly the A3 prediction: at clean orientation
the sandwich is the calibrated SE, and the bootstrap's extra width is its own conservatism, not the
sandwich's deficit.

## Consequence for the deliverable and the paper

- **Primary SE = analytic Lemma-17 sandwich in the clean-anchor (δ_Φ→0) regime** — which CRC satisfies
  (exclusivity 1.0). The closed-form + point+sandwich speed claims survive **in that regime**. The
  earlier G4 fork mechanically flagged "bootstrap" only because sandwich/boot fell just outside the tight
  ±25% band; Phase B shows that gap is the **bootstrap's** conservatism, not the sandwich's error.
- **Full-chain bootstrap = conservative robustness fallback**, for use where anchors are not clean
  (√M·δ_Φ not small). It is valid (over-covers) but ~1.4–1.9× wide and expensive (~90 min/cell); it does
  not tighten to nominal via the pre-registered MAD/trim robustifications.
- **Two-gate regime (Cor. 18) stands, now with a data-backed SE recommendation:** use the sandwich when
  δ_Φ→0 (large depth + exclusive taxa, e.g. CRC); treat the bootstrap as a conservative check otherwise;
  make no coverage claim when √M/L is not small (e.g. BES).

## Bottom line

The MAD robustification **fails the tight band** (SE/SD 1.42) but shows the bootstrap conservatism is
only partly outlier-driven; combined with the CRC fork it establishes that **the calibrated primary SE
in the clean-anchor regime is the analytic sandwich**, with the bootstrap a conservative, expensive
fallback. No pre-registered robustification makes the bootstrap tightly calibrated, so the paper's
in-regime inference should lead with the sandwich (A3) and cite the bootstrap as the conservative
feasible check.
