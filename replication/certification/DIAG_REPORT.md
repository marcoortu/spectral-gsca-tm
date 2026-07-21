# DIAG_REPORT — is the chain's coverage failure curable (A') or structural (B')?

Pure simulation (`sim_dgp.R`, `alpha_beta=0.05`, K=5, P=4, b_max=0.5, σ_ε=0.3, NegBin large L).
Pre-registered in `PREREG_DIAG.md` (+ amendments A1–A6). Metric: perm+sign alignment to Bz,0
(never Procrustes). Scripts: `diag_bias.R` (oracle Part 1), `diag_feasible.R` (A1-corrected
feasible Part 2/3), `diag_figures.R`. Data: `output/diag_P1_oracle.rds`, `output/diag_P2_feasible.rds`.

## Bottom line

**The residual bias is the alternating (joint) drift, not the structural 1/L wall → A' on the bias
axis.** In-regime (L=1e4), oracle-frozen-Φ **at convergence** + the **plain** Lemma-17 sandwich
already covers nominally — no early stop and no jackknife are needed when the orientation is clean.
**But the *feasible* (anchored) chain is not calibrated**: its cross-replicate spread is dominated by
**anchor-orientation variance** that the plain sandwich misses (SE/SD ≈ 0.15) and that the
**split-document jackknife does not repair** (SE/SD 0.13→0.14; even bias-removed coverage fails). The
deliverable calibrated SE therefore requires a **full-chain document bootstrap** over the anchor
stage, not the split-document jackknife — a feasibility-of-anchors task, not a structural wall.

## A3 lead readout — oracle-frozen@convergence (k=10) + PLAIN sandwich, L=1e4

| M | bias norm | SE/SD | true coverage |
|---|---|---|---|
| 1000 | 0.0063 | 0.964 | 0.934 |
| 2000 | 0.0034 | 1.030 | 0.950 |
| 5000 | 0.0035 | 1.003 | 0.946 |

Coverage ∈ [.90,.97] and SE/SD ≈ 1 at every M, **even at k=10 (convergence)**. This is exactly the
estimator Lemma 17 is designed for, and it covers nominally → **not the suspected-bug branch** (bug
gate: SE/SD≈1 with nominal coverage; perm+sign verified self/relabel → 1e-33; sandwich unit-tested
vs the explicit ⊗ loop). The plain sandwich is correct and calibrated when orientation is clean.

## Part 1 — oracle bias isolation (N=500)

Frozen-Φ (true Φ₀) k-curve vs joint@convergence, L-contrast. Figure: `output/figures/diag_kcurve.pdf`.

| L | M | frozen@1 | frozen@10 | joint@conv | | SE/SD f10 | cov f10 | cov joint |
|---|---|---|---|---|---|---|---|---|
| 1e3 | 1000 | 0.0129 | 0.0173 | 0.0661 | | 0.955 | 0.925 | 0.762 |
| 1e3 | 2000 | 0.0104 | 0.0147 | 0.0875 | | 0.984 | 0.931 | 0.526 |
| 1e3 | 5000 | 0.0122 | 0.0174 | 0.1528 | | 0.991 | 0.863 | 0.236 |
| 1e4 | 1000 | 0.0063 | 0.0063 | 0.0608 | | 0.964 | 0.934 | 0.709 |
| 1e4 | 2000 | 0.0031 | 0.0034 | 0.0925 | | 1.030 | 0.950 | 0.457 |
| 1e4 | 5000 | 0.0030 | 0.0035 | 0.1619 | | 1.003 | 0.946 | 0.188 |
(bias norm columns; ‖Bz,0‖≈1.1)

**A4 L-contrast (the discriminator).** The frozen bias-vs-k curve is **flat at L=1e4**
(0.0063→0.0063 at M=1000; 0.0031→0.0034 at M=2000) and mildly rising at L=1e3
(0.0104→0.0147). The residual frozen bias is the **1/L incidental term, cured by depth** — so
**in-regime no early stop is needed** (frozen@1 ≈ frozen@10).

**The joint drift is the culprit and it is structural to *alternating*, not to depth.** joint@conv
bias is **L-insensitive** (0.088 at L=1e3 vs 0.093 at L=1e4, M=2000) and **grows with M**
(0.061→0.093→0.162 at L=1e4), so its coverage **collapses** as √M amplifies it (0.71→0.46→0.19).
This is exactly the Prop-19 alternating drift — re-estimating Φ each sweep accumulates per-document
noise into Φ — and it is what sank the earlier G2c oracle runs (which used `refine="joint"`).

## Part 2/3 — feasible chain (anchors), frozen@k*=1, N=200, L=1e4 (A1-corrected)

Point estimate is always the chain B̂z (A1); jackknife contributes only `var_add`. `B_jk`-centered
is a labelled SECONDARY column (outside the verdict). 0/180 monotonicity violations.

| M | SE variant | bias norm | SE/SD | true cov | bias-removed cov |
|---|---|---|---|---|---|
| 1000 | plain | 0.230 | 0.134 | 0.366 | 0.244 |
| 1000 | jk-inflated | 0.230 | 0.138 | 0.374 | 0.255 |
| 1000 | *sec (B_jk)* | *0.134* | *0.201* | *0.480* | *0.406* |
| 2000 | plain | 0.037 | 0.154 | 0.455 | 0.381 |
| 2000 | jk-inflated | 0.037 | 0.158 | 0.464 | 0.389 |
| 2000 | *sec (B_jk)* | *0.029* | *0.251* | *0.510* | *0.509* |
| 5000 | plain | 0.039 | 0.164 | 0.379 | 0.448 |
| 5000 | jk-inflated | 0.039 | 0.167 | 0.385 | 0.451 |
| 5000 | *sec (B_jk)* | *0.042* | *0.216* | *0.389* | *0.506* |

**Finding.** The feasible chain's sandwich SE is **~6× too small** vs the empirical cross-replicate
SD (SE/SD 0.13–0.17), and the split-document **jackknife inflation barely moves it** (0.134→0.138):
its two halves share the *fixed* anchored Φ̂, so `var_add` captures only count/refinement noise, not
the **anchor re-estimation variance** that dominates when the anchor stage is re-run on a fresh
sample. **Bias-removed coverage also fails** (0.24–0.45), so the deficit is **variance, not bias** —
an anchor-orientation variance the closed-form SE and split-jackknife both omit. (The M=1000 bias
0.230 is inflated by a few bad-anchor replicates at small M/N; at M≥2000 bias is small ≈0.04.)

## §5 verdict — deciding numbers

- **N1** oracle-frozen@k* bias (L=1e4, M=2000) = **0.0031** (tiny); decreases with M
  (0.0063→0.0031→0.0030), and the L-contrast vs joint is decisive: joint bias is L-insensitive
  (0.087↔0.093) and M-*increasing* → the residual is the **alternating drift**, not 1/L.
- **N2** oracle-frozen@k* true coverage (L=1e4) = **0.954** ∈ [.90,.97] ✓.
- **N3** feasible-frozen@k*+jk **bias-removed** coverage (L=1e4) = **0.255 / 0.389 / 0.451**
  (M=1e3/2e3/5e3) ✗ — variance not calibrated.
- **N4** feasible-frozen@k*+jk **true** coverage (L=1e4) = **0.374 / 0.464 / 0.385** ✗.

**Mapping to §5:** the bias question is settled — **A' (curable, alternating drift)**: N1 small, N2
nominal, oracle-frozen@k* clean, k-curve flat at L=1e4. But **A' is *not* fully viable via the
promoted split-document jackknife**, because N3 and N4 fail: this is §5's **"third finding"** —
oracle-frozen@k* is clean while N4 is blocked by a **feasible anchor-stage residual that is variance,
not bias**. Per A5, the trend over M does **not** close it at N=200 (true cov flat ≈0.37–0.46 through
M=5000; SE/SD 0.13→0.16). The block is a *feasibility-of-anchors / SE-estimator* issue, not a
structural wall.

## A6 conclusion — is early stop / jackknife needed in-regime?

- **At clean (oracle) orientation, in-regime (L=1e4): neither is needed.** Frozen-Φ **at
  convergence** + the **plain** Lemma-17 sandwich covers nominally across M. Early stopping is a
  small-L artifact (the k-curve is flat at L=1e4), and jackknife inflation is unnecessary because the
  plain sandwich is already calibrated (SE/SD≈1). The G2c failure was entirely the `refine="joint"`
  default's alternating drift.
- **For the *feasible* deliverable, the missing piece is the SE, not the point estimate or early
  stopping.** The anchor-orientation variance is real and large (SE/SD≈0.15) and is **not** captured
  by the plain sandwich or the split-document jackknife. The calibrated deliverable SE is a
  **full-chain document bootstrap** (resample documents → re-run anchors + orientation + frozen-Φ
  refinement → empirical SD, perm+sign aligned). That, plus `refine="frozen_phi"` (already the
  needed default over `joint`), is the A' object to re-certify.

## One-line bottom line

Alternating-drift → **A' viable** with inferential object = **frozen-Φ (in-regime: at convergence, no
early stop) + a full-chain document-bootstrap SE**; the plain sandwich is calibrated only at clean
orientation, and the split-document jackknife does **not** repair the feasible anchor-orientation
variance.
