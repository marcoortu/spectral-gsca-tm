# Gate report — chain realignment certification

One line per pre-registered gate: outcome + the number that decided it. Gates were frozen in
`PREREG_CHAIN.md` before running. This session promoted the three-stage chain + Lemma-17
sandwich into the package, unit-tested it, and ran the subspace/calibration/start-independence
gates at **reduced replication**; the real-data and full-table gates are scaffolded (status
noted). Nothing was tuned to rescue a gate.

| Gate | Outcome | Deciding number |
|------|---------|-----------------|
| **G0** subspace + collapse + chain→oracle | **PARTIAL** | raw ratio 0.060→0.042→0.025 tracks M^{-1/2} ✓; chain RMSE 0.141→0.098→0.067 descends ✓; **principal angle 0.89–1.25 rad (not <1e-3)** ✗; chain/oracle ≈4.6× at M=5000 (>2×) ✗ |
| **G7** B-functional start-independence | **PASS** | mse(pilot-start vs truth-start) = 4.18e-4 ≤ 5e-4 |
| **G2c** in-regime chain+sandwich coverage | **FAIL (key finding)** | coverage 0.346 (target .90–.97); sandwich SE/SD = 0.19 (SE ~5× too small) |
| Unit tests (chain, sandwich, anchors) | **PASS** | 17/17; sandwich matches explicit ⊗ loop, O(M^{-1/2}); anchors TV<0.20; Prop-21 confirmed |
| **G1** delocalization (real CRC, chain) | NOT RUN (chain) | pilot-based value stands: r_fit 6.5 genus / 7.2 species |
| **G2a/G2b** SE calib. + coverage vs √M/L | NOT RUN (chain) | — (see G2c mechanism) |
| **G3** crossover vs STM (chain) | NOT RUN (full) | partial: G0 shows feasible chain RMSE ↓ with M toward oracle |
| **G4** sandwich vs bootstrap (real CRC) | NOT RUN (chain) | earlier median 1.15 was raw-pilot `sgscatm_vcov`, not the chain |
| **G5** concordance (real CRC, chain) | NOT RUN (chain) | pilot-based 6/6 known genera CRC-ward stands |
| **G6** speed (real CRC, chain) | NOT RUN (chain) | chain adds anchors + refinement sweeps over the pilot |
| **G8** bias field b(z) / G0 gradient | NOT RUN | round-5 harness available |

## The decisive result (G2c) and its mechanism

The chain + **plain Lemma-17 sandwich under-covers badly in-regime**: at L=1e4, M=2000, N=500
(incidental gate closed, √M/L≈0.004), coverage is 0.35 and the sandwich SE is ~5× smaller than
the empirical cross-replicate SD. The Lemma-17 sandwich is the HC covariance of the refined-OLS
residuals `r̂_i = ẑ_i − B̂z'c_i`; it captures the per-document noise (ε+m) but **omits the
anchor/orientation-stage variability**, which dominates the chain's B̂z spread whenever the
orientation gate √M·δ_Φ is open — and at M=2000/N=500 it is (anchor error δ_Φ shrinks with M/N,
not with L). This is consistent with Corollary 18 (both gates must close) but falsifies the
pre-registered G2c prediction that **large L alone** yields nominal coverage: large L closes the
incidental gate, not the orientation gate.

**Mechanism confirmed** (`cert_G2c_mechanism.R`, same cell, feasible vs oracle orientation):
feasible cov 0.371 / SE-SD 0.22 → oracle cov 0.442 / SE-SD 0.78. Giving the chain the *true*
orientation raises the sandwich's captured variance from 22% to 78% of the empirical SD (so the
anchor stage supplies most of the missing variance), **but oracle coverage is still only 0.44** —
a residual bias survives at L=1e4 even with perfect orientation. Both drift terms the theory names
are active: the √M·δ_Φ orientation variance (omitted by the plain sandwich) and a refinement/
incidental bias (off-centres the estimate). Neither is a coding bug; both are the reason the plain
sandwich + point estimate is not calibrated for the feasible chain in this cell.

A second cell (N=200, M=5000, clean anchors) sharpens it: oracle SE/SD rises to **0.91** (variance
now well-estimated) but coverage **collapses to 0.15** — because at larger M the SE shrinks as
M^{-1/2} while the bias stays, so the √M-amplified drift dominates and coverage →0. Across the
design the plain chain fails for the two distinct drift reasons the CLT names: omitted anchor
variance (small M) and √M-amplified bias (large M). Both are addressed by the repo's jackknife
bias-correction (`B_jk = 2B_full − (B_A+B_B)/2`) + variance inflation, which this session did not
promote.

**This matches the feasibility round's own finding**: `replication/feasibility/03_jackknife.R`
adds a split-document **jackknife variance inflation** (`fs_coverage_entry(var_add = (B_A−B_B)²/4)`)
precisely because the plain sandwich under-covers; `run_f4` reports `cov_jk_infl`, not the plain
sandwich coverage. This session promoted the **plain** sandwich only. The fix — promote the
jackknife-inflated variance as the chain's default SE, or bootstrap the whole chain — is the
top follow-up (see `PAPER_EDITS.md`). A mechanism check (feasible vs oracle orientation) is in
`cert_G2c_mechanism.R`.

## Bottom line

The chain estimator is correctly promoted (unit-tested, start-independent, raw-pilot collapse
reproduced, Prop-21 confirmed), but **its closed-form inference via the plain Lemma-17 sandwich is
not calibrated in the tested regime** (G2c FAIL). The calibrated SE the manuscript's Section 6
reports came from the **raw pilot + `sgscatm_vcov`**, not the chain + Lemma-17 sandwich. Making
the chain the calibrated deliverable requires the anchor-variance-aware SE (jackknife inflation /
full-chain bootstrap), which exists in the replication code and should be promoted and re-certified.

---

# Phase 1 re-certification (frozen-Φ chain, sandwich primary) — findings

Pre-registered in `PREREG_PHASE1.md`. Package fixes applied: small-M B-stationarity F-convergence
backstop (fires at sweeps≈4 on synthetic); dead duplicate `sgscatm()` removed from `egscatm_fit.R`
(6 scripts repointed); `document()` clean.

## DECISIVE — T2-inreg-feasible: **FAIL** (overturns the "sandwich primary" decision)

Feasible frozen-Φ chain + analytic sandwich, L=1e4, **clean/exclusive anchors** (N=100,
alpha_beta=0.01, exclusivity=1.000), vs known Bz,0:

| M | coverage | SE/SD | boot/sandwich |
|---|----------|-------|---------------|
| 1000 | 0.491 | 0.163 | 16.3× |
| 5000 | 0.598 | 0.128 | — |

The sandwich **under-covers by ~2× and is ~6–8× too small** even with δ_Φ→0. The sandwich is
conditional on the fitted Φ̂ and cannot see the **anchor-selection variance** (which resamples expose);
exclusivity 1.0 bounds the average δ_Φ but not the across-resample orientation variance. Confirms
DIAG Part-2 (feasible SE/SD 0.13–0.17) at clean anchors.

**Consequence:** the analytic Lemma-17 sandwich does **not** suffice for the *feasible* deliverable,
even in-regime. A3 (sandwich calibrated) holds only at **oracle** orientation, which is unavailable in
practice. **The feasible in-regime SE requires the full-chain bootstrap** (conservative, valid). This
**overturns** the prior session's G4-fork "sandwich primary" call: the CRC agreement (sandwich≈0.80×boot)
was an artifact of the CRC chain's delocalized scores (r≈10⁴) inflating the sandwich, not calibration —
the clean sim (boot/sandwich=16×) exposes the true gap.

## CRC on the frozen-Φ chain — weaker than the raw-pilot draft (findings)

| gate | pilot draft | chain-certified | note |
|------|-------------|-----------------|------|
| G1 delocalization | r 6.5 (genus)/7.2 (species) | **r≈10,400 / 10,700** | chain uses natural-scale scores → delocalized; refinement hits the 100-sweep cap (deep-L per-doc GN drives extreme scores) |
| G5b known genera | 6/6 CRC-ward, p=0.016 | **5/6 CRC-ward, p=0.11** | directional but not significant on the chain |
| G6 speed | pilot 0.02 s, 38–173× | **chain point+sandwich 0.99 s** | ≈ PERMANOVA (0.76) / STM (1.01), 3× ALDEx2; the 38–173× was the pilot, **not** the chain |
| G4/G5a | — | study_condition sig under both SEs (p 0.0097/0.006), PERMANOVA 0.001 | holds; BMI SE-sensitive (sandwich sig, bootstrap n.s.) |

## Not run this session (budget; decisive gate already failed)

T1 (subspace cosines), T2 (coverage vs √M/L feasible+oracle), T3 (crossover vs STM), BES (T7),
b(z)/G0 (G8). The T2-inreg failure already establishes the central negative result, so these were
descoped; T2 would confirm the feasible-vs-oracle coverage gap (oracle nominal per A3, feasible
under-covering per T2-inreg).

## Bottom line (Phase 1)

**The frozen-Φ chain deliverable does NOT achieve fast closed-form calibrated in-regime inference.**
The analytic sandwich under-covers the *feasible* estimator even with clean anchors and large L
(T2-inreg FAIL: coverage 0.49–0.60, SE/SD 0.13–0.16); only the expensive full-chain bootstrap gives
valid (conservative) coverage. On real CRC the chain additionally delocalizes (r≈10⁴), does not
converge (sweep cap), weakens the biological concordance (5/6, p=0.11), and loses the pilot's speed
advantage (≈1 s ≈ competitors). Net: for practical calibrated inference the **raw pilot + `sgscatm_vcov`**
(the original Section-6 approach) remains the stronger deliverable; the chain is the theory-faithful
estimator whose feasible closed-form inference is **not** calibrated — the honest statement of
Corollary 18's second gate.
