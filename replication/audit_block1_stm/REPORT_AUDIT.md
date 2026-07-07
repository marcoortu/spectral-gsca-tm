# Block 1 forensic audit + fair STM comparison

Follow-up to `replication/basin_check/REPORT.md` (§1 alignment audit).
All numbers reproducible via

```
Rscript replication/audit_block1_stm/01_audit_block1.R    # A2 + A3  (~25 s)
Rscript replication/audit_block1_stm/02_stm_fair.R        # B1 + B2  (~95 s)
Rscript replication/audit_block1_stm/03_regime_map.R      # Task C   (~10 s)
Rscript replication/audit_block1_stm/04_report.R
```

## Executive summary

- **H1 (RMSE plateau): confirmed.** The published Block 1 RMSE
  (0.2387 / 0.2450 / 0.2495 at M = 500 / 1000 / 2000) is reproduced to
  four decimals on the same seeds and sits at **92 / 94 / 96 %** of the
  zero-estimator level `sqrt(mean(Bz0^2)) = 0.2598` — mildly
  *increasing* in M, because `‖fit$Bz‖/‖Bz0‖` shrinks like 1/√M
  (0.086 → 0.042). The published curve cannot scale as M^{-1/2}
  because it measures the distance of an (essentially) zero matrix
  from `Bz0`.
- **H2 (coverage 1.000): confirmed, mechanism identified.** The
  coverage check compares a CI centred on the scale-collapsed aligned
  estimate (offset from truth ≈ 0.21–0.22 per entry) with half-widths
  of **47.8 / 26.0 / 11.7** (median): 50–220× the distance under test.
  Two compounding code-level causes, quoted in §A1: `.rotate_se()`
  omits the 1/√M factor (its comment asserts "/M is already in
  ilr_se_analytical", but `ilr_se_analytical` applies /M only to its
  own `se` slot, not to the returned `Sigma_Bz`), and `Sigma_Bz`
  itself is O(1)–O(10²) because it is computed for the unit-norm-score
  convention on a spectrum with near-degenerate gaps. Even the
  /√M-corrected counterfactual still gives 1.000 / 0.988 / 0.650 —
  both errors matter.
- **A3 (corrected two-step): order-of-magnitude accuracy gain, but both
  pre-registered A3 predictions FAIL, with one diagnosed cause.** RMSE
  drops from ≈ 0.24 (published pipeline) to 0.036 / 0.030 / 0.023, but
  the log-log slope is −0.295 (95 % CI [−0.355, −0.235], excludes
  −0.5) and entrywise coverage is 0.867 / 0.765 / 0.712 (below the
  registered [0.90, 0.98]). Decomposition: **bias² is 71–78 % of the
  MSE at every M** and the mean |bias|/SE grows 0.91 → 1.26, while the
  √variance component alone scales with slope −0.409. This is the
  finite-L incidental-parameter bias floor at L = 200 (basin_check E5
  measured the same floor from *truth-started* refinement, so it is
  not a start artifact). The planned two-step CLT needs
  an explicit small-bias condition (√M/L → 0; here √M/L ≈ 0.22 at
  M = 2000) or a per-document bias correction.
- **B (fair STM): the coordinate-mismatch hypothesis is refuted in
  effect size, and the published STM column stands.** The old worker
  never touched STM's `gamma` — it regressed ILR-mapped `theta` on C,
  which is coordinate-consistent. Running both extractions on paired
  data: STM-native-gamma (ALR→ILR-mapped, unit-tested to 1e-15) and
  the old theta-based functional give *identical* MSE (weak 0.0084 vs
  0.0085; strong 0.0220 vs 0.0220), matching published Table 3
  (0.009 / 0.021). **Refined ≤ STM everywhere**: weak 0.0002 vs 0.0084
  (42×), strong 0.0021 vs 0.0220 (10×), M = 5000 strong 0.0013 vs
  0.0136 (10×). New diagnostic: STM's coefficient norm is *inflated*
  1.5–2.1×, and in the weak regime STM's MSE (0.0084) is slightly
  *above* the zero-estimator level (0.0075). One real defect found in
  the old worker: its `Bz0_r` was drawn from an unseeded RNG in a
  fresh subprocess, so Table 3's sgscatm and STM columns were computed
  on **different data with different truths** (distribution-level
  comparison only, and non-reproducible).
- **C (regime map): no optimizer breakdown anywhere in
  b_max ∈ [0.15, 1.5]** (0 Armijo failures, Levenberg ν never left its
  1e-6 floor, 40/40 monotone). The pilot's oracle-subspace quality is
  flat to b_max = 1.0 (mse 0.0002) and degrades ~10× at b_max = 1.5
  (0.0024) as saturated documents appear (0.4 % at the truth). The
  refined estimator's *relative* error grows only mildly (4 % → 7 % of
  signal size); degradation near the simplex boundary is statistical
  (bias), not algorithmic.

---

## A1. Trace of the published Block 1 pipeline

**Bz0.** Fixed matrix, `run_simulation.R` lines 88–93:

```r
set.seed(2026)
Bz0_TRUE <- matrix(c(
   0.40, -0.20,  0.10,  0.30,   # covariate 1 → 4 ILR components
  -0.15,  0.35, -0.25,  0.05,   # covariate 2
   0.20,  0.10,  0.40, -0.30    # covariate 3
), nrow = P_COV, ncol = K_TOPICS - 1L, byrow = TRUE)
```

(The `set.seed(2026)` is vestigial — the entries are hard-coded.)
`sqrt(mean(Bz0^2)) = 0.2598`; true squared row norms
`(0.30, 0.21, 0.30)`. Replicate data: line 125–128,
`sim_dgp(..., Bz0 = Bz0_TRUE, seed = 10000L * which(M_VALUES_B1 == M) + r)`.

**RMSE.** Lines 132–142: `fit <- sgscatm(dat$W, dat$C, K, lambda = 1,
rotate = TRUE)`; `pa <- procrustes_align(fit$Bz, dat$Bz0)`; per-replicate
`mse = pa$mse = mean((fit$Bz %*% R - Bz0)^2)` with R **orthogonal only**
(`sim_utils.R` 30–36: `svd(crossprod(Bz_hat, Bz0))`, `R = u v'`, no
scaling). Table value: `rmse = sqrt(mean(x$mse[ok_mse]))` (line 178).
So the RMSE is measured between the raw-scale truth and the
unit-norm-score-convention estimate (basin_check §1: `fit$Z` has
`Z'Z = I`, `fit$Bz = (C'C)^{-1}C'Z`), which no orthogonal map can
rescale.

**Coverage.** Lines 140 and 152–158 — the analytical (eigenvector
perturbation) SE, **not** the bootstrap:

```r
se_res <- tryCatch(ilr_se_analytical(fit), error = function(e) NULL)
...
se_rot <- .rotate_se(res$se_res, res$R, P_COV, K_TOPICS - 1L)
z_q    <- qnorm(1 - (1 - CONF) / 2)
ci_lo  <- res$Bz_al - z_q * se_rot
ci_hi  <- res$Bz_al + z_q * se_rot
covers <- (dat$Bz0 >= ci_lo) & (dat$Bz0 <= ci_hi)
```

The entrywise CI is centred at the **Procrustes-aligned estimate**
(`Bz_al`, norm ≈ 0.04–0.09 of the truth's) and compared against the
**raw `Bz0`** (no transformation of the truth). The SEs come from
`.rotate_se` (`sim_utils.R` 167–176):

```r
RkI <- kronecker(t(R), diag(P))
Sigma_aligned <- RkI %*% se_res$Sigma_Bz %*% t(RkI)
se_vec <- sqrt(pmax(diag(Sigma_aligned), 0) / nrow(se_res$Sigma_Bz))
# Note: the /M factor is already in ilr_se_analytical
# Recompute from diagonal directly
se_vec <- sqrt(pmax(diag(Sigma_aligned), 0))
```

The first `se_vec` (divided by `nrow(Sigma) = P(K-1) = 12`, itself not
M) is **overwritten** by `sqrt(diag(Sigma_aligned))` with no
normalisation, on the stated belief that /M is already inside. But
`ilr_se_analytical` (`ilr_se.R` 179–182) divides by M only in its own
`se` slot and returns `Sigma_Bz` raw:

```r
se_vec  <- sqrt(pmax(diag(Sigma_Bz), 0) / M)
list(se = matrix(se_vec, P, Kp), Sigma_Bz = Sigma_Bz)
```

So Block 1's CIs use SEs that are √M times the intended convention, on
top of a `Sigma_Bz` whose scale already reflects the unit-norm-score
parametrisation and near-degenerate eigengap denominators
(`1/(d_l d_lp)` terms, `ilr_se.R` 155–167; the eigengap warning at
`tol_gap = 0.05` fires routinely in this design).

**Bootstrap (`ilr_se`, not used by Block 1).** For completeness: it
resamples documents (rows of W and C jointly) with replacement and
refits `sgscatm` per replicate (`ilr_se.R` 51–63), taking entrywise
SDs of the **unaligned** `Bz` across replicates — bootstrap draws are
not rotated to a common orientation, so sign/rotation flips inflate
those SDs; any future use should align each draw first.

## A2. Mechanical reproduction (20 reps per M, published seeds r = 1..20)

| M | published RMSE (50 reps) | reproduced RMSE | zero level | RMSE/zero | ‖fit$Bz‖/‖Bz0‖ | coverage (published SE) | coverage (SE/√M) | med CI half-width | med \|center−truth\| |
|---|---|---|---|---|---|---|---|---|---|
| 500  | 0.2387 | 0.2387 | 0.2598 | 0.92 | 0.086 | 1.000 | 1.000 | 47.80 | 0.214 |
| 1000 | 0.2450 | 0.2449 | 0.2598 | 0.94 | 0.061 | 1.000 | 0.988 | 26.01 | 0.217 |
| 2000 | 0.2495 | 0.2496 | 0.2598 | 0.96 | 0.042 | 1.000 | 0.650 | 11.74 | 0.220 |

- **H1 verdict: PASS.** Within 4–8 % of the zero-estimator level
  (registered band ~15 %), flat-to-increasing in M; the norm ratio
  falls as ~1/√M (the unit-norm score convention), which is exactly
  why the curve *approaches* the zero level from below instead of
  decaying as M^{-1/2}.
- **H2 verdict: PASS — mechanism, one referee-checkable sentence:**
  *the CI half-widths (median 11.7–47.8) are 50–220 times the distance
  they are asked to cover (≈ 0.22), because `.rotate_se` returns
  `sqrt(diag(Sigma_aligned))` without the 1/√M factor that
  `ilr_se_analytical` applies only to its own `se` slot, and `Sigma_Bz`
  itself is on the unit-norm-score scale — so every entrywise CI
  covers by construction.* The counterfactual column shows the missing
  √M is not the whole story: with SEs/√M, coverage is still 1.000 and
  0.988 at M = 500 and 1000 (half-widths ≈ 2.1 and 0.82 vs distance
  0.21) and collapses to 0.650 at M = 2000 — a mis-centred CI with a
  mis-scaled covariance never yields calibrated coverage; it merely
  flips from always-cover to under-cover as M grows.

## A3. Corrected diagnostics: two-step estimator (pilot → GL → k = 5 GN sweeps → OLS)

**Estimator start — an important negative finding.** The GL alignment
of the pilot uses `Z_true` (oracle), inherited from basin_check.
Feasible truth-free starts were tested and **fail**: starting the
identical 5-sweep refinement from the raw pilot `(fit$Z, fit$Phi)`
reaches an *equal or lower* objective value (F = 4.7310 vs 4.7311 at
30 sweeps) with mse_Bz = **0.061** (zero-estimator level) instead of
0.001; a simplex-projected-Φ start (44 % of pilot `Phi` entries are
negative) gives mse 0.03–0.19. The exact λ = 0 criterion's near-flat
floor contains points with arbitrarily bad B: **the objective does not
identify the orientation/scale of the score space — the aligned start
supplies it.** Any feasible version of the two-step theorem must add
an identification device (anchor conditions, or a pilot that outputs
an oriented basis); with the oracle alignment these diagnostics
prototype the *variance* of the theorem, as intended by the brief.

**Sandwich unit test** (heteroscedastic, row-correlated noise, M = 200,
P = 2, K−1 = 3, 3000 sims): mean sandwich vs empirical covariance of
vec(B̂), relative Frobenius error **0.040**. PASSED (gate 0.25).

| M | RMSE (two-step) | RMSE (pilot, published conv.) | norm ratio | entrywise coverage | row-norm coverage | Armijo fails | monotone |
|---|---|---|---|---|---|---|---|
| 500  | 0.0358 | 0.2386 | 1.103 | 0.867 | 0.700 | 0 | 50/50 |
| 1000 | 0.0305 | 0.2450 | 1.093 | 0.765 | 0.427 | 0 | 50/50 |
| 2000 | 0.0234 | 0.2495 | 1.070 | 0.712 | 0.373 | 0 | 50/50 |

Log-log slope of per-replicate RMSE on M: **−0.295**, 95 % CI
[−0.355, −0.235]. Figures: `results/rmse_vs_M.png`,
`results/coverage.png`, `results/qq_M2000.png`.

**Bias/variance decomposition (aligned scale):**

| M | bias² | variance | bias²/mse | mean \|bias\|/SE |
|---|---|---|---|---|
| 500  | 9.20e-04 | 3.71e-04 | 0.71 | 0.91 |
| 1000 | 7.41e-04 | 1.93e-04 | 0.79 | 1.16 |
| 2000 | 4.32e-04 | 1.20e-04 | 0.78 | 1.26 |

- **Registered predictions: FAIL, both.** The slope CI excludes −0.5,
  and entrywise coverage is below [0.90, 0.98] and *degrades* with M;
  row-norm coverage is worse (0.70 / 0.43 / 0.37; the largest-norm row
  drops to 0.10 at M = 2000, as squared norms aggregate the bias).
- **Diagnosis (single cause, demonstrated):** the √variance component
  alone scales with slope −0.409 ≈ −0.5, while bias² barely moves —
  the L = 200 incidental-parameter bias floor (each ẑ_i is estimated
  from one 200-word document through a nonlinear softmax; the bias
  does not average out in B). basin_check E5 measured the same floor
  from *truth-started* refinement (mse ≈ 4e-4 at L = 200 falling to
  ≈ 1e-4 at L = 1000), so it is a property of the criterion, not of
  the oracle start or k = 5 truncation. The QQ plot at M = 2000 shows
  the corresponding over-dispersion of pooled standardised errors
  (spread ≈ 1.4× normal, from entry-specific bias offsets of ±1.3 SE),
  not distributional pathology.
- **Consequence for the theorem:** the two-step CLT must either assume
  √M/L → 0 (here √M/L = 0.11 / 0.16 / 0.22 — visibly not small) or
  include a per-document bias correction (analytic second-order, or
  split-document jackknife). As registered, we report and do not
  correct. Independent corroboration: the earlier corrected rerun in
  `sim_rerun/tables/block1.csv` (different estimator convention) shows
  the same signature — coverage drifting 0.944 → 0.873 as M grows
  500 → 4000.

## B1. ALR → ILR map, unit test, and audit of the old worker

**Derivation.** STM: θ = softmax([η; 0]), reference topic K, so
η = alr(θ). With L = log θ: alr = L_{1..K−1} − L_K and ilr = V′L.
Since V′1 = 0, adding L_K·1 changes nothing:
**ilr = V′[alr; 0] = M_map·alr, M_map = t(V[1:(K−1), ])**. For score
matrices (rows = documents): `ILR = ALR %*% t(M_map)`; for prevalence
coefficients (STM's `gamma` is (P+1) × (K−1) with η_row = x′Γ):
**`Gamma_ilr = Gamma_alr[-1, ] %*% t(M_map)`** (intercept row dropped;
C is centred). Inverse: alr = A_map·ilr with
A_map = V[1:(K−1),] − 1·V[K,].

**Unit test (hard gate, PASSED):** on 200 synthetic compositions with
the package's own `V = ilr_contrast(5)`: max errors
ilr-from-alr 2.7e-15, alr-from-ilr 1.9e-15, `M_map %*% A_map − I`
1.1e-16, DGP round-trip `Bz0 → Gamma_alr → Bz0` 5.0e-16 (tol 1e-12).

**What the old worker actually did** (`block3_stm_worker.R` 73–89): it
**never used STM's `gamma`**. It took `fit_stm$theta`, mapped it to ILR
(`log(pmax(theta, 1e-10)) %*% V_ilr`), regressed on C, and
Procrustes-aligned — i.e. the same functional our pipeline uses, in
ILR coordinates throughout. **No ALR/ILR mismatch existed**, and the
paired rerun confirms it numerically: the theta-based and
gamma-native-mapped extractions agree to the third decimal in MSE
(below). The published Table 3 STM column was *not* computed in
mismatched coordinates.

**A real defect found instead:** worker lines 41–42 draw `Bz0_r` with
`runif` from the **unseeded** ambient RNG of a fresh subprocess, while
the sgscatm loop in `run_simulation.R` draws its own `Bz0_r` from the
main process's RNG stream. Same `sim_dgp` seeds, different `Bz0` ⇒
different Z, W. Table 3's sgscatm and STM columns therefore compare
**different datasets with different truths** (valid only as a
distribution-level comparison, and not reproducible run-to-run). The
rerun below is exactly paired instead.

## B3. Fair STM comparison (paired with basin_check E2 data; means (sd))

| regime | method | mse_Bz paper (Procrustes) | norm ratio ‖B̂‖/‖Bz0‖ | GL/oracle mse | time (s) | EM its |
|---|---|---|---|---|---|---|
| weak | pilot (published pipeline) | 0.0050 (0.0015) | 0.193 | 0.0001 (0.0001) | 0.9 | — |
| weak | pilot + refined (k=5, λ=0) | **0.0002 (0.0001)** | 0.891 | 0.0003 (0.0001) | 1.4 | — |
| weak | STM (native γ, ALR→ILR) | 0.0084 (0.0034) | 2.086 | — | 6.0 | 23 |
| weak | STM (θ→ILR→OLS, old worker) | 0.0085 (0.0035) | 2.090 | 0.0002 (0.0001) | 6.0 | 23 |
| strong | pilot (published pipeline) | 0.0768 (0.0177) | 0.055 | 0.0001 (0.0001) | 0.7 | — |
| strong | pilot + refined (k=5, λ=0) | **0.0021 (0.0012)** | 1.135 | 0.0024 (0.0012) | 1.1 | — |
| strong | STM (native γ, ALR→ILR) | 0.0220 (0.0153) | 1.497 | — | 7.6 | 34 |
| strong | STM (θ→ILR→OLS, old worker) | 0.0220 (0.0153) | 1.498 | 0.0002 (0.0001) | 7.6 | 34 |
| strong, M=5000 | pilot | 0.0753 (0.0186) | 0.024 | 0.0002 (0.0001) | 0.7 | — |
| strong, M=5000 | pilot + refined | **0.0013 (0.0010)** | 1.088 | 0.0014 (0.0011) | 2.8 | — |
| strong, M=5000 | STM (native γ) | 0.0136 (0.0042) | 1.411 | — | 33.7 | 28 |
| strong, M=5000 | STM (θ→ILR→OLS) | 0.0136 (0.0042) | 1.412 | 0.0000 (0.0000) | 33.7 | 28 |

Published Table 3 references: weak 0.009, strong 0.021 (M = 1000);
zero-estimator levels: weak b²/3 = 0.0075, strong 0.0833.

Findings and registered-prediction verdicts (B):

- **"Coordinate-corrected STM differs from the published column" —
  FAILS.** Both extractions coincide (γ-native vs θ-based identical to
  0.0001), and both reproduce the published values on paired data.
  The published STM column survives the coordinate audit.
- **"Refined ≤ STM in both regimes" — PASSES**, by 42× (weak) and 10×
  (strong; 10× again at M = 5000). The published "strong-regime gap"
  (pilot 0.077 vs STM 0.021) does not reappear: it inverts
  (refined 0.0021 vs STM 0.0220).
- **Norm-ratio diagnostic — as predicted for pilot and refined, new
  finding for STM:** pilot 0.19 / 0.055 / 0.024 (scale collapse,
  ~1/√M); refined 0.89–1.14 (≈ 1 ✓); **STM 1.4–2.1 — systematically
  inflated**, and in the weak regime STM's MSE (0.0084) sits *above*
  the zero-estimator level (0.0075): at this signal strength STM's
  prevalence coefficients are mostly noise scaled ~2× too large. Its
  oracle-GL column (0.0002) shows the *subspace* of its θ estimates is
  fine — like the pilot, STM suffers a scale/orientation problem under
  this metric, just in the opposite direction (inflation vs collapse).
- STM cost: 6–8 s per fit at M = 1000 (23–34 EM iterations), 34 s at
  M = 5000 — vs 1.1–2.8 s total for pilot + refined.

**Why a permutation-only (label-switching) metric is not applicable.**
A permutation metric presumes each method identifies topics up to
relabelling, so that the only nuisance is which label goes where. None
of the three estimators is identified that finely without anchor-type
conditions: the spectral pilot returns scores up to a general-linear
transform of the ILR factor space (varimax fixes a rotation
representative, but varimax has no population target for Gaussian
scores — the Gaussian factor model is rotationally invariant, so
simple-structure rotation is not consistent orientation recovery); the
refined estimator inherits whatever orientation its start supplies
(§A3); and STM's Γ lives in ALR coordinates whose reference topic is
itself an arbitrary label, so permutations act on Γ through a
non-orthogonal change of basis rather than a column shuffle. A
permutation-only alignment would therefore be unattainable for every
method and would measure label conventions, not estimation error.
Orthogonal Procrustes (paper metric) and oracle GL bracket the
meaningful range.

## Task C. Operating-regime map (M = 1000, 10 reps/cell, exploratory)

| b_max | sat(θ_true) | sat(θ_end) | pilot GL mse (oracle) | pilot paper mse | refined paper mse | Armijo fails | max ν | monotone |
|---|---|---|---|---|---|---|---|---|
| 0.15 | 0.000 | 0.000 | 0.0001 (0.0000) | 0.0052 (0.0011) | 0.0003 (0.0002) | 0 | 1e-06 | 10/10 |
| 0.50 | 0.000 | 0.000 | 0.0001 (0.0001) | 0.0906 (0.0242) | 0.0036 (0.0015) | 0 | 1e-06 | 10/10 |
| 1.00 | 0.000 | 0.001 | 0.0002 (0.0001) | 0.2806 (0.0699) | 0.0218 (0.0093) | 0 | 1e-06 | 10/10 |
| 1.50 | 0.004 | 0.016 | 0.0024 (0.0019) | 0.7327 (0.1996) | 0.0523 (0.0259) | 0 | 1e-06 | 10/10 |

(sat = share of documents with max_k θ_ik > 0.99.) The per-document
Newton machinery never struggles in this range: no Armijo failures, ν
pinned at its floor, all runs monotone. The pilot subspace (oracle GL
column) is essentially perfect through b_max = 1.0 and degrades ~10×
at b_max = 1.5, where boundary saturation appears. Refined error
relative to signal size (mse ÷ b²/3) grows only 4 % → 7 %: approaching
the simplex boundary raises the statistical difficulty (larger |z|,
stronger finite-L bias), not the algorithmic one. Figure:
`results/regime_map.png`.

## Deviations from the brief

1. **A3 start uses the oracle GL alignment** (basin_check machinery, as
   the brief's estimator definition implies). Feasible truth-free
   starts (raw pilot; simplex-projected Φ) were tried and demonstrably
   fail (§A3) — kept as a headline negative finding rather than hidden
   in the estimator definition.
2. A2 uses 20 reps per M (brief: "20 is enough"), taken as the first
   20 seeds of the published 50-seed sequence, which is why the
   reproduction matches to 4 decimals.
3. `01_audit_block1.R` was extended (before the final run, same seeds)
   to store per-replicate `B_tilde` and SE matrices for the
   bias/variance decomposition.
4. Task C tracks the Levenberg ν by orchestrating `bc_z_step` /
   `bc_phi_step` directly (5 sweeps); `bc_refine` does not expose ν.
   No other duplication of basin_check machinery.
5. STM fit in-session (sanctioned; the dual-R subprocess was
   unnecessary — `stm` loads natively under R 4.5.1; note it was built
   under R 4.5.2, warning only). Optional M = 5000 strong block was
   run (budget held).
6. Windows platform: PSOCK workers (10) instead of `mclapply`, as in
   basin_check.

## Runtime and seeds

- Wall-clock (10 PSOCK workers, 22-core machine): sandwich + map unit
  tests < 30 s; A2 4.1 s; A3 15.8 s; B2 39.4 s; B2+ (M = 5000) 53.1 s;
  Task C 9.4 s; report ≈ 5 s. **Total ≈ 2.5 min** (budget 1.5 h);
  development-time side experiments (feasible starts) ≈ 3 min extra.
- Seeds: A2 replicates reuse the published Block 1 scheme
  `10000·M_index + rep` (M_index 1, 2, 3 for 500, 1000, 2000; rep
  1..20). A3: `60000 + M_index·1000 + rep`, rep 1..50. B2:
  basin_check's `90000 + regime_index·1000 + rep` (weak = 1,
  strong = 2), rep 1..20, so the pilot/refined columns are computed on
  the identical datasets as basin_check E2 (M = 5000 uses the same
  scheme at rep 1..10). Task C: `40000 + bmax_index·1000 + rep`, rep
  1..10. Unit tests: 7001 (sandwich), 7101 (map).
- Files: `01_audit_block1.R`, `02_stm_fair.R`, `03_regime_map.R`,
  `04_report.R`; `results/*.rds` (raw per-replicate outputs),
  `results/summary_audit.csv`, `results/summary_stm.csv`,
  `results/tables_audit.md`, figures `rmse_vs_M.png`, `coverage.png`,
  `qq_M2000.png`, `regime_map.png`, logs `run01-03.log`.
- No package function or pre-existing script was modified.
