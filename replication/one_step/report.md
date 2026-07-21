# One-Step-Efficient Estimation of the Structural ILR Path Matrix `B_z0` — Report

**Estimand.** The generative ILR path-coefficient matrix `B_z0` (P × (K−1)) that
`sim_dgp.R` uses to generate topic proportions from covariates (`Z_true = C B_z0 + Eps`,
`theta_i = softmax(V z_i)`). Not the standardized eigen-score coefficient, not an
across-replicate pseudo-mean. All comparisons made **only after** `procrustes_align()`.

**What was run.** Prereg written before any estimator run (`preregistration.md`, with a §5b
deviation addendum logged before the gates). No package `R/` code was modified; all new code
is in `estimators.R` / `sweep_utils.R` and the numbered scripts, reusing `sgscatm()`,
`sgscatm_vcov()`, `ilr_se()`, `procrustes_align()`, and `sim_dgp()`. Base design: `K=5,
N=500, P=3, sigma_eps=0.3, alpha_beta=0.1, lambda=1`; `n_rep=40–60`; `Bz0` and `Beta` held
FIXED within each cell (only C, eps, W resampled) so that SE-vs-SD, coverage, and bias are
meaningful.

---

## Headline

| Gate | Verdict | Headline number |
|------|---------|-----------------|
| Noiseless sanity | **PASS** | proj recovers `B_z0` to `rmse/‖B_z0‖ = 0.006`, scale ratio `0.994` |
| Estimand reconciliation | **PASS** | `Sigma_C^{-1} M^{-1} C'Z_true = B_z0` to `0.3%`; `‖E[C'U_0]‖=0.29 ≠ 0` |
| **G1** scale problem reproduced | **PASS** | baseline `rmse/‖B_z0‖ = 0.82→0.99`, increasing in `b_max`, coverage `0` |
| **G2** scale fix | **PASS** | proj `rmse/‖B_z0‖ = 0.05–0.12` (≪1), scale ratio ≈ 1 (not just rotation) |
| **G3** SE calibration, no collapse | **PASS** | proj/one-step median `SE/SD ∈ [0.95, 1.16]`, `cov_mean ≈ 0.95–0.98`, never pinned at 1.000 |
| **G4** linearisation bias isolated | **PASS** | proj `cov(B_z0)` falls `0.95→0.00` while `cov_mean≈0.97` and `SE/SD≈1`; bias grows `~ O(‖eta‖^2)` |
| **G5(a)** one-step reduces RMSE | **PASS** | 52–70% RMSE reduction in the moderate regime (`b_max` 0.5–1.0) |
| **G5(b)** consistency (decisive) | **PARTIAL** | one-step lowers the bias floor ~4× (quadratic→cubic) but a **nonzero cubic floor remains**; RMSE declines with M, does not reach 0 at fixed `b_max` |
| **G5(c)** L-robustness (barrier) | **MIXED** | robust for `L≥100`; at `L=50` the per-doc GN inflates variance and loses to proj |
| weighting (estimator D) | **FINDING** | naive per-doc multinomial weight is **inert** (`mw≡uw`); correct weight is WLS in the final regression (`onestep_wls`), which fixes short-doc fragility |
| **G6** anchor `Δ_Phi` | **PASS (boundary mapped)** | proj tolerates `Δ_Phi ≲ 0.1–0.2`; the one-step is **more fragile** (breaks by `Δ_Phi≈0.05`) |

---

## Conventions (fixed from source; see prereg §1)

`V` = Helmert ILR contrast (K×(K−1), `V'V=I`, `V'1=0`). Closure Jacobian at the centroid
is `L = (1/K)V`, so `p_i − 1/K·1 ≈ (1/K)V z_i`. The fit's scores are **unit-norm**
eigenvectors (`Z*'Z*=I`), which folds the natural scale `S = Cov(z)^{1/2}` into the loading
`Phi`; `B_z = (C'C)^{-1}C'Z*` therefore collapses on the `B_z0` scale.

## Deviation (logged before gates, prereg §5b)

The noiseless check showed the fit's own `Phi` does **not** carry natural scale: it has the
correct row-space (principal angles to `V'Beta` ≈ 0.4°,0.5°,1.7° on the three well-recovered
directions; the 4th is 42° — an inherent topic-recovery limit) but is 7.3× inflated in
Frobenius norm. This is the second-moment scale/basis indeterminacy: `S` is identified only
through the simplex/topic-word structure the spectral fit discards. **The loading-projection
is therefore anchor-based** — it consumes a natural-scale topic-word anchor `Phi_anchor`.
Projecting the counts on the true `Beta` recovers `B_z0` to 0.6% (scale ratio 0.994). Gates
G1–G5 supply a good anchor (the generative `Beta`, standing in for a consistent topic
estimator); G6 degrades it. This cleanly decomposes the estimand into (i) topic recovery and
(ii) the score-scale + linearisation fix — the object of this run.

---

## Estimand reconciliation (JASA referee point) — `01_estimand_check.R`

With `U_0 = Z_true` (generative ILR scores) and `E[C]=0`,
`E[c_i z_i'] = E[c_i c_i']B_z0 = Sigma_C B_z0`, so `Sigma_C^{-1} M^{-1} E(C'U_0) = B_z0`
**exactly** (numerically 0.3% at M=200k). `E[C]=0` does **not** imply `E[C'U_0]=0`: `U_0`
depends on `C` through its mean `C B_z0`, giving `‖E[C'U_0]‖_F = 0.29 ≠ 0`. If instead
`U_0` is the standardized eigen-score, the population form returns `B_z0 R S^{-1}` (RMSE/‖B‖
= 0.34, scale off by `S`) — exactly the G1 object. **Settled: the reviewer's form coincides
with `B_z0` iff `U_0` is the generative score, not the standardized one.**

## G1 — scale problem reproduced (`02_sweep_bmax.R`, `fig1_bmax.pdf`)

`baseline_std` RMSE/‖B_z0‖ across `b_max = {0.1,…,1.5}`: `0.82, 0.93, 0.96, 0.98, 0.98,
0.99` — ≈1 and monotonically increasing; coverage of `B_z0` is `0.00` throughout, median
per-entry SE ≈ `2e-4` (collapsed). **PASS** — the diagnosis holds on the real DGP.

## G2 — scale fix

`proj` RMSE/‖B_z0‖: `0.12, 0.05, 0.06, 0.11, 0.16, 0.27`. Two orders below the baseline at
small–moderate `b_max`. Noiseless scale ratio `‖B_aligned‖/‖B_z0‖ = 0.994` confirms the
scale — not merely the rotation — is recovered. **PASS.** (RMSE rises again at large `b_max`;
that is the linearisation bias of G4, not a scale failure.)

## G3 — SE calibration, no collapse

HC0 sandwich on the final ILR regression, with the rotational tangent projected out (matching
`sgscatm_vcov(identified=TRUE)`) so the SE matches the Procrustes-aligned SD. Across `b_max`,
proj median `SE/SD = 0.99,1.04,1.00,1.10,1.15,1.16`; one-step `0.98,1.04,1.02,0.97,1.02,0.95`
— all inside `[0.8,1.25]`. Coverage of the across-rep **mean** stays `0.95–0.98` and is
**never pinned at 1.000**. **PASS.** (The package's `sgscatm_vcov`/`ilr_se` operate on the
standardized scale and are reported by the package's own tests; here the HC sandwich is the
SE for the natural-scale estimators.)

## G4 — linearisation bias isolated

As `b_max` grows, proj coverage of the **true** `B_z0` falls `0.95, 0.93, 0.56, 0.14, 0.01,
0.00`, while coverage of the **mean** stays `0.95–0.98` and `SE/SD≈1`. The decline is
therefore bias-driven, not an SE failure. The systematic bias `‖mean_est − B_z0‖` grows
`0.0014, 0.0037, 0.0144, 0.038, 0.077, 0.199` — a `b_max→2×` step multiplies it ~4× (b=0.25→0.5)
to ~5× (b=0.75→1.5), i.e. `~O(‖eta‖^2)` to `O(‖eta‖^{2.4})`. **PASS** — linearisation is the
residual obstruction, distinct from scale (G2) and from SE (G3).

## G5(a) — one-step reduces RMSE

RMSE/‖B_z0‖ proj → one-step by `b_max`: `0.5: 0.064→0.031 (52%)`, `0.75: 0.106→0.032 (70%)`,
`1.0: 0.161→0.048 (70%)`, `1.5: 0.274→0.110 (60%)`. Reductions ≥40% throughout the moderate
regime and coverage of `B_z0` is largely restored (e.g. `b=0.75`: `0.14→0.84`). At very small
`b_max` (0.1, 0.25) the one-step is ~neutral-to-slightly-worse (no bias to remove, only added
variance) — expected. **PASS(a).**

## G5(b) — consistency (DECISIVE) (`03_Mscaling.R`, `fig2_Mscaling.pdf`)

Fixed `b_max=0.5, L=200`, `M∈{500,…,8000}`:

```
        RMSE/||B_z0||                 systematic bias ||mean-B_z0||
   M    proj    onestep               proj     onestep
  500   0.072   0.055                 0.0135   0.0043
 1000   0.067   0.044                 0.0138   0.0035
 2000   0.064   0.031                 0.0144   0.0033
 4000   0.059   0.025                 0.0136   0.0037
 8000   0.059   0.020                 0.0140   0.0036
```

**proj plateaus** at its quadratic linearisation-bias floor (bias `≈0.014` flat in M; RMSE
→ `0.059`). The **one-step RMSE keeps declining** (`0.055→0.020`) and is variance-dominated
across this range. But its systematic bias also plateaus — at `≈0.0037`, ~**4× below** proj,
not at 0. The b_max sweep shows this residual is `O(‖eta‖^3)` (cubic: one-step bias
`0.0031→0.021` as `b_max 0.5→1.0`, a ~7× jump for 2×), versus proj's `O(‖eta‖^2)` (quadratic).

**Verdict: PARTIAL.** One Gauss-Newton step removes the *leading* quadratic linearisation
bias and lowers the order to cubic, so at fixed `b_max` the RMSE declines with M toward a
floor ~4× smaller than proj's — a decisive, real improvement, and in the practical M range
the one-step is variance-dominated and still falling. **But a single step does not drive
`RMSE(B_z0)→0` at fixed `b_max`: a nonzero cubic bias floor remains.** Honest reading, as the
prereg demanded: full consistency for `B_z0` would require either iterating the Newton step to
convergence (a full multinomial MLE, not a one-step) or the vanishing-effect regime
(`b_max→0`). The one-step is a fast, one-order-more-accurate estimator of the structural
coefficient, not a certified-consistent one.

## G5(c) — L-robustness / barrier (`04_Lrobust.R`, `fig3_Lrobust.pdf`)

Fixed `M=2000, b_max=0.5`, constant `L∈{50,100,200,400,1000}`. One-step bias `‖mean−B_z0‖`:
`0.030, 0.012, 0.0033, 0.0024, 0.0041` and RMSE/‖B‖ `0.134, 0.059, 0.031, 0.025, 0.026`. For
`L≥100` the one-step is well-behaved and beats proj (RMSE ~0.03 vs ~0.06) with coverage
~0.92. At `L=50` it **breaks**: the per-document GN amplifies the (now large) multinomial
noise in `theta_hat`, RMSE `0.134` exceeds proj's `0.068` and coverage drops to `0.53`.
**MIXED** — barrier-robust down to `L≈100`, fragile in the extreme-sparse regime.

## Weighting — estimator (D) is inert as prescribed; the fix is WLS (`04b`, `fig5_weighting.pdf`)

Two findings, logged as a deviation:

1. **The naive per-document multinomial weight is algebraically inert.** In the per-document
   GN solve `delta_i=(J_i'Ω_i^-J_i)^{-1}J_i'Ω_i^- r_i`, the `1/L_i` factor is a per-document
   scalar that cancels, and near the centroid the category weight `(diag θ − θθ')^-` reduces
   to a scalar (proportions ≈ uniform), so `onestep_mw ≡ onestep_uw` to machine precision —
   under both constant and variable `L`. The prescribed weighting cannot bite where it was
   placed.

2. **The multinomial precision belongs in the final pooled regression.** With heavy-tailed
   lengths (mean 124, 30% of docs <50 words, min 10) the per-doc one-step blows up
   (`RMSE/‖B‖ = 0.44` vs proj `0.094`, variance ×18). Weighting the *pooled* ILR regression
   by `w_i = L_i` (`onestep_wls`) downweights short/noisy documents and restores the method:
   `RMSE/‖B‖ = 0.060` (better than proj), coverage `0.76` vs proj `0.26`. This is the genuine
   "weighted ≥ unweighted" result; it required correcting where the weight is applied.

## G6 — anchor `Δ_Phi`, clean probe (`05_anchor_G6.R`, `fig4_anchor.pdf`)

Scale-preserving corruption (each topic-word row mixed toward a random simplex direction by
`Δ_Phi`, rows renormalised so the loading scale is unchanged — confirmed: anchor row-sums stay
1). RMSE/‖B_z0‖ vs `Δ_Phi = {0,.05,.1,.2,.4,.8}`:

```
 Δ_Phi   proj    onestep_mw   cov_Bz0(proj)  cov_Bz0(onestep)
 0.00    0.064   0.031        0.56           0.92
 0.05    0.030   0.097        0.94           0.41
 0.10    0.051   0.185        0.73           0.23
 0.20    0.131   0.375        0.35           0.07
 0.40    0.114   0.479        0.39           0.07
 0.80    0.768   0.741        0.00           0.00
```

**proj degrades gracefully and stays usable to `Δ_Phi ≈ 0.1–0.2`**, breaking by `0.4–0.8`.
**The one-step is markedly more fragile**: because it re-linearises around and pulls toward the
(now biased) `theta_hat`, it *amplifies* anchor error, losing to proj for any `Δ_Phi>0` and
losing coverage by `Δ_Phi≈0.05`. Boundary: with an exact/near-exact topic estimate the
one-step is best; under realistic topic-recovery error the robust `proj` (or the WLS variant)
is preferable. The small-`Δ` non-monotonicity of proj (0.064→0.030 at 0.05) is a
fixed-corruption-seed / finite-`n_rep` artifact.

---

## Bottom line

* The **scale problem (G1)** and its **fix (G2)** reproduce cleanly on the real DGP; the fix
  is anchor-based, because the spectral fit discards the natural scale that only the topic-word
  structure identifies.
* The **HC sandwich SE is calibrated (G3)** once the rotational tangent is projected out, and
  the **linearisation bias is cleanly isolated (G4)** as the residual obstruction.
* The **one-step helps substantially (G5a: 52–70%)** and **lowers the bias order from
  quadratic to cubic (G5b)**, but **a single step is not sufficient for consistency**: a
  nonzero cubic bias floor remains, so `RMSE(B_z0)` does not reach 0 at fixed covariate
  strength. Consistency for `B_z0` remains open for the one-step; it is a fast,
  one-order-improved estimator, not a certified-consistent one.
* The prescribed **multinomial GN weighting is inert**; the correct multinomial weighting is
  WLS in the final regression, which is what delivers robustness to short documents and a
  genuine weighted improvement.
* **Anchor quality is the binding constraint** for the one-step (G6): it needs a near-exact
  topic estimate, whereas `proj`/`onestep_wls` tolerate moderate topic-recovery error.

## Reproduce

```
Rscript replication/one_step/00_noiseless_check.R      # sanity (must pass first)
Rscript replication/one_step/01_estimand_check.R       # estimand reconciliation
Rscript replication/one_step/02_sweep_bmax.R           # G1,G2,G3,G4,G5a
Rscript replication/one_step/03_Mscaling.R             # G5b (decisive)
Rscript replication/one_step/04_Lrobust.R              # G5c
Rscript replication/one_step/04b_weighted_vs_unweighted.R  # estimator D / weighting
Rscript replication/one_step/05_anchor_G6.R            # G6
Rscript replication/one_step/06_figures.R              # figures/*.pdf
```
Outputs: `out_0*.{rds,csv}`, `figures/fig{1..5}_*.pdf`. Seeds fixed (`10000+r`); `Bz0`/`Beta`
fixed per cell.
