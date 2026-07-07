# `sim_rerun/` — corrected SE, clean-DGP re-run & paper artefacts

Validates the **corrected analytical standard error** for the ILR path
coefficients (new package function `sgscatm_vcov()`), re-runs Block 1 on a
clean design that removes the `P<K-1` near-degeneracy, and produces the
lambda-scaling diagnostic. Run from the package root with R 4.5.1:

```
Rscript sim_rerun/run_block.R            # Block1, lambda diag, Blocks 2-5, verdict
Rscript sim_rerun/run_block.R 1          # Block 1 only
Rscript sim_rerun/run_block.R 0          # lambda diagnostic only
SGSCATM_REDUCED=1 Rscript sim_rerun/run_block.R
```

## Package change (Part 1) — `R/vcov.R`
`sgscatm_vcov(fit, W, C, rotation = NULL, identified = TRUE)` returns the
asymptotic covariance and SEs of `vec(B_z)` on the standardized scale. The
influence function has **three terms** — (A) regression, (B) eigenvector
fluctuation, (C) **eigenvalue fluctuation (new)** — all on the O(1)
eigenvalue scale `rho_k = eig((1/M) H'H)`, with the previously missing final
`/M` so `SE = O(M^{-1/2})`. `lambda` does not enter at first order (the
covariate block perturbs the leading eigenspace at `O(lambda/M)`).

Because the spectral estimand is identified only up to rotation, when the
estimate is aligned to a reference the aligned estimator has no variance
along the rotational tangent `{vec(B Ω): Ω' = -Ω}`; `identified = TRUE`
projects that tangent out of the covariance. **This projection is required
for the SE to match the empirical SD / jackknife** (without it the raw
covariance overstates the SE ~5x on this design). `sgscatm(..., se = TRUE)`
attaches `$B_z_se`. Unit tests: `tests/testthat/test-vcov.R`.

## Clean DGP (Part 2)
`P = K-1 = 4`, `N = 500`, `L = 200`, `Phi ~ Dir(0.1)`. Score variances
`d = (1, .7, .5, .35)` (all gaps `> 0.1`), `sigma_eps^2 = .15`,
`B0 = R0 diag(sqrt(d - sigma_eps^2))`, so `eig(Cov(z)) = d`. The
**pseudo-true estimand** `B_star` is the standardized spectral fit on one
large pilot corpus (`M_pilot`). Block 1 reports three variance estimates —
analytic SE (`sgscatm_vcov`), across-replicate empirical SD (gold standard),
and delete-block jackknife — plus RMSE, coverage, and `ratio_to_null`.

## Lambda rule & diagnostic (Part 3)
Data-driven, truth-free: `lambda_A = ` the `(K-1)`-th eigenvalue of the word
Gram `W~ W~^T` (== `sigma_{K-1}`), so the covariate-augmentation block is
commensurate with the smallest retained word-signal eigenvalue. This exact
rule is used by the paper and the BES application. `lambda_diag.R` re-runs
the grid at `lambda = 1` and `lambda = lambda_A` and prints the LAMBDA
VERDICT (penalty asymptotically inert ⇒ Section 4's lambda-free variance is
safe).

## Outputs
- `data/` — `B0.rds`, `R0.rds`, `score_d.rds`, `B_star.rds`, `lambda_A.rds`,
  `block{1..5}.rds`, `lambda_diag.rds` (raw per-replicate, saved before aggregation).
- `tables/` — `block{1..5}.csv/.tex`, `lambda_diagnostic.csv`.
- `imgs/` — `block1_{rmse_vs_M,coverage_vs_M,se_vs_M,qqplot}.pdf`,
  `block2_linearisation.pdf`, `block3_{mse_boxplot,timing}.pdf`,
  `block5_{disagree_vs_snr,bestgain_vs_M,ndistinct}.pdf`.

## Verdict
`SE_SHRINKS` (slope of analytic SE ∈[−.65,−.35]), `SE_MATCHES_EMPIRICAL`
(∈[.85,1.20]), `SE_MATCHES_JACKKNIFE` (∈[.8,1.25]), `COVERAGE` (∈[.92,.97]);
OVERALL "SE FORMULA VALIDATED" iff those four pass. Plus `RMSE_DECREASES`,
`NOT_COLLAPSED`, `DEGENERACY_GONE`, the LAMBDA VERDICT, and a full
figure/table inventory with byte sizes.
