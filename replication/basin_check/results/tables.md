## Table A — Pilot accuracy under the three alignments

| regime | reps | mse_Bz paper (Procrustes on fit$Bz) | mse_Bz GL-aligned Z | mse_Bz OP-aligned Z | rho_pil GL (Z / Phi) | rho_pil OP |
|---|---|---|---|---|---|---|
| weak | 20 | 0.0050 (0.0015) | 0.0001 (0.0001) | 0.0051 (0.0014) | 17.52 (17.51 / 0.59) | 20.35 |
| strong | 20 | 0.0768 (0.0177) | 0.0001 (0.0001) | 0.0769 (0.0176) | 16.06 (16.06 / 0.16) | 35.85 |

rho values are means over replicates of
sqrt(||Z_pil - Z_true||_F^2 + ||Phi_pil - Phi_true||_F^2).

## Table E1 — Endpoint coincidence (refine from pilot vs from truth)

| regime | lambda | same basin (strict) | same basin (gauge-aware) | med dZ_rel | med dFit_rel | med dEtaPerp_rel | med dF_rel | med sweeps (pilot) | converged p/t |
|---|---|---|---|---|---|---|---|---|---|
| weak | 0 | 0/20 | 4/20 | 5.86e-01 | 3.80e-04 | 2.52e-01 | 9.50e-08 | 37 | 17/2 |
| weak | 1 | 0/20 | 0/20 | 5.44e-01 | 3.50e-06 | 5.56e-05 | 5.02e-06 | 600 | 0/0 |
| strong | 0 | 0/20 | 0/20 | 3.48e-01 | 1.08e-03 | 1.99e-01 | 6.47e-07 | 100 | 0/0 |
| strong | 1 | 0/20 | 0/20 | 4.93e-02 | 1.81e-04 | 1.31e-04 | 3.82e-05 | 600 | 0/1 |

Strict criterion: dZ_rel < 1e-3 AND dF_rel < 1e-8 (pre-registered). Gauge-aware: dFit_rel < 1e-3 AND dF_rel < 1e-8, where dFit compares the gauge-invariant fitted matrices Theta(Z)Phi.

Monotone F decrease held in **all** refinement runs: yes.

## Table E2 — mse_Bz, mean (sd) over replicates

| regime | estimator | paper metric (Procrustes) | GL / direct metric |
|---|---|---|---|
| weak | pilot (paper alignment) | 0.0050 (0.0015) | 0.0001 (0.0001) |
| weak | pilot + refined (lambda = 0) | 0.0003 (0.0002) | 0.0003 (0.0002) |
| weak | pilot + refined (lambda = 1) | 0.0071 (0.0018) | 0.0071 (0.0018) |
| weak | refined from truth (lambda = 0, oracle floor) | 0.0002 (0.0001) | 0.0003 (0.0001) |
| weak | refined from truth (lambda = 1) | 0.0071 (0.0018) | 0.0071 (0.0018) |
| weak | STM (published Table 3, M = 1000) | 0.0090 | — |
| strong | pilot (paper alignment) | 0.0768 (0.0177) | 0.0001 (0.0001) |
| strong | pilot + refined (lambda = 0) | 0.0004 (0.0002) | 0.0006 (0.0003) |
| strong | pilot + refined (lambda = 1) | 0.0652 (0.0292) | 0.0669 (0.0275) |
| strong | refined from truth (lambda = 0, oracle floor) | 0.0004 (0.0002) | 0.0006 (0.0003) |
| strong | refined from truth (lambda = 1) | 0.0623 (0.0321) | 0.0641 (0.0307) |
| strong | STM (published Table 3, M = 1000) | 0.0210 | — |

Pilot GL column: OLS of the GL-aligned scores on C, direct entry-wise MSE (the GL alignment already absorbs rotation and scale). Refined estimates live in model coordinates, so their two columns differ only by the final Procrustes rotation.

## Table E3 — Hessian spectrum diagnostics (5 reps per regime)

| regime | lambda | anchor | lambda_max | gamma (perp at l=0 / raw min at l=1) | max |QtHQ| rel | max ||H q||/lmax | cosines > 0.99 (of 20) | sep. ratio |
|---|---|---|---|---|---|---|---|
| weak | 0 | truth eta0 | 400.0 | -1.09e-03 | 4.64e-07 | 7.66e-05 | 0.0 | -6.3 |
| weak | 0 | M-est. eta* | 400.1 | 5.81e-06 | 7.79e-11 | 3.77e-09 | 5.6 | 330.3 |
| weak | 1 | truth eta0 | 400.0 | 7.86e-04 | 5.00e-03 | 4.97e-03 | 0.0 | 0.0 |
| weak | 1 | M-est. eta* | 1790.5 | -1.02e-05 | 2.90e-04 | 2.79e-04 | 6.2 | -0.0 |
| strong | 0 | truth eta0 | 400.9 | -6.06e-04 | 2.96e-07 | 7.68e-05 | 0.0 | -5.2 |
| strong | 0 | M-est. eta* | 401.4 | 1.20e-07 | 1.51e-10 | 4.99e-09 | 1.4 | 3.0 |
| strong | 1 | truth eta0 | 400.9 | 2.46e-04 | 4.99e-03 | 4.78e-03 | 0.0 | 0.0 |
| strong | 1 | M-est. eta* | 1120.5 | -1.18e-05 | 2.52e-03 | 2.27e-03 | 2.0 | -0.0 |

QtHQ = 20 x 20 compression of H onto the gauge basis; its eigenvalues measure curvature along the gauge orbit directly (ARPACK cannot resolve the 20-fold degenerate cluster; principal cosines of the raw Lanczos vectors undercount the null space for that reason). sep. ratio = gamma / max|gauge eigenvalue|: how cleanly the non-gauge curvature separates from the gauge block.

## Table E4 — Kantorovich diagnostic (same 5 reps)

| regime | lambda | L_H | rho_perp (from eta0) | rho_perp (from eta*) | gamma(eta0) | gamma(eta*) | r = 2 L_H rho_perp/gamma (eta0) | r (eta*) |
|---|---|---|---|---|---|---|---|---|
| weak | 0 | 0.731 (0.027) | 9.04 | 15.93 | -1.09e-03 | 5.81e-06 | -1.30e+04 | 4.53e+06 |
| weak | 1 | 0.731 (0.027) | 9.04 | 11.92 | 7.86e-04 | -1.02e-05 | 1.70e+04 | -1.13e+07 |
| strong | 0 | 0.660 (0.110) | 12.59 | 20.65 | -6.06e-04 | 1.20e-07 | -3.29e+04 | -3.09e+09 |
| strong | 1 | 0.660 (0.110) | 12.59 | 20.79 | 2.46e-04 | -1.18e-05 | 6.91e+04 | -8.47e+07 |


## Table E5 — refined-from-truth (lambda = 0) vs document length
(weak regime)

| doc_length | mse_Bz paper | mse_Bz direct | mse ratio vs 1/L ratio |
|---|---|---|---|
| 50 | 0.0010 (0.0006) | 0.0021 (0.0012) | 6.56 vs 4.00 |
| 200 | 0.0002 (0.0001) | 0.0004 (0.0003) | 1.00 vs 1.00 |
| 1000 | 0.0001 (0.0000) | 0.0001 (0.0000) | 0.42 vs 0.20 |

