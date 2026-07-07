## Table A2 — published pipeline: reproduction and mechanism

| M | published RMSE (50 reps) | reproduced RMSE (20 reps) | zero-est. level | RMSE/zero | norm ratio ||Bz||/||Bz0|| | coverage (published SE) | coverage (SE/sqrt(M)) | med CI half-width | med |center-truth| |
|---|---|---|---|---|---|---|---|---|---|
| 500 | 0.2387 | 0.2387 | 0.2598 | 0.92 | 0.086 | 1.000 | 1.000 | 47.800 | 0.214 |
| 1000 | 0.2450 | 0.2449 | 0.2598 | 0.94 | 0.061 | 1.000 | 0.988 | 26.011 | 0.217 |
| 2000 | 0.2495 | 0.2496 | 0.2598 | 0.96 | 0.042 | 1.000 | 0.650 | 11.742 | 0.220 |


### Bias/variance decomposition (two-step, aligned scale)

| M | bias^2 | variance | bias^2/mse | mean(|bias|/SE) |
|---|---|---|---|---|
| 500 | 9.20e-04 | 3.71e-04 | 0.71 | 0.91 |
| 1000 | 7.41e-04 | 1.93e-04 | 0.79 | 1.16 |
| 2000 | 4.32e-04 | 1.20e-04 | 0.78 | 1.26 |

Log-log slope of the sqrt(variance) component on M: **-0.409** (3 points) — the variance obeys the M^{-1/2} law; the bias component is M-independent (finite-L floor at L = 200).

## Table A3 — corrected two-step estimator (k = 5, lambda = 0)

| M | RMSE (two-step) | RMSE (pilot, published conv.) | norm ratio | entrywise coverage | row-norm coverage | Armijo failures | monotone |
|---|---|---|---|---|---|---|---|
| 500 | 0.0358 | 0.2386 | 1.103 | 0.867 | 0.700 | 0 | 50/50 |
| 1000 | 0.0305 | 0.2450 | 1.093 | 0.765 | 0.427 | 0 | 50/50 |
| 2000 | 0.0234 | 0.2495 | 1.070 | 0.712 | 0.373 | 0 | 50/50 |

Log-log slope of per-replicate RMSE on M: **-0.295** (95% CI [-0.355, -0.235]); M^{-1/2} reference = -0.5.

## Table B3 — fair STM comparison (paired with basin_check E2 data)

| regime | method | mse_Bz paper (Procrustes) | norm ratio ||B||/||Bz0|| | GL/oracle mse | time (s) | EM its |
|---|---|---|---|---|---|---|
| weak | pilot (published pipeline) | 0.0050 (0.0015) | 0.193 | 0.0001 (0.0001) | 0.9 | — |
| weak | pilot + refined (k=5, l=0) | 0.0002 (0.0001) | 0.891 | 0.0003 (0.0001) | 1.4 | — |
| weak | STM (native gamma, ALR->ILR) | 0.0084 (0.0034) | 2.086 | — | 6.0 | 23 |
| weak | STM (theta -> ILR -> OLS, old worker) | 0.0085 (0.0035) | 2.090 | 0.0002 (0.0001) | 6.0 | 23 |
| strong | pilot (published pipeline) | 0.0768 (0.0177) | 0.055 | 0.0001 (0.0001) | 0.7 | — |
| strong | pilot + refined (k=5, l=0) | 0.0021 (0.0012) | 1.135 | 0.0024 (0.0012) | 1.1 | — |
| strong | STM (native gamma, ALR->ILR) | 0.0220 (0.0153) | 1.497 | — | 7.6 | 34 |
| strong | STM (theta -> ILR -> OLS, old worker) | 0.0220 (0.0153) | 1.498 | 0.0002 (0.0001) | 7.6 | 34 |
| strong M=5000 | pilot (published pipeline) | 0.0753 (0.0186) | 0.024 | 0.0002 (0.0001) | 0.7 | — |
| strong M=5000 | pilot + refined (k=5, l=0) | 0.0013 (0.0010) | 1.088 | 0.0014 (0.0011) | 2.8 | — |
| strong M=5000 | STM (native gamma, ALR->ILR) | 0.0136 (0.0042) | 1.411 | — | 33.7 | 28 |
| strong M=5000 | STM (theta -> ILR -> OLS, old worker) | 0.0136 (0.0042) | 1.412 | 0.0000 (0.0000) | 33.7 | 28 |

Published Table 3 references: weak 0.009, strong 0.021 (M = 1000).

## Table C — operating-regime map (M = 1000, 10 reps per cell)

| b_max | sat(theta_true) | sat(theta_end) | pilot GL mse (oracle) | pilot paper mse | refined paper mse | Armijo fails | max nu | monotone |
|---|---|---|---|---|---|---|---|---|
| 0.15 | 0.000 | 0.000 | 0.0001 (0.0000) | 0.0052 (0.0011) | 0.0003 (0.0002) | 0 | 1.00e-06 | 10/10 |
| 0.50 | 0.000 | 0.000 | 0.0001 (0.0001) | 0.0906 (0.0242) | 0.0036 (0.0015) | 0 | 1.00e-06 | 10/10 |
| 1.00 | 0.000 | 0.001 | 0.0002 (0.0001) | 0.2806 (0.0699) | 0.0218 (0.0093) | 0 | 1.00e-06 | 10/10 |
| 1.50 | 0.004 | 0.016 | 0.0024 (0.0019) | 0.7327 (0.1996) | 0.0523 (0.0259) | 0 | 1.00e-06 | 10/10 |

