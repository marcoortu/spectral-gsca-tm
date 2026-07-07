## Table F3 — MSE(B_hat) by sweep count k (oracle start, unconstrained)

| regime | k=3 | k=5 | k=10 | k=20 | k=50 | k=100 | rule sweeps (med) | rule MSE |
|---|---|---|---|---|---|---|---|---|
| weak | 0.00017 | 0.00021 | 0.00024 | 0.00026 | 0.00026 | 0.00026 | 15 | 0.00026 |
| strong | 0.00216 | 0.00211 | 0.00174 | 0.00116 | 0.00056 | 0.00038 | 50 | 0.00056 |

### F3 diagnosis (strong): sweeps to reach 1.1x final MSE vs replicate traits

- slow tercile mean(relgap, rho_GL, nBz0): 0.022, 15.849, 1.022 vs fast tercile: 0.027, 16.372, 0.984
- correlations of sweeps-to-1.1x with (relgap, rho_GL, sat, nBz0): 0.03, -0.57,  NA, 0.26

### Criterion pathology — simplex-constrained descent FROM THE TRUTH

| regime | F(truth) | F sweep 1/10/50/100 | mse sweep 1/10/50/100 | monotone |
|---|---|---|---|---|
| weak | 4.9710 | 4.7571 / 4.7294 / 4.7268 / 4.7266 | 0.00013 / 0.00235 / 0.01193 / 0.01472 | 5/5 |
| strong | 4.9643 | 4.7765 / 4.7634 / 4.7622 / 4.7622 | 0.00075 / 0.01800 / 0.06428 / 0.07763 | 5/5 |

F decreases monotonically while mse_Bz rises by two orders: the exact frequency-LS criterion (constrained or not) does not identify B at finite L — the estimator must be k-step, not argmin.

## Table F1c — feasible variants vs oracle reference (M = 1000, 20 reps)

| regime | variant | mse paper | mse permutation | norm ratio | sweeps (med) | time (s) |
|---|---|---|---|---|---|---|
| weak | V1 | 0.04710 (0.01845) | 0.05333 (0.02148) | 3.54 | — | 0.8 |
| weak | V2 | 0.02869 (0.01139) | 0.03033 (0.01191) | 2.98 | 50 | 4.3 |
| weak | V2_jk | 0.01791 (0.00922) | 0.01959 (0.00984) | 2.53 | — | — |
| weak | V3 | 0.07381 (0.03740) | 0.07703 (0.03813) | 4.16 | 50 | 4.8 |
| weak | V3_jk | 0.05647 (0.03398) | 0.05988 (0.03451) | 3.71 | — | — |
| weak | V4 | 0.02008 (0.01064) | 0.02352 (0.01230) | 2.62 | 50 | 4.4 |
| weak | V4_jk | 0.01396 (0.00954) | 0.01715 (0.01117) | 2.30 | — | — |
| weak | oracle_k5 | 0.00021 (0.00009) | — | 0.89 | — | 0.4 |
| weak | oracle_rule | 0.00026 (0.00015) | — | 0.88 | 15 | 1.2 |
| weak | STM (native gamma; audit B3) | 0.00844 (0.00343) | — | 2.09 | — | 6.0 |
| weak | anchored Phi TV (context) | 0.342 (0.038) | | | | |
| strong | V1 | 0.07488 (0.07865) | 0.09446 (0.09352) | 1.78 | — | 0.7 |
| strong | V2 | 0.09039 (0.03642) | 0.09474 (0.04406) | 2.02 | 50 | 4.1 |
| strong | V2_jk | 0.05333 (0.03105) | 0.05819 (0.04005) | 1.77 | — | — |
| strong | V3 | 0.43181 (0.21465) | 0.44605 (0.22915) | 3.23 | 50 | 4.7 |
| strong | V3_jk | 0.35242 (0.19688) | 0.36692 (0.21163) | 2.99 | — | — |
| strong | V4 | 0.05119 (0.04259) | 0.06259 (0.05111) | 1.70 | 50 | 4.3 |
| strong | V4_jk | 0.03397 (0.03148) | 0.04497 (0.03991) | 1.53 | — | — |
| strong | oracle_k5 | 0.00211 (0.00119) | — | 1.13 | — | 0.4 |
| strong | oracle_rule | 0.00056 (0.00032) | — | 1.05 | 50 | 3.1 |
| strong | STM (native gamma; audit B3) | 0.02197 (0.01530) | — | 1.50 | — | 7.6 |
| strong | anchored Phi TV (context) | 0.210 (0.075) | | | | |

## Table F1d — identification boundary (strong regime, M = 1000)

| alpha_beta | true exclusivity | anchor TV | V2 mse paper | V2 mse perm |
|---|---|---|---|---|
| 0.05 | 1.000 | 0.232 (0.078) | 0.0954 (0.0510) | 0.1008 (0.0562) |
| 0.10 | 1.000 | 0.228 (0.060) | 0.0933 (0.0480) | 0.0995 (0.0601) |
| 0.30 | 0.964 | 0.252 (0.045) | 0.0816 (0.0244) | 0.0849 (0.0252) |
| 1.00 | 0.744 | 0.386 (0.040) | 0.0447 (0.0142) | 0.0616 (0.0186) |

## Table F2 — split-document jackknife (oracle start, Block 1 grid)

| M | RMSE uncorrected | RMSE jackknifed | cov entry (unc) | cov entry (jk) | cov entry (jk, inflated) | cov rownorm (unc) | cov rownorm (jk) | cov rownorm (jk, infl) |
|---|---|---|---|---|---|---|---|---|
| 500 | 0.0272 | 0.0242 | 0.942 | 0.940 | 0.970 | 0.913 | 0.853 | 0.933 |
| 1000 | 0.0205 | 0.0237 | 0.935 | 0.878 | 0.938 | 0.900 | 0.687 | 0.780 |
| 2000 | 0.0139 | 0.0255 | 0.942 | 0.707 | 0.823 | 0.860 | 0.193 | 0.307 |
| 4000 | 0.0135 | 0.0242 | 0.788 | 0.507 | 0.595 | 0.640 | 0.080 | 0.160 |

Log-log RMSE slopes: uncorrected **-0.369** [-0.422, -0.317]; jackknifed **0.021** [-0.030, 0.072] (target -0.5).

### F2 bias/variance decomposition (aligned scale)

| M | bias^2 share (uncorrected) | bias^2 share (jackknifed) |
|---|---|---|
| 500 | 0.56 | 0.44 |
| 1000 | 0.55 | 0.69 |
| 2000 | 0.45 | 0.82 |
| 4000 | 0.59 | 0.86 |

## Table F2-L — jackknife vs L (weak regime, M = 1000, oracle start)

| L | mse uncorrected | mse jackknifed | reduction |
|---|---|---|---|
| 50 | 0.00202 | 0.00286 | 0.71x |
| 200 | 0.00026 | 0.00037 | 0.71x |
| 1000 | 0.00008 | 0.00008 | 1.00x |

## Table F4(ii) — feasible chain (init V4) on the Block 1 grid

| M | RMSE (jk) | RMSE perm (jk) | cov entry | cov entry (infl) | cov rownorm | cov rownorm (infl) |
|---|---|---|---|---|---|---|
| 500 | 0.3291 | 0.3554 | 0.380 | 0.420 | 0.007 | 0.013 |
| 1000 | 0.2175 | 0.2411 | 0.388 | 0.425 | 0.040 | 0.053 |
| 2000 | 0.1245 | 0.1429 | 0.453 | 0.497 | 0.160 | 0.193 |

Feasible-chain RMSE slope: **-0.828** [-0.950, -0.706].

