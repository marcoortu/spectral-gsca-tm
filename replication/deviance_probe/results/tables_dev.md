## Table P1 — truth-start pathology probe: deviance vs LS on identical data

| regime | criterion | mse sweep 1/10/50/100 | norm ratio @100 | gauge @100 | tilt (perp) @100 | monotone |
|---|---|---|---|---|---|---|
| weak | deviance | 0.0001 / 0.0001 / 0.0007 / 0.0037 | 1.67 | 16.7 | 49.1 | 5/5 |
| weak | LS (constrained) | 0.0001 / 0.0040 / 0.0128 / 0.0150 | 2.41 | 32.8 | 75.8 | 5/5 |
| weak | (mse at sweep 0 = truth) | 0.00006 | | | | |
| strong | deviance | 0.0002 / 0.0052 / 0.0391 / 0.1169 | 2.22 | 60.7 | 103.0 | 5/5 |
| strong | LS (constrained) | 0.0008 / 0.0216 / 0.0657 / 0.0782 | 2.01 | 42.8 | 78.8 | 5/5 |
| strong | (mse at sweep 0 = truth) | 0.00005 | | | | |

Gate P1 (dev mse@100 <= 2x mse@0): evaluated per replicate in the run log; ratios are reported there and in REPORT_DEV.md.

## Table P2a — oracle-start k curves: deviance vs LS (identical data)

| regime | criterion | k=1 | k=3 | k=5 | k=10 | k=20 | k=50 | k=100 | rule stop (med) | rule MSE |
|---|---|---|---|---|---|---|---|---|---|---|
| weak | deviance | 0.00406 | 0.00709 | 0.00929 | 0.01167 | 0.00962 | 0.00896 | 0.01779 | 100 | 0.01779 |
| weak | LS (F3)  | 0.00006 | 0.00017 | 0.00021 | 0.00024 | 0.00026 | 0.00026 | 0.00026 | 15 | 0.00026 |
| strong | deviance | 0.00856 | 0.02697 | 0.05166 | 0.10060 | 0.12367 | 0.17692 | 0.27911 | 100 | 0.27911 |
| strong | LS (F3)  | 0.00140 | 0.00216 | 0.00211 | 0.00174 | 0.00116 | 0.00056 | 0.00038 | 50 | 0.00056 |

## Table P3 — anchor polish + feasible deviance chain (M = 1000)

| regime | anchor TV | polished TV | chain mse paper | chain mse perm | norm ratio | coverage | time (s) |
|---|---|---|---|---|---|---|---|
| weak | 0.342 (0.038) | 0.247 (0.020) | 0.0582 (0.0276) | 0.0617 (0.0287) | 3.80 | 0.417 | 52.8 |
| weak | (refs) STM 0.0084 | LS-V4+jk 0.0140 | oracle-LS rule 0.0003 | | | | |
| strong | 0.210 (0.075) | 0.143 (0.052) | 0.1260 (0.0755) | 0.1334 (0.0792) | 2.16 | 0.300 | 53.4 |
| strong | (refs) STM 0.0220 | LS-V4+jk 0.0340 | oracle-LS rule 0.0006 | | | | |
| strong M=5000 | 0.079 (0.021) | 0.062 (0.012) | 0.0015 (0.0019) | 0.0020 (0.0020) | 1.07 | 0.608 | 132.4 |
| strong M=5000 | (ref) STM 0.0136 | | | | | | |

