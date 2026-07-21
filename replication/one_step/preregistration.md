# Preregistration — One-Step-Efficient Estimation of the Structural ILR Path Matrix `B_z0`

**Status:** written BEFORE any estimator run. Gates below are fixed. Failures are
findings, not bugs to be patched away. Deviations from this document are logged in
`report.md` with a timestamped rationale.

**Estimand.** The *generative* ILR path-coefficient matrix `B_z0` (P × (K−1)) that
`replication/simulation/sim_dgp.R` uses to generate topic proportions from
covariates — **not** the across-replicate pseudo-mean, and **not** the standardized
eigen-score coefficient. Consistency for `B_z0` is the point of this run.

No package `R/` code is modified. All new estimator code lives in
`replication/one_step/estimators.R`, marked experimental, and reuses `sim_dgp()`,
`sim_utils.R::procrustes_align()`, `sgscatm()`, `sgscatm_vcov()`, `ilr_se()`,
`ilr_se_analytical()`.

---

## 1. Conventions fixed from the package source (read first)

All symbols below are taken verbatim from `sim_dgp.R`, `sgscatm_fit.R`,
`ilr_contrast.R`, `vcov.R`. Everything downstream must match these.

### 1.1 ILR contrast matrix `V`
`ilr_contrast(K)` returns a **K × (K−1)** Helmert-style matrix with

* `V' V = I_{K−1}` (orthonormal columns),
* `V' 1_K = 0` (columns span `1_K^⊥`, the tangent space of the simplex at the centroid).

### 1.2 Generative model (`sim_dgp`)
For document `i` (`c_i` = P-vector of **column-centred** covariates, `E[C]=0` enforced
by `scale(..., center=TRUE)`):

```
z_i   = B_z0' c_i + eps_i,        eps_i ~ N(0, sigma_eps^2 I_{K-1})     # ILR scores  (M×(K-1) = Z_true)
theta_i = softmax(V z_i)          # topic proportions  (closure of exp(V z_i))   (M×K = Theta_true)
pi_i   = Beta' theta_i            # term probabilities  (Beta = K×N, rows sum to 1)
w_i   ~ Multinomial(L_i, pi_i)    # word counts  (W = M×N)
```

so `Z_true = C B_z0 + Eps` and `theta_i = softmax(V z_i)`. `eta_i := B_z0' c_i` is the
**structural (noise-free) ILR mean**; `p_i := theta_i` is the topic proportion vector.
`Sigma_C := E[c c']` (identity unless a wrapper overrides it). Document length `L_i`
(here written `L`) is the multinomial size; `M` documents, `N` vocabulary, `K` topics.

### 1.3 Estimator internals (`sgscatm`)
* `W` optionally row-normalised to term frequencies (`scale_W=TRUE`); `w_bar = colMeans(W)`;
  `W_tilde = W − 1 w_bar'` (M×N, centred).
* Truncated SVD → augmented factor `H = [U_r diag(sig_r), sqrt(lambda) Q_C]`;
  eigendecompose `H'H`.
* **Scores** `Z* = H E_top diag(s_top^{-1/2})`, so **`Z*' Z* = I_{K-1}`** — the scores are
  *unit-norm* eigenvectors (columns have sum-of-squares 1, i.e. per-document scale
  `O(M^{-1/2})`). This is the object stored as `fit$Z`.
* `Psi = K · Z*' W_tilde` ((K−1)×N); `Phi = V Psi + 1_K w_bar'` (K×N, **rows sum to 1**;
  `Psi = V' Phi`).
* `B_z = (C'C)^{-1} C' Z*` (P×(K−1)) — regression of *unit-norm* scores on covariates.
* Proportions `Pi = softmax(Z* V')`.
* Optional varimax rotates `Z*` and `Psi` by an orthogonal `R*` (identified only up to rotation, Thm 3).

`sgscatm_vcov()` reports SE on the **standardized** scale `Z~ = sqrt(M) Z*`
(`M^{-1} Z~'Z~ = I`), with `B` built from the O(1) word-Gram eigenpairs
`rho = eig(M^{-1} W_tilde' W_tilde)`, `ztil = (v_l' w~_i)/sqrt(rho)`.

### 1.4 The scale gap (why the baseline fails on `B_z0`)
`Z*` has total (not per-document) variance normalised to 1: `Cov` of a `Z*` column is
`~ 1/M`, whereas the true scores `z_i` have `Cov(z) = B_z0' Sigma_C B_z0 + sigma_eps^2 I`
of `O(1)`. Hence `Z* ≈ Z_true R S^{-1}` with `S = Cov(z)^{1/2}` (times an `M^{-1/2}`
factor absorbed by the normalisation), and

```
B_hat_std = (C'C)^{-1} C' Z*  ≈  B_z0 R S^{-1} · (scale ~ M^{-1/2}),
```

which collapses toward 0 on the `B_z0` scale and worsens as covariate strength (`b_max`)
grows `S`. This is the object of gate **G1**.

---

## 2. The linearised generative model (basis for `proj` and the one-step)

**Closure Jacobian at the centroid.** `d softmax(V eta)/d eta` at `eta=0` equals
`(1/K)(I_K − 1 1'/K) V = (1/K) V` (using `V'1=0`). So the loading from `eta` to the
**centred** proportion vector is `L = (1/K) V` (K×(K−1)), and near the centroid

```
p_i − 1/K·1_K  ≈  (1/K) V z_i.                                   (LIN)
```

**Proportion recovery is exact-linear (scale-correct), independent of L.** Term
probabilities are *exactly* linear in the proportions: `pi_i = Beta' theta_i`. With the
fitted topic-word matrix `Phi` (K×N, rows sum to 1) as the loading, the ordinary
right-inverse projection of the observed term frequencies `f_i = w_i/L_i` (row `i` of the
row-normalised `W`) is

```
theta_hat_i = (Phi Phi')^{-1} Phi f_i            (P_hat = F Phi' (Phi Phi')^{-1}, M×K)
```

which returns `theta_i` **exactly** whenever `f_i = Phi' theta_i` and `Phi` has rank K —
no linearisation, natural scale carried by `Phi`. Multinomial sampling injects
mean-zero noise of order `L^{-1/2}`, so `theta_hat_i` is unbiased for `theta_i` at any `L`.

**Loading-projection scores (`proj`).** Centre `P_tilde = P_hat − colMeans`, then invert
(LIN) by its centroid linearisation:

```
H_proj = P_tilde · L (L'L)^{-1} = K · P_tilde · V           (M×(K-1)),
B_proj = (C'C)^{-1} C' H_proj.
```

In the linearised, noiseless limit `H_proj` recovers `z_i` and `B_proj` recovers `B_z0`
(verified numerically before use — §5). The **only** approximation left is (LIN)'s
centroid linearisation, whose error is `O(||z_i||^2)` and grows with `b_max`. This is the
residual obstruction gate **G4** isolates.

**One-step Gauss-Newton correction (`onestep_*`).** Take exactly ONE Newton step from the
`sqrt(M)`-consistent start `eta_hat_i = H_proj_i`, re-linearising the closure map about
`eta_hat_i` (not the centroid). With `s_i = softmax(V eta_hat_i)`, closure Jacobian
`J_i = (diag(s_i) − s_i s_i') V` (K×(K−1)) and closure residual
`r_i = theta_hat_i − s_i` (K-vector),

```
onestep_uw:  eta_i^+ = eta_hat_i + (J_i' J_i)^{-1} J_i' r_i                 (Omega = I)
onestep_mw:  eta_i^+ = eta_hat_i + (J_i' Omega_i^- J_i + ridge)^{-1}
                                     J_i' Omega_i^- r_i
             Omega_i = (1/L_i)(diag(theta_hat_i) − theta_hat_i theta_hat_i')  # multinomial working cov
```

`Omega_i^-` = ridge-stabilised generalized inverse (drop the null direction `1_K`).
Then `B_onestep = (C'C)^{-1} C' H_new`, `H_new = [eta_i^+]`. This is a *k-step /
one-step-efficient* estimator: exactly ONE step, no argmin, barrier-robust (a document
whose `J_i` is near-singular contributes a bounded, ridge-damped update).

**Standard errors (proj, onestep_uw, onestep_mw).** White/HC0 sandwich on the FINAL ILR
regression `H = C B + resid`. Per-document residual `e_i = H_i − B' c_i` ((K−1)-vector):

```
vcov(vec B_hat) = (C'C)^{-1} ⊗ I  ·  [ sum_i (c_i c_i') ⊗ (e_i e_i') ]  ·  (C'C)^{-1} ⊗ I
SE = sqrt(diag(...)) reshaped to P×(K-1).
```

All comparisons to `B_z0` and all coverage statements are made **only after**
`procrustes_align()` (hard requirement); the sandwich covariance is rotated to the aligned
basis by `(R' ⊗ I_P)` exactly as `sgscatm_vcov()` does. `sgscatm_vcov()` (standardized
scale) and `ilr_se()` (bootstrap) are reported alongside where meaningful.

---

## 3. Estimand reconciliation (JASA referee point)

Reviewer's proposed population form: `Sigma_C^{-1} E(C' U_0)`. With `U_0 = Z_true` (the
true ILR scores) and `E[C]=0`:

```
E[c_i z_i'] = E[c_i (B_z0' c_i + eps_i)'] = E[c_i c_i'] B_z0 = Sigma_C B_z0   (eps ⟂ C)
```

so `Sigma_C^{-1} · M^{-1} E(C' U_0) = B_z0` **exactly** — the two coincide when `U_0` is the
generative ILR score. Note `E[C]=0` does **not** imply `E[C'U_0]=0`: `U_0` depends on `C`
through its mean `C B_z0`, so `Cov(C, U_0) = Sigma_C B_z0 ≠ 0`. If instead `U_0` is the
*standardized eigen-score* `Z~`, the population form returns `B_z0 R S^{-1}` (the G1
object), differing from `B_z0` by the scale `S` — precisely the gap this run closes.
**Verified numerically on the DGP** (`00_estimand_check.R`): `Sigma_C^{-1} M^{-1} C'Z_true`
vs `B_z0`, and the standardized version vs `B_z0`.

---

## 4. Gates (pass/fail fixed now)

| Gate | Claim | PASS criterion |
|------|-------|----------------|
| **G1** | scale problem reproduced | `baseline_std`: `RMSE(B_z0)/‖B_z0‖ ≈ 1` and **increasing in `b_max`** |
| **G2** | scale fix | `proj`: `RMSE(B_z0) ≪ ‖B_z0‖` at small `b_max` **and** direct RMSE ≈ rotation-aligned RMSE (scale recovered, not only rotation) |
| **G3** | SE calibration, no collapse | HC sandwich for `proj`: median `SE/SD ∈ [0.8,1.25]`, coverage of the across-replicate mean ≈ nominal, **never pinned at 1.000** |
| **G4** | linearisation bias isolated | coverage of `B_z0` declines with `b_max` while `SE/SD ≈ 1` and the decline tracks `‖mean_est − B_z0‖` growing `~ O(‖eta‖^2)` |
| **G5(a)** | one-step reduces RMSE | `onestep_mw` (and `_uw`) reduce `RMSE(B_z0)` vs `proj` across `b_max` (target **≥40%** in the moderate regime), and weighted ≥ unweighted |
| **G5(b)** | **CONSISTENCY (decisive)** | at fixed `b_max`, fixed `L`, `RMSE(B_z0)` of one-step **→ 0 as M grows** (`M∈{500,1000,2000,4000,8000}`); `proj`-only plateaus at the linearisation-bias floor |
| **G5(c)** | L-robustness (barrier check) | at fixed `M`, one-step residual bias does not grow as `L` shrinks (`L∈{50,100,200,400,1000}`) |
| **G6** | Gate-2 (Δ_Phi) clean probe | scale-**preserving** corruption of `Phi` (vary topic-recovery error without changing loading scale): report `RMSE(B_z0)` and coverage vs `Δ_Phi`, and the boundary |

**G5 PASS = (a) + (b).** (c) is reported as the barrier check. If **G5(b) FAILS** (one-step
does not drive `RMSE(B_z0)→0` as `M` grows) we say so plainly: a single step does not remove
the linearisation bias, consistency for `B_z0` remains open, and the method estimates a
standardized/linearised effect rather than the structural coefficient. That is decisive
information, reported, not hidden.

---

## 5. Pre-run sanity checks (must pass before any gate is trusted)

* **Noiseless recovery** (`sigma_eps→0`, `b_max` small, `L` huge e.g. 1e5, `M` large):
  `H_proj` recovers `Z_true` and `B_proj` recovers `B_z0` to `<1%` relative error after
  Procrustes. If this fails the `proj` derivation is wrong and NO gate is run.
* **Self-consistency:** `validate_dgp()` passes; `Phi` rows sum to 1; `Psi = V'Phi`.

## 5b. Deviation addendum (logged 2026-07-18, before gate runs)

The noiseless check (`00_noiseless_check.R`) revealed that projecting on the fit's OWN
`Phi` does **not** recover natural scale: `Phi` has the correct row-space (principal angles
to `V'Beta` ≈ 0.4°, 0.5°, 1.7° on the three well-recovered directions) but is
scale-inflated 7.3× in Frobenius norm because the unit-norm score normalisation
(`Z*'Z* = I`) folds the scale `S = Cov(z)^{1/2}` into `Phi`. Projecting on it reproduces the
collapsed baseline. This is the fundamental second-moment scale/basis indeterminacy: `S` is
identified only through the simplex / topic-word structure, which the spectral fit discards.

**Confirmed on the DGP:** projecting the counts on the TRUE topic-word matrix `Beta`
recovers `B_z0` to `rmse/‖B_z0‖ = 0.006` with scale ratio `‖B_al‖/‖B_z0‖ = 0.994`. The
loading-projection is therefore an **anchor-based** estimator: it consumes a natural-scale
topic-word anchor `Phi_anchor` (the estimated topics), and the score-scale fix operates
*given* that anchor.

**Consequence for the design (no gate weakened):** this cleanly decomposes the estimand
into (i) topic recovery (a consistent natural-scale `Beta` estimate) and (ii) the
score-scale + linearisation fix. Gates **G1–G5** supply a good anchor (the generative
`Beta`, standing in for a consistent topic estimator) so they isolate (ii) — the object of
this run. Gate **G6** was *designed* to degrade the anchor (`Δ_Phi`); it now directly
measures how topic-recovery error in (i) propagates to `B_z0`. This matches fact 2 ("derive
`L` from the package conventions": `L = (1/K)V` is the closure Jacobian; the topic-word
matrix enters only through the proportion read-out) and G6's stated purpose. The honest
caveat — proj/onestep require a consistent natural-scale topic estimate — is carried in
`report.md`.

## 6. Design grid

* Base: `K=5, N=500, P=3, alpha_beta=0.1, lambda=1, sigma_eps=0.3` unless a gate varies it.
* `b_max` sweep: `{0.1, 0.25, 0.5, 0.75, 1.0, 1.5}` (near-centroid → delocalised transition).
* `M` scaling (G5b): `{500,1000,2000,4000,8000}` at `b_max=0.5, L=200`.
* `L` robustness (G5c): `{50,100,200,400,1000}` at `M=2000, b_max=0.5`.
* `Δ_Phi` (G6): scale-preserving rotation/mixing of `Phi` rows, grid of corruption strength.
* Replicates: `n_rep=100` for RMSE/coverage gates (G1–G4, G6), `n_rep=50` for the heavier
  M- and L-sweeps (G5b, G5c). Seeds `10000+r` (matching `run_simulation`).
* Fixed decisions: HC0 (not HC3) for the primary SE; ridge `1e-6·tr(J'J)/(K-1)` for
  `onestep_mw`; generalized inverse via eigen-thresholding at `1e-8·max`. Deviations logged.
