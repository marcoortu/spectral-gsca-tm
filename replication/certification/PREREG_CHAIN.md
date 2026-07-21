# Pre-registration — Realign package to the Theorem-16 estimator (chain + Lemma 17 sandwich)

**Written:** 2026-07-08, before executing any rerun. Gates and predictions below are **frozen**.
Failed gates are reported as findings with their mechanism; nothing is retro-fitted; the chain is
**not** forced to reproduce the current draft's raw-pilot numbers.

**Package:** `sgscatm`, loaded via `devtools::load_all(".")`.

---

## Objective

Make the package's **primary** estimator and SE coincide with the manuscript's theory
(Sections 3–4): a three-stage estimator
**spectral pilot → anchored general-linear orientation → k-step least-squares Gauss–Newton
refinement (frozen Φ, λ=0, B-stationarity stop)**, with inference from the
**heteroskedasticity-robust sandwich of Lemma 17**. Then re-run and certify every manuscript
number against the chain. The current package ships only the raw spectral pilot
(`sgscatm()`, step 10 `Bz = solve(C'C, C'Z*)` = the degenerate raw-pilot regression of
Corollary 11) plus `sgscatm_vcov()` (an influence SE on the *standardized* score scale): a
different estimator and a different SE object than the theory.

## Non-negotiables (any violation is a bug)

1. z-step refinement criterion is **least squares**, never multinomial deviance. Deviance appears
   in exactly one place: the pooled EM polish of the anchored Φ̂ (a Φ-only step).
2. Orientation alignment is **general-linear** (√M scale + anisotropic (K−1)×(K−1) map), never
   orthogonal Procrustes.
3. Refinement uses **λ = 0** in the exact objective (λ belongs to the pilot only).
4. **Do NOT force the chain to reproduce the current draft's application numbers.** The chain
   produces different numbers; the manuscript is updated to the chain's certified values.
5. **Reuse, do not reimplement.** Promote the validated code in `replication/feasibility/`
   (anchors, GL alignment V4), `replication/basin_check/` and `replication/deviance_probe/`
   (LS k-step refinement, deviance Φ-polish, trust-region guard, B-stationarity rule) into `R/`.

## Estimator (`sgscatm_chain()`, new `R/chain.R`)

- **Stage 1 — pilot** (reuse `sgscatm(..., rotate=FALSE)` internals): W̃ = column-centered
  row-normalized W; truncated randomized SVD; H = [UΣ | √λ Q_C]; top K−1 eigenvectors → Z\*
  (unit-norm); expose O(1)-scale word-Gram eigenpairs ρ_k = eig(M⁻¹W̃′W̃).
- **Stage 2 — anchored orientation** (promote `feasibility/01_anchors.R`): factorial
  co-occurrence Qc; SPA anchors; constrained LS → Φ̂₀; pooled multinomial-EM Φ-polish → Φ̂;
  per-document read-out z⁰ = Vᵀ ln Π_{δ0}((Φ̂Φ̂ᵀ)⁻¹Φ̂ w) then damped per-document GN on
  ‖w − Φ̂ᵀf(z)‖² with trust-region cap ‖δz‖ ≤ 1 → ẑᴬ; general-linear alignment
  A_c = argmin_A ‖√M Z\*A − H_M Ẑ_A‖²_F; oriented pilot Ẑ₀ = √M Z\*A_c + 1_M z̄_Aᵀ.
- **Stage 3 — k-step LS refinement** (promote `basin_check/` + `deviance_probe/`): damped
  per-document GN on the **exact** ‖w − Φ̂ᵀf(z)‖² at frozen Φ̂, λ=0, start Ẑ₀; B-stationarity
  stop (max|ΔB̂z|/rms(B̂z) < η=1e-3 on two consecutive sweeps, cap k̄=100).
- **Output** (`class "sgscatm_chain"`): oriented generative-scale `Bz`, `Theta=F(Ẑ)`, anchored
  `Phi`, refined `Z`, stored `W_tilde`, `C_centred`, residuals; raw-pilot `Bz` and Stage-1
  subspace kept as fields (Table 1 collapse demo).

## Inference (`R/sandwich.R`, `vcov.sgscatm_chain`)

Refined residuals r̂_i = ẑ_i − B̂zᵀc_i, then
Σ̂_B = (M/(M−P))·(I⊗(CᵀC)⁻¹)·[Σ_i (r̂_ir̂_iᵀ)⊗(c_ic_iᵀ)]·(I⊗(CᵀC)⁻¹);
SE = sqrt(diag Σ̂_B); Wald from Σ̂_B. Estimates the variance of B̂z around its own mean
(includes √M/L and √M·δ_Φ drifts), so coverage of Bz,0 is nominal only when both drifts vanish
(Corollary 18). Keep `sgscatm_vcov()` (raw-pilot influence) and `ilr_se()` (bootstrap) as
reference SEs for G4.

## Identification / alignment conventions

- In the estimator: orientation/scale from anchors (general-linear); B̂z is generative-scale and
  oriented up to **topic permutation + sign** (Theorem 8), not continuous rotation.
- Simulation metrics: align B̂z to Bz,0 by **permutation + sign only** — never continuous
  Procrustes, never a general-linear map to truth. (The old `sim_utils::procrustes_align` must
  not be used for the chain's reported RMSE/coverage.)
- Real-data bootstrap (G4): refit chain on resampled units, permutation+sign align each replicate
  to the point estimate, empirical SD.
- Oracle = chain given true orientation (skip stage 2); feasible = full stage-2 anchors.

## Frozen gates (on the CHAIN)

- **G0** subspace recovery + raw collapse (sim). PASS = principal angle < 1e-3 at all M; raw ratio
  ‖B̂z^pilot‖/‖Bz,0‖ ∝ M^{−1/2}; chain RMSE within ~2× oracle at largest M.
- **G1** delocalization (real CRC). r = M·max‖ẑ_i‖²/(K−1) small & stable across prevalence at
  genus, materially larger at species; report both ranks.
- **G2a** sandwich SE / empirical SD ∈ [0.8,1.25] across cells including high-√M/L (variance
  calibrated even where coverage degrades).
- **G2b** coverage of Bz,0 ∈ [0.90,0.97] in small-drift regime and **degrades** as √M/L grows
  (confirmed Corollary-18 prediction = PASS). ≈0.95 at √M/L≲0.22 → ≈0.6 at 0.35. Coverage pinned
  at 1.000 = FAIL (scale collapse).
- **G2c (linchpin)** in-regime chain calibration (sim, LARGE L ∈ {1e4,1e5}, alpha_beta≈0.05,
  M ∈ {500,2000,5000}): chain + sandwich achieve nominal coverage of Bz,0. PASS = coverage
  ∈ [0.90,0.97], SE/SD ∈ [0.9,1.15]. The experiment the repo currently lacks.
- **G3** crossover vs STM (sim): feasible trails STM at small M, beats at large M; oracle
  dominates. Match Table 3 direction.
- **G4 (CRITICAL)** sandwich SE vs document bootstrap (≥200, perm+sign aligned) on real CRC:
  median ratio ∈ [0.8,1.25] AND ≥60% entries within ±25%. Report sandwich, sgscatm_vcov, ilr_se.
- **G5** concordance (real CRC): (a) chain joint Wald agrees with PERMANOVA on study_condition;
  (b) top disease-loading genera overlap known CRC set above chance (hypergeometric p<0.05).
- **G6** speed (real CRC): chain fit+sandwich ≤ 1/10 of ALDEx2, STM, PERMANOVA.
- **G7** B-functional start-independence (sim): chain B̂z from pilot-start vs truth-start agree
  (mse diff ≤ 5e-4).
- **G8** bias field b(z) + fixed-L gradient G0 (sim): MC + Richardson matches closed-form b(z);
  CRN finite-diff matches G0; ‖G0‖_F ≫ MC SE; corr(log‖G0,·j‖, log ω_j) ≈ 0.98.

## Experiments

A. Simulation: Table 1 (recovery+collapse, M∈{1000,2000,5000}); Table 2 (coverage vs √M/L,
L=200, M∈{1000,2000,4000,5000}, ≥50 reps); **Table 2-largeL** (G2c); Table 3 (crossover vs STM);
Fig 1 (SE/SD vs b_max); b(z)/G0 validation.
B. Microbiome CRC (swap raw pilot → chain in 06/07/08): Table 4 (delocalization genus+species),
Table 5 (Wald vs PERMANOVA), Table 6 (speed), Fig 2 (forest, sandwich CI), Fig 3 (SE vs bootstrap),
Fig 4 (known genera).
C. BES (swap → chain): Table 7 / Fig 5 (recovery strong vs weak, chain vs STM, runtime);
out of regime (√M/L≈86) → recovery+speed only, confirm coverage gate fails as predicted.

## Protocol

Minimal edits to existing `R/` (leave `sgscatm()` pilot intact, log changes); chain in new files;
roxygen → NAMESPACE/man; tests in tests/testthat/. Fixed seeds (base + regime·1000 + rep), PSOCK
workers, reproduce-from-scripts. Log every deviation with before/after. Monotone-descent + stability
checks on every refinement run. Data: Zeller 2014 via SIAMCAT (actual source; retry
curatedMetagenomicData optionally). Deliverables per brief §10; bottom line: does chain + sandwich
achieve G2c and G3, and do G1/G4/G5/G6 pass on the chain?
