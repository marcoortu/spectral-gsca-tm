# Basin-condition verification for the Newton refinement (Block 3 setting)

Numerical validation of the two-stage architecture (spectral pilot →
block Gauss–Newton refinement of the exact objective) planned for the
Newton–Kantorovich theorem. Design, code and results live in
`replication/basin_check/`; all numbers below are reproducible via

```
Rscript replication/basin_check/02_run_experiments.R          # ~25 min on 10 workers
Rscript replication/basin_check/03_report.R
```

Full configuration: `M = 1000`, `N = 500`, `K = 5`, `P = 3`,
`sigma_eps = 0.3`, `alpha_beta = 0.1`, `doc_length = 200`, signal levels
`weak: b_max = 0.15`, `strong: b_max = 0.50`; 20 replicates (E1/E2),
5 replicates (E3/E4), 10 replicates (E5).

---

## Executive summary

1. **The refinement works, dramatically.** From the pilot, the exact
   λ = 0 refinement reaches mse(B̂z) = **0.0003 (weak) / 0.0004
   (strong)** under the paper metric — equal to the refined-from-truth
   oracle floor to within 1e-4, and **20–50× better than the published
   STM reference** (0.009 / 0.021). Prediction P2 is exceeded: the
   pilot→STM gap in the strong regime is closed 137 % (pilot 0.0768 →
   refined 0.0004).
2. **The Table 3 pilot–STM gap is an alignment artifact, not an
   estimation gap.** `fit$Z` has unit-norm columns, so `fit$Bz` is
   ~5–18× too small in norm; orthogonal Procrustes cannot rescale, and
   the reported pilot MSE (0.077 strong) is ≈ E‖Bz0‖²/12 = b²/3
   (0.083), i.e. the MSE of the zero matrix. Under the general-linear
   (GL) alignment the pilot subspace is essentially perfect
   (mse ≈ 0.0001). P4 confirmed.
3. **The basin story must be reformulated before proving anything.**
   The exact λ = 0 criterion has (i) an *exactly* 20-dimensional gauge
   null space at the M-estimator η\*, as predicted (P3, confirmed by
   direct compression Q'HQ, |eigs| ≤ 8e-8 vs λmax ≈ 400), but also
   (ii) γ_perp that is *positive yet tiny* (5.8e-6 weak, ~1e-7 strong —
   at strong signal numerically indistinguishable from the gauge
   block), and (iii) H(η₀) at the finite-sample truth is slightly
   **indefinite** (γ_perp ≈ −1e-3). With ρ_pilot ≈ 9–21 (comparable to
   ‖Z_true‖_F!), the Kantorovich ratio r = 2·L_H·ρ/γ is 1e4–1e9 —
   the basin condition in the raw Euclidean (Z, Φ) metric certifies
   nothing, even though the refinement demonstrably converges from the
   pilot to the same model (endpoint fitted matrices agree to ~1e-3–1e-4
   relative, objective values to ~1e-7). The parameters (Z, Φ) are only
   weakly identified; the functional B(Z) = (C'C)⁻¹C'Z is what is
   sharply determined (start-independence of mse to ≤ 5e-4). **The
   theorem should be stated for the gauge quotient and, realistically,
   for the B functional (or the fitted matrix Θ(Z)Φ), not for raw
   (Z, Φ).**
4. **The λ = 1 exact objective is mis-scaled and harmful.** The penalty
   ‖Z − CB‖² at the truth (≈ 360) dominates the reconstruction term
   (≈ 5) by ~70×, so the exact λ = 1 minimiser collapses Z into span(C)
   at essentially zero reconstruction cost (multinomial noise floor)
   and lets B chase noise: refined-λ=1 mse is 0.0071 (weak) / 0.065
   (strong) — *worse than the pilot*. The spectral objective's λ = 1
   and the exact objective's λ = 1 live on different scales; the
   refinement theory should use λ = 0 or a noise-matched
   λ = O(1/(L·σ_ε²·(K−1))).

A failed prediction is reported as such below (P1 strict basin
criterion, P3 at the truth anchor); nothing was tuned to rescue them.

---

## 1. Alignment audit (mandatory deliverable)

Traced verbatim from `R/sgscatm_fit.R` (`sgscatm()`), commit state of
this working tree.

**How `fit$Bz` is produced:**

1. `scale_W = TRUE` (default, used by Block 3): `W ← W / rowSums(W)` —
   the DTM is converted to **relative frequencies**. All basin-check
   objectives use the same convention (`Wf`).
2. `C ← scale(C, center = TRUE, scale = FALSE)`.
3. `W̃ = W − 1 w̄'` (column-centred), `w̄ = colMeans(W)`.
4. Truncated SVD of W̃ with `r = min(M, N, 100)` then the safety cap
   `r ← min(r, M−1, N, K−1+P)`; for the Block 3 design **r = 7**.
5. `H = [U_r diag(σ_r), √λ Q_C]` with `Q_C = qr.Q(qr(C))`
   (M × (r+P) = M × 10).
6. Eigendecomposition of `H'H`; top K−1 = 4 eigenpairs `(E_top, s_top)`.
7. `Z* = H E_top diag(1/√s_top)` — therefore **Z*'Z* = I₄: the score
   columns are unit-norm eigenvectors of S_z = HH'**. They are *not*
   eigenvalue-scaled and *not* √M-scaled. This is the crux of finding 2.
8. `Ψ̂ = K · Z*' W̃`; **varimax** (`normalize = FALSE`) on `t(Ψ̂)` gives
   the orthogonal `R_star`; `Z* ← Z* R_star`, `Ψ̂ ← R_star' Ψ̂`.
   Varimax is applied **before** the regression on C (orthogonal, so it
   does not affect the identification issue).
9. `Φ̂ = V Ψ̂ + 1_K w̄'`, `Π̂ = softmax_rows(Z* V')`.
10. **`Bz = solve(C'C, C' Z*)`** — OLS of the unit-norm, varimax-rotated
    scores on C.

**What `procrustes_align()` does** (`replication/simulation/sim_utils.R`):
`R = U V'` from `svd(Bz_hat' Bz0)`; **orthogonal only — no scaling, no
translation** (reflections allowed; no det +1 constraint);
`mse = mean((Bz_hat R − Bz0)²)`.

**Consequences (measured, Table A):** `‖fit$Bz‖/‖Bz0‖` = 0.19 (weak) /
0.055 (strong). With `Bz_hat ≈ 0`, mse ≈ mean(Bz0²); with
Bz0 ~ U(−b, b), E mean(Bz0²) = b²/3 = 0.0075 (weak), 0.0833 (strong) —
the observed paper-metric pilot values (0.0050, 0.0768) sit essentially
at this "zero-estimator" level. The GL correction A_hat is strongly
**anisotropic** (mean singular values 8.7/6.6/3.7/1.4 weak,
26/17.8/8.9/3.6 strong), so no single rescaling — and no orthogonal map
— absorbs it: the general-linear part is real and eigenvalue-dependent.
The same caveat applies to the STM column of Table 3 (STM's prevalence
coefficients live on their own scale); a scale-invariant metric is
recommended for the revision.

**Table A — pilot accuracy under the three alignments (M = 1000, 20 reps, mean (sd))**

| regime | mse paper (Procrustes on fit$Bz) | mse, GL-aligned Z | mse, OP-aligned Z | ρ_pil GL (Z / Φ parts) | ρ_pil OP |
|---|---|---|---|---|---|
| weak   | 0.0050 (0.0015) | 0.0001 (0.0001) | 0.0051 (0.0014) | 17.5 (17.5 / 0.59) | 20.4 |
| strong | 0.0768 (0.0177) | 0.0001 (0.0001) | 0.0769 (0.0176) | 16.1 (16.1 / 0.16) | 35.9 |

Two readings matter for the theory: (i) GL ≤ OP everywhere (P4 ✓), and
(ii) even after oracle GL alignment **ρ_pil ≈ 16–17 with ‖Z_true‖_F ≈
21**: the pilot recovers the C-predictable component of Z but not the
per-document ILR residuals ε_i — the pilot is *not* a small-ball
estimate of (Z, Φ); it is only close in the B-relevant directions.
(The GL-aligned pilot mse of 0.0001 is *oracle-optimistic* — the
alignment borrows Z_true; the feasible number is the paper-metric one.)

---

## 2. Implementation and verification

- Exact objective, both variants, on `Wf` (frequencies):
  `F0 = ‖Wf − Θ(Z)Φ‖²_F`, and for λ > 0 the **profiled**
  `F1 = F0 + λ‖(I − P_C)Z‖²_F` (B has the closed form
  `B(Z) = (C'C)⁻¹C'Z`; envelope theorem gives the gradient).
- Analytic gradients verified against `numDeriv::grad` on the tiny
  instance (M = 30, N = 40, K = 3): rel. error **4.8e-9 (λ=0),
  1.8e-8 (λ=1)** — gate threshold 1e-6, passed.
- HVPs by central differences, `h = 1e-5·(1+‖η‖)/‖v‖`, verified against
  the full `numDeriv` Hessian on the tiny instance: rel. error ≤ 8e-10.
- Gauge tangency `g·t = 0` verified to 4e-18 (exact invariance of F0
  along the orbit, any point).
- Refinement: sweep = Z-step (2 damped GN iterations/doc, Armijo with
  ≤ 30 halvings, per-doc Levenberg ν from 1e-6, ×10 on failure) →
  Φ-step (ridge 1e-8 LS) → B-step. Convergence: relative F decrease
  < 1e-10 AND max_i ‖∇_{z_i}F‖_∞ < 1e-7·(1+F); cap 100 sweeps (λ = 0).
  **Monotone decrease of F held in all 190 refinement runs** (tolerance
  1e-12·(1+|F|)).
- λ = 1 uses a two-block variant with B profiled *inside* the Z-step
  (Woodbury rank-P(K−1) correction; single global Armijo). Measured
  reason: the brief's three-block descent zig-zags along the Z ≈ CB
  valley at rate 0.9886/sweep (>1100 sweeps to tolerance); on the same
  instance the profiled solver reaches F = 0.9988 in 100 sweeps vs
  1.0827 for the three-block after 400. Cap 600 sweeps (tail rate still
  ≈ 0.98 — the λ = 1 landscape is nearly flat, see §5).

## 3. E1 — operational basin check

**Table E1 (20 reps per cell; "p/t" = from pilot / from truth)**

| regime | λ | same basin (strict) | same basin (gauge-aware) | med dZ_rel | med dFit_rel | med dEtaPerp_rel | med dF_rel | med sweeps (p) | converged p/t |
|---|---|---|---|---|---|---|---|---|---|
| weak   | 0 | 0/20 | 4/20 | 5.9e-01 | 3.8e-04 | 2.5e-01 | 9.5e-08 | 37  | 17/2 |
| weak   | 1 | 0/20 | 0/20 | 5.4e-01 | 3.5e-06 | 5.6e-05 | 5.0e-06 | 600 | 0/0 |
| strong | 0 | 0/20 | 0/20 | 3.5e-01 | 1.1e-03 | 2.0e-01 | 6.5e-07 | 100 | 0/20 |
| strong | 1 | 0/20 | 0/20 | 4.9e-02 | 1.8e-04 | 1.3e-04 | 3.8e-05 | 600 | 0/1 |

- **The strict pre-registered criterion (dZ_rel < 1e-3 ∧ dF_rel < 1e-8)
  fails everywhere — P1 fails as stated** — but for an instructive
  reason, not because the starts land in different valleys: at λ = 0
  the minimiser is a ≥20-dimensional gauge orbit, so endpoint Z's
  cannot coincide; the gauge-invariant fitted matrices Θ(Z)Φ agree to
  **3.8e-4 (weak) / 1.1e-3 (strong) relative**, objective values to
  **1e-7**, and the deliverable B̂z agrees between the two starts to
  |Δmse| ≤ 4.8e-4 (median 1.6e-4 weak, 6.6e-5 strong) — i.e. the same
  model and the same estimator, different gauge/soft-mode
  representative. The formal gauge-aware criterion still mostly misses
  its dF < 1e-8 threshold because the 100-sweep cap leaves rel. F
  residuals of ~1e-7 (from-truth runs converged only 2/20 weak, 0/20
  strong within the cap).
- The remaining parameter discrepancy dEtaPerp ≈ 0.2–0.25 at λ = 0
  survives first-order gauge projection: consistent with E3's finding
  of *additional* near-flat non-gauge directions (γ_perp ≈ 1e-7–6e-6
  vs λmax ≈ 400). The (Z, Φ) parameters are weakly identified; B(Z)
  is not.

## 4. E2 — accuracy

**Table E2 — mse(B̂z), mean (sd), M = 1000, 20 reps**

| regime | estimator | paper metric (Procrustes) | GL / direct metric |
|---|---|---|---|
| weak   | pilot                          | 0.0050 (0.0015) | 0.0001 (0.0001)* |
| weak   | pilot + refined (λ = 0)        | **0.0003 (0.0002)** | 0.0003 (0.0002) |
| weak   | pilot + refined (λ = 1)        | 0.0071 (0.0018) | 0.0071 (0.0018) |
| weak   | refined from truth (λ = 0)     | 0.0002 (0.0001) | 0.0003 (0.0001) |
| weak   | refined from truth (λ = 1)     | 0.0071 (0.0018) | 0.0071 (0.0018) |
| weak   | STM (published Table 3)        | 0.0090 | — |
| strong | pilot                          | 0.0768 (0.0177) | 0.0001 (0.0001)* |
| strong | pilot + refined (λ = 0)        | **0.0004 (0.0002)** | 0.0006 (0.0003) |
| strong | pilot + refined (λ = 1)        | 0.0652 (0.0292) | 0.0669 (0.0275) |
| strong | refined from truth (λ = 0)     | 0.0004 (0.0002) | 0.0006 (0.0003) |
| strong | refined from truth (λ = 1)     | 0.0623 (0.0321) | 0.0641 (0.0307) |
| strong | STM (published Table 3)        | 0.0210 | — |

\* oracle-optimistic (GL alignment uses Z_true); the feasible pilot
number is the paper-metric column.

- **P2: pass, exceeded.** Strong regime: pilot 0.0768 → refined 0.0004
  = 137 % of the pilot→STM gap closed (target ≥ 70 %, "≈ 0.02 = full
  success"); the refined estimator equals the from-truth oracle floor.
  In sweep terms, mse(B̂z) stabilises within 0–3 sweeps (max 50) —
  far below the "≤ 15 median" bar — although *F*-convergence takes
  ~40–100 sweeps (the extra sweeps polish Φ/Z directions irrelevant
  to B).
- **P1: falsified in a favourable direction.** Weak regime was
  predicted "refined ≈ pilot ≈ STM (0.006–0.009)"; actually refined =
  0.0003 ≪ pilot(paper) = 0.0050 ≤ STM = 0.009. Even at weak signal the
  refinement buys ~16×.
- λ = 1 refinements *hurt* (see §5 and Executive summary point 4);
  identical numbers from both starts confirm it is the objective, not
  the optimisation.
- STM reference: published Table 3 values; the dual-R subprocess
  machinery was not run (allowed by the brief). Caveat: those STM
  numbers inherit the same orthogonal-only alignment convention.

## 5. E3 — Hessian spectrum at η₀ and at η\*

**Table E3 (5 reps per regime, means; γ = γ_perp at λ=0, raw min at λ=1)**

| regime | λ | anchor | λmax | γ | max\|Q'HQ\|/λmax | max‖Hq‖/λmax | sep. ratio γ/max\|gauge\| |
|---|---|---|---|---|---|---|---|
| weak   | 0 | truth η₀ | 400.0 | **−1.09e-03** | 4.6e-07 | 7.7e-05 | −6.3 |
| weak   | 0 | M-est. η\* | 400.1 | **+5.81e-06** | 7.8e-11 | 3.8e-09 | **330** |
| weak   | 1 | truth η₀ | 400.0 | +7.9e-04 | 5.0e-03 | 5.0e-03 | ≈ 0.16 |
| weak   | 1 | M-est. η\* | 1790 | −1.0e-05 | 2.9e-04 | 2.8e-04 | ≈ 0 |
| strong | 0 | truth η₀ | 400.9 | **−6.06e-04** | 3.0e-07 | 7.7e-05 | −5.2 |
| strong | 0 | M-est. η\* | 401.4 | **+1.20e-07** | 1.5e-10 | 5.0e-09 | **3.0** |
| strong | 1 | truth η₀ | 400.9 | +2.5e-04 | 5.0e-03 | 4.8e-03 | ≈ 0.05 |
| strong | 1 | M-est. η\* | 1120 | −1.2e-05 | 2.5e-03 | 2.3e-03 | ≈ 0 |

Figure: `results/spectrum_tail.png` (raw tail, deflated tail, gauge
block, and λ = 1 tail at η\*, replicate 1 per regime).

- **Gauge structure (P3), at η\*:** exactly as predicted — the 20
  columns of Q_T satisfy ‖H q‖/λmax ≤ 5e-9 and the compressed gauge
  block Q'HQ has |eigenvalues| ≤ 8e-8 (absolute, λmax ≈ 400): a
  numerically exact 20-dimensional null space aligned with the gauge.
  The pre-registered principal-cosine test is *not* usable as stated:
  ARPACK cannot resolve a 20-fold degenerate cluster from one Krylov
  sequence (it returns 3–6 copies; cosines for those are 1.000, the
  rest are missed). The direct Q'HQ compression is the correct
  instrument and confirms P3's geometric claim.
- **But γ_perp is tiny and, at the truth, negative.** At η₀ the λ = 0
  Hessian is indefinite (γ_perp ≈ −1e-3): the truth is not a stationary
  point at finite sample, and the gauge-orbit curvature term
  (H t = −(∂t)g) bleeds the gradient into the "null" directions
  (‖Hq‖/λmax ≈ 8e-5 at η₀ vs 5e-9 at η\*). At η\*, γ_perp > 0 but of
  order 1e-6 (weak) down to 1e-7–1e-9 (strong; one replicate −1.5e-9,
  i.e. zero to numerical precision): beyond the exact gauge there are
  *additional* near-flat non-gauge directions, more so at strong
  signal. P3's "γ_perp clearly positive and separated" holds cleanly in
  the weak regime at η\* (separation ×330), only marginally at strong
  (×3), and fails at η₀. 
- **λ = 1 lifts the gauge but does not fix identification:** at η₀ the
  gauge block rises to 5e-3·λmax (lifted, as predicted), yet at the
  λ = 1 "minimiser" λmax inflates 3–4× (Φ drifts along soft modes) and
  the smallest eigenvalues remain ≈ 0 (−1e-5): the penalty constrains
  12 = P(K−1) directions of Z but leaves the (Z, Φ) degeneracy intact.

## 6. E4 — Kantorovich diagnostic

**Table E4 (5 reps, means; ρ_perp = gauge-projected pilot displacement)**

| regime | λ | L_H | ρ_perp(η₀) | ρ_perp(η\*) | γ(η₀) | γ(η\*) | r(η₀) | r(η\*) |
|---|---|---|---|---|---|---|---|---|
| weak   | 0 | 0.731 (0.027) | 9.0 | 15.9 | −1.1e-03 | 5.8e-06 | −1.3e+04 | 4.5e+06 |
| weak   | 1 | 0.731 (0.027) | 9.0 | 11.9 | 7.9e-04 | −1.0e-05 | 1.7e+04 | −1.1e+07 |
| strong | 0 | 0.660 (0.110) | 12.6 | 20.7 | −6.1e-04 | 1.2e-07 | −3.3e+04 | −3.1e+09 |
| strong | 1 | 0.660 (0.110) | 12.6 | 20.8 | 2.5e-04 | −1.2e-05 | 6.9e+04 | −8.5e+07 |

The basin ratio r = 2·L_H·ρ_perp/γ is 4–9 orders of magnitude above 1
(or negative where γ < 0). Per the brief this outcome is admissible
(Kantorovich is sufficient, not necessary), but the magnitude carries a
structural message beyond conservatism: **in the raw Euclidean (Z, Φ)
metric the theorem's premises cannot hold in this model class.**
ρ_perp ≈ ‖Z_true‖_F because the pilot carries no information about the
per-document ILR residuals ε_i; γ is pinned near zero by weakly
identified directions (rare-word Φ columns, near-degenerate topics);
and η\* itself sits 20–35 Frobenius units from η₀ (the exact M-estimator
overfits per-document scores at L = 200 — consistent with the E5 1/L
floor). None of this prevents the *B functional* from being determined
to 4 decimals (E1/E2). Recommended reformulations, in decreasing order
of ambition: (i) NK on the gauge quotient with a Hessian-weighted
(Fisher-type) metric; (ii) NK for the profiled/reduced criterion in
B alone, treating (Z, Φ) as an infinite-dimensional-style nuisance with
its own rate; (iii) a direct two-step argument: pilot consistency of
span/projection + one-step (or k-step) Newton estimator theory for B̂,
which matches what E1/E2 actually show and does not need a basin in
(Z, Φ) at all.

## 7. E5 — finite-L bias floor (weak regime, λ = 0, refined from truth)

| doc_length L | mse paper | mse direct | mse ratio vs 1/L ratio (ref L = 200) |
|---|---|---|---|
| 50   | 0.0010 (0.0006) | 0.0021 (0.0012) | 6.6 vs 4.0 |
| 200  | 0.0002 (0.0001) | 0.0004 (0.0003) | 1.0 vs 1.0 |
| 1000 | 0.0001 (0.0000) | 0.0001 (0.0000) | 0.42 vs 0.20 |

Roughly 1/L scaling from 50→200 (slightly faster), flattening toward an
L-independent floor (~6e-5, M-driven) by L = 1000 — consistent with an
`a + b/L` error decomposition. Figure: `results/e5_biasfloor.png`.

## 8. Verdicts on the pre-registered predictions

- **P1 — fail (favourably).** Strict same-basin fraction is 0 %, not
  ~100 % — the criterion is unsatisfiable at λ = 0 (gauge orbit) and
  the gauge-aware version still trips over the 1e-8 dF threshold given
  the 100-sweep cap; the substantive content (start-independence of the
  estimator) holds: |Δmse| between starts ≤ 4.8e-4. "Refined ≈ pilot ≈
  STM" is false: refined (0.0003) beats both by an order of magnitude.
- **P2 — pass, exceeded.** 137 % of the strong-regime gap closed
  (0.0768 → 0.0004 ≪ 0.035; ≈ oracle floor); mse stabilises in ≤ 3
  sweeps (median), though full F-convergence needs 40–100.
- **P3 — pass at η\* for the gauge geometry, fail at η₀ and for
  separation at strong signal.** Null space is exactly 20-dim and
  gauge-aligned at η\* (via Q'HQ; the cosine test is an ARPACK
  artifact); γ_perp > 0 but ~1e-6 (weak) / ~1e-7 (strong, marginal);
  H(η₀) indefinite; λ = 1 lifts the gauge block (×1e4) as predicted but
  inflates λmax and keeps near-zero directions.
- **P4 — pass.** GL pilot mse (0.0001) ≤ Procrustes pilot mse
  (0.0050/0.0768); GL–OP gap enormous; scaling chain documented in §1
  (unit-norm eigenvector scores → varimax → OLS → orthogonal-only
  alignment).

## 9. Deviations from the brief (all logged, none post-hoc tuned)

1. **Platform:** `parallel::makePSOCKcluster` + `parLapplyLB` instead of
   `mclapply` (Windows: mclapply is serial). 10 workers.
2. **Bz0 draw:** drawn *inside* `sim_dgp()` under the replicate seed
   (`Bz0 = NULL`, `b_max` per regime), mirroring Block 3's
   per-replicate random Bz0 but fully reproducible from the seed alone
   (Block 3 draws it from the ambient RNG stream *before* seeding
   `sim_dgp`). Discrepancy flag from the brief stands: the paper text
   describes a fixed Bz0 (`run_simulation.R` Block 1 uses one), Block 3
   code randomises per replicate — we follow the Block 3 code.
3. **Monotonicity violation policy:** flag-and-record (+ warning)
   instead of hard `stop()`, so parallel replicates keep diagnostics.
   No violation occurred in any of the 190 runs.
4. **λ = 1 solver:** profiled two-block Gauss–Newton (Woodbury) instead
   of the brief's three-block; the three-block variant is retained
   (`profile_B = FALSE`) and its measured valley rate (0.9886/sweep) is
   the documented reason. Sweep cap 600 (vs 100) at λ = 1; λ = 0 keeps
   the brief's cap (100; 200 for the η\* anchor runs in E3).
5. **E1 gauge-aware endpoint criteria added** (dFit on Θ(Z)Φ,
   gauge-projected dEta): the pre-registered dZ criterion is
   gauge-blind at λ = 0. Strict results still reported first.
6. **E3 anchored at η\* in addition to η₀**: H(η₀) is indefinite at
   finite sample, so the NK γ is measured where the theorem needs it.
   Added the Q'HQ gauge-block compression after finding that ARPACK
   under-resolves the 20-fold degenerate cluster (affects the
   pre-registered principal-cosine diagnostic).
7. **Lanczos settings:** `RSpectra::eigs_sym` shifted operator,
   tol 1e-7, maxitr 1000, ncv = max(4k, 80); k = 30.
8. **STM:** published Table 3 values (no dual-R subprocess), as the
   brief allows.
9. Gradient tolerance choice: `tol_g = 1e-7·(1+|F|)` ("scaled
   tolerance" was left open by the brief).
10. E5 run with 10 reps (optional budget available).

## 10. Runtime and reproducibility

- Machine: Windows 11, 22 logical cores, R 4.5.1; 10 PSOCK workers.
- Wall-clock: derivative gate 0.5 s; E0 smoke 133 s; E1/E2 824 s;
  E3/E4 453 s; E5 20 s; report ≈ 5 s. **Total ≈ 24 min** (budget 3 h).
- Seed scheme: replicate seed = `BASE + regime_index·1000 + rep` with
  weak = 1, strong = 2; BASE = 90000 (E1–E4; E3/E4 reps 1–5 share data
  with E1/E2 reps 1–5), 80000 (E0 smoke), and E5 uses
  `70000 + dl_index·1000 + rep`. Derivative gate seed 4242. Power
  iterations seed at the replicate seed (+j for E4 segment pairs, +100
  for η\* analyses). Bz0, C, ε, Beta, W are all functions of the
  replicate seed via `sim_dgp()`.
- Files: `01_functions.R` (objective/gradients/GN/HVP/gauge/Lanczos),
  `02_run_experiments.R` (E0–E5, `--quick`, `--no-e5`, `--ncores N`),
  `03_report.R` (tables.md, summary.csv, figures),
  `results/*.rds` (raw per-replicate output incl. full per-sweep
  traces), `results/summary.csv`, `results/tables.md`,
  `results/spectrum_tail.{png,pdf}`, `results/e5_biasfloor.png`,
  `results/full_run.log`.
- No package function or existing script was modified.
