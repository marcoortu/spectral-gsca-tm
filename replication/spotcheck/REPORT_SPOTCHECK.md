# Section-4 algebra spot-check — b(z) and G0

Pre-integration verification of the draft's two hand-derived objects.
Reproduce via

```
Rscript replication/spotcheck/02_run.R            # ~20 min on 10 workers
Rscript replication/spotcheck/02_run.R --only-a   # SC-A rerun path
Rscript replication/spotcheck/03_report.R
```

## Gate verdicts (one line each)

- **SC0 (unit gates): PASS** — H2theta vs numDeriv 1.4e-10 (< 1e-8);
  H2f 2.5e-10 (< 1e-6); optimized-vs-naive objects 5.4e-15; certified
  GN from truth passes (see deviation 1 on the certificate level).
- **Gate A as registered (1/sqrt(L)-fit intercepts + bhat_800,
  max(5%, 3 se)): FAIL** — 10/36 intercepts and 10/36 bhat_800
  components miss, all explainable (below). **Under the corrected
  remainder rate the algebra passes 36/36**: the Richardson
  extrapolation on the clean cells (2·bhat_800 − bhat_400) matches
  closed-form b(z) at every component of every test point
  (e.g. A_u6 comp 2: 23.03 vs 23.036). **The Lemma's b(z) algebra is
  CONFIRMED; the draft's O(L^{-1/2})-remainder claim is WRONG — the
  remainder is O(1/L)** (multinomial third cumulants are O(L^{-2});
  independently re-derived and empirically confirmed by the L-grid).
- **Gate B1 (||E[G0]||_F > 10 MC se): PASS** — ratios 748 (weak) and
  198 (strong). **But the draft's localization claim is inverted**:
  ||G0[, j]|| scales WITH word frequency (corr of logs = 0.98); for
  unit directions confined to the 50 lowest-frequency columns,
  <G0, U> is ~300× smaller than for generic directions. The
  step-criterion lemma's "action on low-frequency columns" must be
  corrected.
- **Gate B2 as registered (every direction within max(10%, 3 se),
  both h, both L): FAIL at 1 of 24 cells** (direction 4, L = 200:
  Dhat = −2.1e-6 ± 3.2e-6 vs formula +1.21e-5; wrong sign, 4.4 se).
  All other cells pass, and after the control-variate sharpening the
  informative cells agree at ratios 0.87–1.37. The one failure is
  quantitatively consistent with an O(L^{-2}) remainder with a large
  direction-4 coefficient (c_U ≈ −0.56 fitted at L = 200 predicts the
  observed L = 400 value within 1 se, where the cell PASSES at ratio
  0.71). **G0's formula and the envelope derivation are validated
  asymptotically; the draft should state the finite-L remainder as
  O(L^{-2}) with direction-dependent constants.**
- **P-C (optional cosine > 0): FAIL/inconclusive** — cosines
  (−0.013, 0.023, −0.005, 0.034, 0.008), mean 0.009, at the random
  noise floor 1/sqrt(K·N) ≈ 0.02 for these 2500-dim objects. The
  predicted alignment of the sweep-10 LS Phi-drift with −G0 is not
  detectable.

**Bottom line for the manuscript:** both hand-derived objects are
algebraically right; two *claims around them* are wrong and now get
fixed pre-referee — the b(z) remainder order (O(1/L), not
O(L^{-1/2}), with the expansion outside its validity range at
L ≲ 100 for ||z|| ≥ 1), and the G0 low-frequency-localization
statement (the action sits on high-frequency columns; corr 0.98).

## SC-A — bias field b(z)

Full componentwise table in `results/tables_spotcheck.md`; figure
`results/a_bias_extrapolation.png`. Summary of the four estimates per
component (36 components = 9 test points × 4):

| estimate | gate passes | comment |
|---|---|---|
| 1/sqrt(L)-fit intercept (registered) | 25/36 | misspecified: remainder is O(1/L) |
| 1/L-fit intercept (full grid)        | 33/36 | residual misses from L ≤ 200 contamination |
| Richardson 2·b800 − b400 (clean)     | **36/36** | the correct estimate under O(1/L) |
| bhat_800 raw                         | 26/36 | fails only where c/800 > 5% of b (||z|| ≥ 1) |

- **Localization** (means of |bhat_800 − x|/se): full formula b —
  0.5–7 se; b1 alone — 18–140 se; −0.5·H⁻¹a2 alone — 16–71 se. Both
  halves of the algebra are needed and both are right; nothing to fix
  in the formula.
- **Where the registered gate broke, mechanically:** (i) the
  L·bias expansion's remainder is c(z)/L with c growing steeply in
  ||z|| (c/b ≈ 3–40 at ||z|| = 1–1.5 in units of 1/L — at L = 800
  that is a genuine 5–18% offset in bhat_800); (ii) at L = 100 the
  certificate excludes up to 12% of draws at ||z|| = 1.5 (extreme
  tails), conditioning the sample and bending the low-L end of any
  full-grid fit; (iii) the 1/sqrt(L) regressor then extrapolates
  through the wrong curvature.
- **Brute-force cell** (K = 3, N = 30, L = 50, ||z|| = 1, R = 300k,
  no CV): raw (15.75, −13.86) and CV (18.97, −15.44) estimates agree
  in order but differ by ~15 se — 6.3% of draws are excluded even at
  the 1e-8 certificate, and E[m | certified] ≠ 0 shifts the two
  estimators differently. Both sit at ~2.2× the closed form: at
  L = 50 the second-order term is as large as b itself. The
  registered brute agreement gate fails **for the same reason the
  corrected-rate analysis succeeds** — L = 50 at ||z|| = 1 is outside
  the first-order expansion's radius. The draft should carry an
  explicit validity condition (empirically: c(z)/L ≪ b(z), violated
  below L ≈ 100–200 at ||z|| ≥ 1).

## SC-B — gradient direction G0

Route 1 (`results/tables_spotcheck.md`): ||E[G0]||_F = 0.0694 (weak) /
0.0757 (strong), MC se ratios 748 / 198 — decisively nonzero (P-B1
PASS). Top columns by norm have word marginals ≈ 0.02 vs corpus
median 0.00088; corr(log ||G0[, j]||, log p_j) = 0.98 — the
localization claim inverts (see verdicts).

Route 2 (figure `results/b_fd_vs_formula.png`): central FD of the
profiled criterion with common random numbers, sharpened post hoc by
the analytic control variate c_r = −2⟨e_r, θ(z_r)'U⟩ (deterministically
recomputable from the seeded base draws; reduces se 3–5×, from
~1e-5 — at which level the raw test was vacuous, se ≈ |formula| — to
1.4e-6–3.2e-6). Results: h-stability exact to 3 digits between
h = 1e-3 and 1e-4; informative Gaussian directions agree at ratios
0.87–1.37; low-frequency directions have formula values (≈4e-8) below
the sharpened se (≈1e-7) — consistent, low power, and itself the
demonstration that G0 has no low-frequency action. The 1/L-scaling
ratios are ≈2 where the formula dominates the remainder (dirs 1, 3)
and noisy elsewhere. The single hard failure (dir 4, L = 200) is
attributed to the O(L^{-2}) term as quantified in the verdicts.

## SC-C — optional cosine check

The feasibility RDS does not store Phi paths, so the sweep-10 LS
pathology displacement was recomputed on the same seeds (5 strong
reps, 10 constrained-LS sweeps from the truth). Cosines with −E[G0]
are at the noise floor (mean 0.009). Not necessarily in tension with
SC-B: G0 is the gradient of the population per-document-refit
profiled criterion at Phi0, while the probe's displacement reflects
10 sweeps of empirical two-GN-step block descent with simplex
projection at M = 1000 — path curvature and the projection evidently
decorrelate the two long before sweep 10. The draft should not lean
on this alignment for intuition.

## Deviations (all logged, none post-hoc tuned)

1. **GN certificate 1e-12 → 1e-8** (registered 1e-12). Two measured
   reasons: (i) BFGS confirms rows stalling at gmax ≈ 7e-11 are at the
   optimum to 7 digits — 1e-12 is below float granularity at loss
   ≈ 1e-2, so it excludes converged rows; (ii) exclusions bias the
   estimators through E[m | certified] ≠ 0 (at 1e-10 the brute cell
   excluded 9.8% and raw-vs-CV disagreed by 50 se; at 1e-8 exclusions
   drop to ≤ 0.7% at L ≥ 400). z-accuracy at gmax = 1e-8 is ~1e-5,
   three orders below the smallest measured bias increment.
2. **Additional SC-A analyses beyond the registered fit**: the 1/L
   fit and the clean-cell Richardson extrapolation, motivated by the
   independent re-derivation of the remainder order before any data
   were seen (noted in `01_formulas.R`'s header). The registered
   1/sqrt(L) gate is still reported and evaluated literally.
3. **SC-B2 control variate added at the analysis stage** (no rerun:
   the base draws are seeded and deterministically reproducible).
   Raw CRN columns are reported alongside.
4. **SC-C recomputed** rather than loaded (Phi paths not stored in
   feasibility RDS); same seeds, 10 sweeps.
5. `--only-a` flag added to rerun SC-A/brute after the certificate
   fix without repeating SC-B/SC-C.
6. Windows: PSOCK workers (10), as in all previous rounds.

## Runtime and seeds

- Wall-clock (10 workers): SC0 0.2 s; SC-A 361 s (+ 346 s certificate
  rerun); brute 152 s (+ 168 s rerun); SC-B1 1.4 s; SC-B2 681 s;
  SC-C 2.6 s; report ≈ 30 s. **Total ≈ 28 min including the rerun**
  (budget 1 h).
- Certificates: 720,000 SC-A solves with 7,168 exclusions (1.0%,
  concentrated at L = 100, ||z|| ≥ 1 — reported per cell in the
  results); 300,000 brute solves, 6.3% excluded; SC-B2 exclusions
  handled pairwise under CRN.
- Seeds: Phi0-A 77001, Phi0-B 77002, test directions 77010, SC-A
  draws 77000 + point_idx·100 + L_idx, brute 77777 / 77700 + chunk,
  SC-B route 1 78000 + regime, SC-B2 base draws 78200 + L (z) and
  78300 + L (counts), directions 78100, SC0 77000, SC-C
  90000 + 2000 + rep (feasibility E2, strong).
- Files: `01_formulas.R` (objects, naive reference, certified batch
  solver, SC0), `02_run.R` (`--quick`, `--only-a`, `--ncores`),
  `03_report.R`; `results/*.rds`, `results/summary_spotcheck.csv`,
  `results/tables_spotcheck.md`, figures `a_bias_extrapolation.png`,
  `b_fd_vs_formula.png`, logs `full_run.log`, `rerun_a.log`.
- No package function or pre-existing script modified.
