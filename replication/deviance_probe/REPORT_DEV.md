# Deviance probe — does multinomial weighting cure D1?

Architecture-deciding round, building on
`replication/feasibility/REPORT_FEAS.md` (D1–D3). Reproduce via

```
Rscript replication/deviance_probe/02_run.R      # ~11 min on 10 workers
Rscript replication/deviance_probe/03_report.R
```

## Gate verdicts

- **Gate P1 (deviance descent from truth does not degrade B): FAIL.**
  Per-replicate mse@100 / mse@0 ratios: 42–121 (weak), 1384–7735
  (strong); required ≤ 2. The drift is slower than LS early (weak:
  4× lower at sweep 100; strong: crosses *above* LS after sweep ~80)
  but it does not stop, and the near-converged deviance solution sits
  12–200× above the oracle reference (0.0037/0.117 vs 0.0003/0.0006).
- **P-P1 (gauge absorbs the displacement; tilt shrinks ≥ 10× vs LS):
  FAIL.** Under both criteria the tilt (gauge-orthogonal) component
  dominates the gauge component at sweep 100 (deviance: 49 vs 17 weak,
  103 vs 61 strong; LS: 76 vs 33, 79 vs 43). The tilt shrinks only
  ~1.5× at weak and *grows* at strong.
- **Gate P2a (strong-regime transient disappears; rule stops < 50):
  FAIL — inverted.** Under deviance the oracle-start k-curve is
  *monotonically increasing* in both regimes (weak 0.0041 → 0.0178,
  strong 0.0086 → 0.279 over k = 1…100); the rule never triggers
  (median stop = cap 100). The LS k-step is 20–500× better at every k
  on identical data. P2b skipped per protocol (gate P1 failed).
- **Gate P3 soft (feasible-deviance ≤ 1.5× STM at M = 1000; ≤ STM at
  M = 5000): FAIL at M = 1000** (0.058/0.126 vs STM 0.0084/0.0220),
  **PASS at M = 5000** (0.0015 vs STM 0.0136 — 9× better, norm ratio
  1.07, coverage still low at 0.61).
- **P-P3 (polish TV < 0.15; norm-ratio < 1.3): partial.** The
  likelihood polish genuinely improves Φ̂ everywhere — TV 0.342→0.247
  (weak), 0.210→0.143 (strong, < 0.15 ✓), 0.079→0.062 (M = 5000) — but
  misses 0.15 at the weak working point, and the norm-ratio inflation
  at M = 1000 remains (3.8/2.2); it collapses to 1.07 only at M = 5000.

## P1 — the decision experiment

Figure `results/p1_pathology_overlay.png`; gauge-drift table in
`results/tables_dev.md`. On identical data (E2 seeds, reps 1–5):

| regime | criterion | mse at sweep 1 / 10 / 50 / 100 | tilt@100 / gauge@100 |
|---|---|---|---|
| weak   | deviance | 0.0001 / 0.0001 / 0.0007 / 0.0037 | 49 / 17 |
| weak   | LS       | 0.0001 / 0.0040 / 0.0128 / 0.0150 | 76 / 33 |
| strong | deviance | 0.0002 / 0.0052 / 0.0391 / 0.1169 | 103 / 61 |
| strong | LS       | 0.0008 / 0.0216 / 0.0657 / 0.0782 | 79 / 43 |

(mse at sweep 0 ≈ 5e-5 both regimes; all 20 runs monotone in their own
objective.) **Interpretation:** weighting removes the *heteroscedastic*
component of the D1 pathology (the early drift slows 4–30×) but not
its core: with L = 200 tokens per document, the per-document scores
are incidental parameters, and the joint (Z, Φ) likelihood descent is
Neyman–Scott-inconsistent — each Φ-update absorbs the per-document
ML noise, each z-update re-fits to the distorted Φ, and B̂ drifts
without bound. The deviance argmin is not a usable estimator either.

## P2a — the twist: deviance is not even the better step criterion

Figure `results/p2_k_curves.png`. From the *same oracle-GL start* on
the same data, the first deviance sweep already multiplies B̂'s error
by ~60× relative to the first LS sweep (weak k = 1: 0.0041 vs 0.00006),
and it deteriorates from there, while LS holds 0.0003 (weak) and
*improves* to 0.0004 (strong, k = 100). Mechanism (consistent with
D2): given a slightly-blended Φ, the per-document deviance MLE is
strongly inflated, and damped Fisher scoring goes there aggressively;
the LS z-step's "mis-weighting" acts as implicit damping that keeps the
iterate near the (information-bearing) start. Two stability guards
were required even to run the deviance blocks from imperfect starts —
a Φ-first EM ordering and a trust-region cap ‖δz‖ ≤ 1 (without the cap,
documents jump to softmax saturation where gradients die; measured
mse 2–9). **For the k-step estimator, the step criterion stays LS**
(possibly with a noise-matched structural shrinkage on z — the
correctly-scaled version of the λ-penalty — as the next refinement).

## P3 — feasible chain: the M-regime result

| cell | anchor TV | polished TV | chain mse (paper) | norm ratio | vs STM |
|---|---|---|---|---|---|
| weak, M = 1000   | 0.342 | 0.247 | 0.0582 | 3.80 | 0.0084 (STM wins) |
| strong, M = 1000 | 0.210 | 0.143 | 0.1260 | 2.16 | 0.0220 (STM wins) |
| strong, M = 5000 | 0.079 | 0.062 | **0.0015** | **1.07** | 0.0136 (**chain wins 9×**) |

The deviance polish is the one component that unambiguously helps
(30% TV reduction at every cell — the "likelihood-based anchor
recovery" flagged in the feasibility report works). The chain's
binding constraint remains anchor/Φ̂ quality: at M = 1000 the residual
Φ̂ flatness still inflates ẑ (norm ratios 2.2–3.8) and the chain loses
to STM; by M = 5000 the anchors are clean (TV 0.06), the inflation is
gone (1.07), and the fully feasible chain beats STM 9× — confirming
the F4(ii) extrapolation (slope −0.83) and locating the crossover
between M = 1000 and M = 5000 for this design. Coverage remains low
(0.61 at M = 5000): the 1/L bias question is untouched by this round
(D3's joint correction is still open).

## End-state decision (the closing paragraph)

**The data select End-state A, in a sharper form than pre-registered.**
The k-step theorem around an oriented pilot is the definitive framing:
P1 shows the deviance argmin inherits the D1 drift (Neyman–Scott, not
heteroscedasticity — gate ratios 42–7735, tilt-dominated displacement),
so no criterion in this family supports an argmin/M-estimator theorem
on the anchored slice at fixed L; and P2a shows the deviance is not
even the preferred *step* criterion — its exact per-document fitting
amplifies the incidental-parameter noise that the LS step implicitly
damps (first-sweep error ×60, monotone deterioration thereafter). The
theorem should therefore be stated for the **LS k-step refinement with
the B-stationarity rule (cap 100) around an oriented pilot**, with
deviance retained in exactly one place: the **EM polish of the anchored
Φ̂** (TV −30% at every cell), which is a Φ-only step and immune to the
incidental-parameter mechanism. The operating-regime statement writes
itself from P3: with feasible anchors the chain is STM-dominated at
M = 1000 but beats STM 9× at M = 5000 with norm ratio 1.07 — the
oracle gap is a finite-M anchor-quality gap, closed by corpus size,
not a structural impossibility. Open for the paper: the joint 1/L
bias correction (D3) and the noise-matched structural shrinkage of the
z-step.

## Deviations

1. Two stability guards added to the deviance blocks after measured
   failures (both logged with before/after numbers, both applied
   uniformly to all runs): Φ-first EM ordering in `dv_refine`
   (clipped starts otherwise produce p ≈ 1e-8 on counted cells) and a
   trust-region cap ‖δz‖ ≤ 1 on the Fisher step (prevents saturation
   stranding; the truth-start trajectories are unaffected by either).
2. Gate P1's literal denominator (mse at sweep 0, i.e. of B(Z_true) ≈
   5e-5) makes the 2× gate unpassable even for an estimator at the
   oracle floor (0.0003 ≈ 6× mse0); the verdict above therefore also
   reports the brief's secondary reading (distance to the oracle
   reference) — the conclusion is the same under both (12–200× above).
3. The LS truth-start probe was re-run (rather than reused from F3)
   to add the gauge-drift decomposition with identical recording;
   mse paths agree with the F3 pathology numbers on the shared seeds.
4. LS oracle-start k-curve references reused from feasibility
   `f3_results.rds` (identical seeds/data); STM references from audit
   `b2_results.rds` / `b2_m5000_results.rds` (identical seeds/data).
5. P2b skipped per protocol (gate P1 failed).
6. Windows: PSOCK workers (10), as in all previous rounds.

## Runtime and seeds

- Wall-clock (10 workers): unit tests 0.2 s; P1 76 s; P2a 108 s; P3
  218 s; P3 crossover 237 s — **total ≈ 10.7 min** (budget 1 h);
  development-time stabilization probes ≈ 5 min (scratchpad).
- Seeds: P1/P2a/P3 `90000 + regime·1000 + rep` (weak = 1, strong = 2;
  P1 reps 1–5 = the F3/D1 probe data; P3 M = 5000 uses regime strong,
  reps 1–10, paired with audit B2+); unit tests 66001. Thinning: not
  used this round (no jackknife).
- Monotonicity: 100% of runs monotone in their own criterion (the P1
  truth-start runs by design — that is the finding).
- Files: `01_deviance_blocks.R` (blocks + unit tests), `02_run.R`
  (`--quick`, `--ncores`), `03_report.R`; `results/*.rds`,
  `results/summary_dev.csv`, `results/tables_dev.md`, figures
  `p1_pathology_overlay.png`, `p2_k_curves.png` (`p2_coverage.png`
  not produced — P2b skipped), log `results/full_run.log`.
- No package function or pre-existing script modified.
