# Feasibility round — identified estimator, bias correction, k rule

Final pre-theory verification, building on `replication/basin_check/REPORT.md`
and `replication/audit_block1_stm/REPORT_AUDIT.md`. Reproduce via

```
Rscript replication/feasibility/04_run.R          # ~7 min on 10 workers
Rscript replication/feasibility/05_report.R
```

## Gate verdicts (one line each)

- **Gate F1 (best feasible within 2× oracle): FAIL** — best feasible
  (V4+jk) is 0.0140 / 0.0340 vs 2×oracle 0.0004 / 0.0042 (weak/strong):
  a 33×/8× shortfall.
- **Gate F2-i (jackknifed RMSE slope CI contains −0.5): FAIL** — slope
  +0.021 [−0.030, 0.072]; the jackknife *flattens* RMSE at ~0.024 and is
  worse than uncorrected for M ≥ 1000.
- **Gate F2-ii (entrywise coverage ≥ 0.88, non-degrading): FAIL** —
  jk coverage 0.940 / 0.878 / 0.707 / 0.507 at M = 500…4000 (degrading;
  uncorrected is better: 0.942 / 0.935 / 0.942 / 0.788).
- **Gate F2-iii (row-norm coverage improves over 0.70/0.43/0.37): FAIL**
  — jk row-norm coverage 0.853 / 0.687 / 0.193 / 0.080 (better than
  uncorrected only at M = 500).
- **Gate F4-i (feasible+jk ≤ STM, both regimes/metrics): FAIL** — V4+jk
  0.0140 / 0.0340 vs STM 0.0084 / 0.0220; STM beats every feasible
  variant (the oracle-start refined estimator still beats STM 10–40×).
- **Gate F4-ii (Block 1 coverage in [0.88, 0.98], slope CI ∋ −0.5):
  FAIL** — coverage 0.38–0.45; slope −0.828 [−0.950, −0.706] (steeper
  than −0.5 because anchor quality itself improves with M).

**All four pre-registered predictions fail.** The round's purpose was
to decide whether the theory gets written as planned; the verdict is
that it must be redirected. The failures are informative and coherent —
three structural discoveries below.

## The three structural discoveries

**D1. The exact frequency-LS criterion does not identify B at finite L —
constrained or not, from any start.** The definitive probe (F3
pathology, 5 reps/regime): simplex-constrained block descent started
**at the truth** decreases F monotonically (weak: 4.757→4.727; strong:
4.777→4.762; monotone 10/10) while mse(B̂z) rises **two orders of
magnitude** (weak 0.00013→0.0147; strong 0.00075→0.0776). The simplex
constraint on Φ removes the scale/gauge freedom but does not make the
criterion statistically identify B: unweighted LS on frequencies lets
small-cell multinomial noise (var ∝ p/L, heteroscedastic across four
orders) purchase objective decreases at the price of distorted (Z, Φ).
Every good result in all three rounds came from an oriented start plus
**early stopping**. Consequence: the theorem must be a **k-step (one-step
efficiency style) estimator theorem**, not an argmin/M-estimator
theorem; and the natural repair of the criterion itself is multinomial
weighting (deviance / quasi-likelihood), which is the first thing the
redesigned theory should examine.

**D2. Feasible identification via anchors transfers its Φ error into a
~1.5–3× scale bias on ẑ, which no downstream step removes.** The
anchor pipeline passes its unit test (TV = 0.126 ≤ 0.15 at
α_β = 0.05, M = 2000) but delivers TV ≈ 0.21–0.34 at the Block 3
working point (α_β = 0.1, M = 1000). A too-flat Φ̂ forces inflated
per-document ẑ (norm ratios 1.7–3.5 across V1–V4), the inflation
enters any orientation estimate built from those ẑ (V4's feasible GL
map has singular values ~2–3× the oracle's), and refinement cannot fix
it because scale lives in the criterion's flat directions (D1).
Measured chain: anchored Φ TV 0.21–0.34 → ẑ scale bias ×1.5–3 →
mse(B̂) floor 0.014–0.053, i.e. 30–100× the oracle reference. The F1d
boundary scan confirms the mechanism rather than the registered
prediction: as α_β rises to 1.0, true exclusivity falls (1.00→0.74)
and anchor TV worsens (0.23→0.39), yet V2's mse *improves*
(0.093→0.045) — the error tracks the *flatness mismatch* between Φ̂
and the true Φ, not anchor TV per se (P-F1's monotone-degradation
prediction fails in the inverted direction).

**D3. The split-document jackknife as designed is not a valid 1/L bias
correction for this estimator.** On the Block 1 grid it injects an
M-independent error (~0.024 RMSE) that dominates for M ≥ 1000, and on
the L-grid it *inflates* MSE 1.4× at L = 50 and L = 200 (ratio 0.71×)
and is neutral at L = 1000. The design held Φ̂ fixed at the full-data
refined value and refit only ẑ on each half: that assumed the
document-level ẑ|Φ bias is the dominant term and exactly 2× at L/2.
Neither holds — the relevant bias is a joint (Z, Φ) object (Φ̂ has
already absorbed the full-data ẑ biases through the Φ-steps), so
`2·B_full − (B_A+B_B)/2` cancels the wrong quantity and adds a
systematic offset. Where the error is dominated by the *anchored scale
bias* (feasible variants, F1c) the jackknife helps ~25–40% (V4
0.0512→0.0340, V2 0.0904→0.0533) — consistent with that bias being
roughly L-linear — but it never approaches a full correction. The
analytic second-order per-document correction (optional in the brief)
was not attempted; given D3 it becomes the *required* alternative, and
it must be derived for the joint estimator, not for ẑ|Φ.

## F3 — k curves and the recommended rule

| regime | k=3 | k=5 | k=10 | k=20 | k=50 | k=100 | rule stop (med) | rule MSE |
|---|---|---|---|---|---|---|---|---|
| weak   | 0.00017 | 0.00021 | 0.00024 | 0.00026 | 0.00026 | 0.00026 | 15 | 0.00026 |
| strong | 0.00216 | 0.00211 | 0.00174 | 0.00116 | 0.00056 | 0.00038 | 50 (cap) | 0.00056 |

Figure `results/f3_k_curves.png` (with the D1 pathology overlay).

- The two regimes behave oppositely. Weak: flat in k (the small early
  rise is the oracle-optimism of the GL start washing out, cf. the
  audit's "sweeps to beat STM = 0" artifact); anything k ≥ 3 is
  equivalent. Strong: mse still *improving* at k = 100 — the audit's
  k = 5 vs convergence gap (0.0021 vs 0.0004) is a **regime-wide slow
  transient, not a minority of slow replicates**: the median rule stop
  sits at the cap and the mse path decreases throughout in essentially
  every replicate. P-F3's "minority with identifiable signature" fails:
  the strongest correlate of sweeps-to-1.1×-final is ρ_GL with
  *negative* sign (−0.57; slower runs start closer), eigengap and
  saturation correlate at |r| ≤ 0.27.
- **Recommended rule (one line, used verbatim in F4):** *stop when
  max|ΔB|/rms(B) < 1e-3 on 2 consecutive sweeps, cap 100* — the cap at
  50 measurably truncates the strong regime (0.00056 vs 0.00038 at
  k = 100), so the cap moves to 100; with D1 in hand the rule's virtue
  is that it stops B-stationarity long before the criterion's bad
  region (which the oracle-start unconstrained path approaches only
  far beyond k = 100, if at all).

## F1c — feasible variants vs oracle (M = 1000, 20 reps, E2 seeds)

Full table in `results/tables_feas.md`; essentials (paper metric,
mean mse; permutation metric within 15% of paper everywhere):

| variant | weak | strong | norm ratio (w/s) |
|---|---|---|---|
| V1 anchored-only                | 0.0471 | 0.0749 | 3.5 / 1.8 |
| V2 anchored + refined           | 0.0287 | 0.0904 | 3.0 / 2.0 |
| V2 + jackknife                  | 0.0179 | 0.0533 | 2.5 / 1.8 |
| V3 pilot-NNLS-mapped + refined  | 0.0738 | 0.4318 | 4.2 / 3.2 |
| V3 + jackknife                  | 0.0565 | 0.3524 | 3.7 / 3.0 |
| V4 anchor-oriented pilot (added)| 0.0201 | 0.0512 | 2.6 / 1.7 |
| **V4 + jackknife (best feasible)** | **0.0140** | **0.0340** | 2.3 / 1.5 |
| oracle k = 5                    | 0.0002 | 0.0021 | 0.9 / 1.1 |
| oracle rule                     | 0.0003 | 0.0006 | 0.9 / 1.1 |
| STM (audit B3, same data)       | 0.0084 | 0.0220 | 2.1 / 1.5 |

- **Pilot-value measurement (explicit):** with anchors in place, the
  pilot buys a 22–36% mse reduction (V4 vs V2: 0.0179→0.0140 weak,
  0.0533→0.0340 strong, jk versions) at equal wall time (~4.4 s vs
  4.3 s; both hit the 50-sweep cap; V1's z-init uses 10 per-document GN
  iterations either way). The brief's specified pilot mapping (V3,
  NNLS through near-uniform reconstructed compositions) *destroys* the
  pilot's information (0.43 strong — worse than no pilot); the added
  V4 mapping (feasible 4×4 GL of the pilot scores onto the anchored
  scores) is the version that captures it. So: the pilot helps
  measurably but modestly in the feasible chain — its irreplaceable
  role remains inside the oracle frame, which no feasible device in
  this round reproduces.
- Note the norm-ratio column: every feasible variant is inflated
  (1.5–3.5×), STM's own 1.5–2.1× inflation (audit finding) puts it in
  the same family; the oracle variants sit at 0.9–1.1.

## F1d — identification boundary (strong, M = 1000, 10 reps/cell)

| α_β | true exclusivity | anchor TV | V2 mse (paper) |
|---|---|---|---|
| 0.05 | 1.000 | 0.232 | 0.0954 |
| 0.10 | 1.000 | 0.228 | 0.0933 |
| 0.30 | 0.964 | 0.252 | 0.0816 |
| 1.00 | 0.744 | 0.386 | 0.0447 |

Inverted relative to the registered expectation (see D2): B error
tracks the Φ̂-vs-Φ flatness mismatch, not exclusivity. Figure
`results/f1_alpha_boundary.png`.

## F2 — jackknife tables

Block 1 grid (oracle start, adaptive rule, 50 reps/M):

| M | RMSE unc | RMSE jk | cov unc | cov jk | cov jk+infl | rownorm unc | rownorm jk |
|---|---|---|---|---|---|---|---|
| 500  | 0.0272 | 0.0242 | 0.942 | 0.940 | 0.970 | 0.913 | 0.853 |
| 1000 | 0.0205 | 0.0237 | 0.935 | 0.878 | 0.938 | 0.900 | 0.687 |
| 2000 | 0.0139 | 0.0255 | 0.942 | 0.707 | 0.823 | 0.860 | 0.193 |
| 4000 | 0.0135 | 0.0242 | 0.788 | 0.507 | 0.595 | 0.640 | 0.080 |

Slopes: uncorrected **−0.369 [−0.422, −0.317]** over the full grid
(−0.48 over M = 500…2000; the L-floor bites at M = 4000, where the
RMSE flattens at 0.0135 and coverage drops to 0.788), jackknifed
**+0.021 [−0.030, 0.072]**. L-grid (weak, M = 1000): jk/unc mse ratios 1.41× at
L = 50, 1.41× at L = 200, 1.00× at L = 1000 — the correction *adds*
error exactly where it was supposed to help. Note also how much the
adaptive-rule fit improved the *uncorrected* estimator over the
audit's k = 5 (RMSE 0.036/0.030/0.023 → 0.027/0.021/0.014, coverage
0.87/0.77/0.71 → 0.94/0.94/0.94 for M ≤ 2000): a large share of what
the audit attributed to the L-floor was still optimization error at
k = 5; the genuine floor shows at M = 4000 (coverage 0.788). Figures
`results/f2_slope.png`, `results/f2_coverage.png`. The "inflated SE"
variant (sandwich + entrywise (B_A−B_B)²/4) recovers part of the
coverage but cannot repair a mis-centred correction; **the theory
should adopt neither** — see D3.

## F4 — dress rehearsal

(i) E2 seeds (from F1c, jackknifed, sandwich CIs): table above; STM
wins against all feasible variants; both metrics agree. (ii) Block 1
grid with the V4 chain (auto-selected): RMSE 0.329/0.218/0.124,
entrywise coverage 0.38–0.45, slope −0.828 [−0.950, −0.706] — anchor
quality improves with M faster than √M, another sign the binding
constraint is Φ̂ estimation error, not sampling noise in B.

## Deviations from the brief (all logged, none tuned post-gate)

1. **Anchor pipeline hardening** (within F1a's "simplified robust
   variant" allowance; all changes made to *pass the pre-registered
   unit test*, before any downstream experiment ran): candidate cutoff
   docfreq ≥ 5% of M (a ≥10-docs cutoff selects 11–16-document words
   whose profiles are noise); affine-span successive projection;
   rank-K eigenspace denoising of Qbar; **exact support-enumeration QP**
   for the simplex regressions (projected gradient silently stalls on
   the ill-conditioned anchor Gram matrix — C error 0.25); posterior
   threshold 0.10 (kills the noise-sprinkle on the ~85% true-zero Φ
   entries). Exclusivity re-selection round implemented; it re-picks
   the same anchors (no-op). Unit test after hardening: TV 0.126.
2. **Φ-step**: FISTA (300 iterations) replaces 10-step PGD, which froze
   the Φ block. Important negative control: FISTA did **not** change
   the constrained-run outcomes — D1 is a property of the criterion,
   not an optimization artifact.
3. **V4 added** (anchor-oriented pilot, unconstrained k-rule
   refinement) after V3-as-specified failed; V3 retained and reported.
   F4's init chooser picks among V2/V3/V4 by jackknifed mse (picked V4).
4. F2's full fit uses the adaptive rule rather than k = 5 ("standard
   fit" read as the final architecture; the k = 5 sensitivity is
   visible in the audit tables and F3 curves).
5. `fs_sandwich` duplicates 15 lines from the audit script (that file
   is a runner; sourcing it would launch the audit).
6. STM columns reused from audit B3 (identical seeds/data), not refit.
7. F3 carries the D1 pathology probe (constrained-from-truth paths,
   reps 1–5 per regime).
8. Analytic second-order bias correction: not attempted (see D3 — it
   is the recommended next step, but it is not "trivially within
   budget" once it must target the joint bias).
9. Windows: PSOCK workers (10), as in previous rounds.

## Runtime and seeds

- Wall-clock (10 workers): unit test 0.5 s; F3 42 s; F1c 65 s; F1d
  17 s; F2 grid 163 s; F2 L-grid 11 s; F4(ii) 99 s — **total ≈ 6.6 min**
  (budget 2 h); development-time diagnostics (anchor debugging, valley
  probes) ≈ 15 min extra, all preserved in the scratchpad scripts.
- Seeds: F3/F1c/F4(i) `90000 + regime·1000 + rep` (weak = 1,
  strong = 2; identical data to basin_check E2 / audit B3);
  F1d `30000 + α_index·1000 + rep`; F2 grid and F4(ii)
  `60000 + M_index·1000 + rep` (M = 500,1000,2000,4000 → 1..4; paired
  with audit A3 for M ≤ 2000); F2 L-grid `70000 + dl_index·1000 + rep`
  (paired with basin_check E5); thinning seed = replicate seed + 500;
  unit tests 55001 (anchors), 7101/7001 inherited from the audit.
- Monotonicity: all refinement runs monotone (flags recorded per
  replicate; the D1 pathology runs are monotone *in F* by design —
  that is the point).
- Files: `01_anchors.R`, `02_constrained_refine.R`, `03_jackknife.R`,
  `04_run.R` (`--quick`, `--ncores`), `05_report.R`; `results/*.rds`,
  `results/summary_feas.csv`, `results/tables_feas.md`, figures
  `f1_alpha_boundary.png`, `f2_slope.png`, `f2_coverage.png`,
  `f3_k_curves.png`, log `results/full_run.log`.
- No package function or pre-existing script modified.

## What the theory should now be

1. **A k-step theorem around an oriented pilot** (one-step efficiency
   style): conditions = pilot subspace consistency + an orientation
   oracle with o(1) error + k fixed or B-stationarity stopping. The
   argmin framing is dead (D1).
2. **Reweight the refinement criterion** (multinomial deviance /
   quasi-likelihood) and re-test D1: if the weighted criterion's
   descent from the truth no longer degrades B, the M-estimator
   framing can be revived with the identification devices of this
   round; this is the single highest-value experiment left.
3. **Orientation is the bottleneck for feasibility**: anchors at
   realistic separation transfer Φ̂ flatness into a multiplicative B
   bias (D2). Either strengthen the Φ̂ estimator (likelihood-based
   anchor recovery), assume stronger separation, or state the theorem
   with the orientation error entering B's bias linearly (and report
   the oracle gap as irreducible without it).
4. **Bias correction must target the joint (Z, Φ) bias** (D3): derive
   the analytic second-order correction for the k-step estimator, or
   use L-extrapolation across genuinely refit (not Φ̂-frozen) halves.
