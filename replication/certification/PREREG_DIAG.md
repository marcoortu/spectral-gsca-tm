# Pre-registration — Decisive diagnostic: is the chain's coverage failure curable (A') or structural (B')?

**Written before running.** Numbers map to A'/B' by the §5 rule below; nothing is tuned.
**Package:** `sgscatm`, `devtools::load_all(".")`. Script: `replication/certification/diag_bias.R`.

## Question

G2c showed chain + plain Lemma-17 sandwich under-covers in-regime, with two separable causes:
**(V)** omitted anchor variance (curable; jackknife inflation exists in `feasibility/03_jackknife.R`);
**(B)** a residual bias surviving even at oracle orientation (oracle coverage 0.44 at L=1e4,M=2000;
0.15 at M=5000 with SE/SD→0.91). The whole question is whether (B) is the **alternating/joint
refinement drift** (Prop-19, L-insensitive, curable by frozen-Φ early stop) or the **structural
1/L incidental bias** (Prop-21, ≈1e-4 at L=1e4, uncurable). The G2c oracle runs used
`refine="joint"`. Hypothesis: (B) is the joint drift, and frozen-Φ stopped early avoids it.

## Orientation modes (precise)

- **oracle** = freeze **true** Φ₀ (= `dat$Beta`), skip anchors; start = general-linear alignment of
  the pilot to the **true** centred scores. Removes (V) and δ_Φ; isolates the refinement.
- **feasible** = anchored Φ̂ + anchor GL alignment (`sgscatm_chain` stage 2).

## Design cell (in-regime)

`alpha_beta=0.05`, K=5, P=4 (2 continuous, 2 binary), `b_max=0.5`, `sigma_eps=0.3`, NegBin large L.
- **L ∈ {1e3, 1e4}** (both close √M/L; the L-contrast is the discriminator: structural bias ∝1/L,
  joint drift L-insensitive).
- **M ∈ {1000, 2000, 5000}** (SE∝M^{-1/2}; fixed bias becomes √M-amplified).
- **N=500** at oracle (Φ true, N only affects noise); **N=200** at feasible (clean anchors, small δ_Φ).
- R = 50 per cell (80 at the decisive oracle & feasible cells L=1e4,M=2000). Seeds base+regime·1000+rep.

## Parts

- **Part 1 (oracle Φ₀, N=500):** frozen_phi at k∈{1,2,3,5,10} (degradation curve) + joint@conv,
  over {L=1e3,1e4}×{M=1000,2000,5000}. Pick **k\*** = min-bias early stop at L=1e4,M=2000. Report
  bias-norm scaling in L and M for oracle-frozen@k\* vs oracle-joint@conv.
- **Part 2 (feasible, N=200):** feasible-frozen@k\* × {plain, jackknife-inflated SE} × M, L=1e4.
- **Part 3 (feasible, N=200, L=1e4):** feasible-frozen@k\* + jackknife → true coverage of Bz,0.

## Metrics (perm+sign align each B̂z to Bz,0 first; never Procrustes)

bias norm ‖mean(B̂z)−Bz,0‖_F; empirical SD; SE/SD per SE variant; RMSE²=bias²+var; bias²-share;
**true coverage** = mean 1[|B̂z−Bz,0|≤1.96 SE]; **bias-removed coverage** = mean 1[|B̂z−B̄z|≤1.96 SE]
(centre at replicate mean B̄z — measures variance calibration alone). k-curve (bias & RMSE vs k), k\*.
Jackknife: B_jk = 2B_full−(B_A+B_B)/2, var_add=(B_A−B_B)²/4 (split-document halves, Φ fixed).

## §5 Decision rule (frozen)

N1 = bias norm oracle-frozen@k\* (L=1e4,M=2000) + L/M scaling; N2 = true coverage oracle-frozen@k\*
(L=1e4); N3 = bias-removed coverage feasible-frozen@k\*+jackknife (L=1e4); N4 = true coverage
feasible-frozen@k\*+jackknife in-regime (L=1e4).

- **A' VIABLE** iff: N1 small AND not growing with M faster than SE shrinks (true & bias-removed
  coverage stay together as M grows) AND N2∈[.90,.97] AND N3∈[.90,.97] AND N4∈[.88,.97]; corroborated
  by oracle-frozen@k\* bias ∝1/L while oracle-joint bias L-insensitive. → deliverable inferential
  object = frozen-Φ@k\* + jackknife-inflated SE; proceed to full A' re-certification.
- **B' (structural wall)** iff: N2 stays low despite nominal bias-removed coverage, bias does not
  shrink at smaller k, and is L-insensitive at oracle. → chain for point estimation only; closed-form
  chain inference declared open; keep raw pilot + `sgscatm_vcov`.
- **Ambiguous / third finding:** oracle-frozen@k\* clean but N4 blocked by a *feasible* (anchor-bias)
  residual → report which term blocks and whether larger M/N closes it (feasibility-of-anchors, still
  favours A' with a stated regime). If even oracle-frozen@k\* under-covers with calibrated variance at
  L=1e4 (where 1/L bias≈1e-4) → suspected sandwich/alignment BUG: check perm+sign on a trivial case
  first, do not draw A'/B'.

## Protocol

Minimal edits (leave `sgscatm()`/`sgscatm_chain` core intact; log any). Reuse jackknife + DGP verbatim.
Fixed seeds, reproduce-from-scripts. A FAIL is the decisive finding selecting B'.

---

## Binding amendments (applied post-hoc, before the corrected feasible run)

- **A1 contamination guard:** point estimate for every N3/N4 coverage number is the chain's
  `B̂z` at frozen@k\*; the split-document jackknife contributes **only** `var_add=(B_A−B_B)²/4`.
  `B_jk = 2B_full−(B_A+B_B)/2` is **never** the point estimate (D3: split-document bias correction
  injects an M-independent error); reported only as a labelled **secondary** column, outside the verdict.
- **A2 identity check:** oracle freezes the topic–term matrix Φ₀ (K×N). Assert `dim(dat$Beta)==(K,N)`;
  abort otherwise.
- **A3 lead readout:** oracle-frozen@convergence (k=10) + PLAIN sandwich, L=1e4, across M — the
  estimator Lemma 17 is designed for; must cover nominally. If not (cov low, SE/SD≈1) → suspected-bug
  branch (re-check perm+sign + sandwich vs explicit ⊗ loop before any verdict).
- **A4 L-contrast:** report the frozen bias-vs-k curve for L=1e3 and L=1e4 side by side; flat at 1e4
  = residual frozen bias is the 1/L term, cured by depth (no early stop in-regime).
- **A5 N4 over M:** report feasible-frozen@k\*+`var_add` across M∈{1000,2000,5000}; over-coverage
  (>0.97) acceptable; only under-coverage is a failure. Verdict names the M at which coverage reaches
  nominal and the trend.
- **A6 explicit conclusion:** state whether, in-regime, early stop / jackknife are needed at all, or
  frozen-Φ@convergence + plain sandwich suffices.
