# Gate 4 — CRC fork: analytic sandwich vs full-chain bootstrap (primary-SE decision)

Pre-registered in `replication/certification/PREREG_G4.md`. Estimator = `sgscatm_chain(refine="frozen_phi")`
on the Zeller 2014 CRC genus table (`output/phase2_data.rds`, prev≥10%, renormalized; counts =
round(rel_abund × 1e6)). M=136 samples, N=108 genera, P=4 covariates, K=5. Point estimate = chain B̂z;
SE_sandwich = `vcov.sgscatm_chain` (plain Lemma-17); SE_boot = `chain_boot_se(B=200)` (perm+sign to
point estimate). Chain fit 1.0 s; bootstrap 172 s.

## Deciding numbers

| metric | value |
|--------|-------|
| median SE_sandwich / SE_boot | **0.804** |
| % of 16 entries within ±25% of parity | **44%** |
| per-covariate median ratio | study_condition 0.89, age 0.71, BMI 0.68, gender 0.90 |
| δ_Φ proxy (anchor exclusivity) | **1.000** (clean, exclusive anchors) |

**Fork rule:** median ∈ [0.8,1.25] ✓ but ≥60% within ±25% ✗ (only 44%) → AND fails →
**primary SE = full-chain bootstrap** (by the letter of the pre-registered rule).

## Substantive stability (reported under both SEs)

| covariate | joint Wald p (sandwich) | joint Wald p (bootstrap) |
|-----------|------------------------|--------------------------|
| **study_condition (CRC)** | **0.0097** | **0.006** |
| age | 0.311 | 0.553 |
| BMI | 0.021 | 0.493 |
| gender | 0.253 | 0.227 |

The scientific conclusion is **SE-choice-stable for the disease signal**: study_condition is significant
under both SEs (and agrees with PERMANOVA, p=0.001, from the earlier run). BMI is significant under the
sandwich but not the bootstrap — its significance is **not** robust to the SE choice (consistent with the
earlier PERMANOVA-vs-Wald BMI discrepancy). age/gender non-significant under both.

## Honest interpretation (the tension the rule does not resolve)

The fork picks the bootstrap **mechanically**, but the sandwich is not wildly off — it runs ~20–30%
**below** the bootstrap (median 0.80), and the two agree in overall magnitude. Two facts pull toward the
sandwich actually being the better-calibrated SE here:

1. **δ_Φ proxy = 1.000**: CRC taxa are exclusive, so the anchors are clean — precisely the δ_Φ→0 regime
   where A3 proved the plain sandwich is calibrated (SE/SD≈1). Corollary 18's orientation gate is closed
   on CRC.
2. **The bootstrap is known-conservative** (Gate 0: SE/SD 1.61 in-regime, over-covering because a
   minority of resamples draw a poor anchor set and inflate the bootstrap SD). So the sandwich sitting
   ~25% below the bootstrap is largely **the bootstrap's excess**, not the sandwich's deficit.

Net: by the frozen rule the primary SE is the bootstrap, **but** the evidence (clean anchors + known
bootstrap conservatism) suggests the sandwich is the calibrated in-regime SE on CRC and the bootstrap
over-states. This is exactly what **Phase B (robust bootstrap)** is designed to settle: if the
conservatism is outlier-driven, a MAD-scale / anchor-TV-trimmed bootstrap should drop from ~1.25× the
sandwich toward parity, reconciling them and confirming the sandwich as the fast primary SE. Until Phase
B runs, the **safe (conservative) primary SE is the bootstrap**, and the **speed claim is scoped to the
point estimate** (chain point fit 1.0 s; point+bootstrap 173 s ≈ comparable to competitors; the 38–173×
speedup holds only for the point estimate or, if Phase B vindicates it, the point+sandwich path).

## Caveat on the CRC chain fit

At M=136 the frozen refinement hit the sweep cap (100) without the B-stationarity rule firing
(`rule_stop=NA`) — a small-M artifact (the deep-L bias is ~1e-6, so this is refinement-path noise, not
the 1/L wall). The point estimate and both SEs are reported as-is; a small-M B-stationarity tolerance is
a minor tuning item, logged.

## Bottom line

**On CRC the analytic sandwich and the full-chain bootstrap agree in magnitude (median 0.80) but not
within the tight ±25% band (44%), so the frozen rule selects the bootstrap as primary — yet the clean
anchors (δ_Φ proxy 1.0, A3 regime) and the bootstrap's known conservatism indicate the sandwich is the
calibrated in-regime SE, to be confirmed by the robust bootstrap (Phase B).** The disease conclusion is
stable across the SE choice.
