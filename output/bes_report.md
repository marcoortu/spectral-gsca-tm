# BES application — sgscatm vs STM in strong and weak covariate-signal regimes

**Data:** British Election Study Wave 25, "most important issue" (MII) open-text responses.
18,836 documents × 298 stemmed terms — **ultra-short texts** (mean 1.6 words/doc, median 2, max 12).
Covariates: `age_std`, `female`, `educ_std`, `lr_std` (left–right), `leave` (Brexit vote).
**Scripts:** `replication/application/bes_strong_weak.R`, `bes_base_and_figures.R`.

## What was tested

Whether **sgscatm recovers the covariate→content effect better and faster than STM**, and how
that advantage depends on signal strength, by **selecting documents into two regimes**:

- **Strong scenario** — politically committed respondents (`|lr_std| ≥ 1`; 6,690 docs), where
  Brexit stance strongly structures topical content.
- **Weak scenario** — centrists (`|lr_std| ≤ 0.3`; 4,325 docs), where Brexit stance weakly
  structures content.

**Ground truth is model-agnostic:** the observed Leave−Remain word-frequency shift on a large
held-out portion of each pool. Each method's *recovery* is the correlation of its **implied**
covariate→word shift (Δ topic-prevalence × topic-word matrix) with that observed shift. 25
resamples per scenario (train n=1,500; the rest held out), both models given the same five
covariates.

## Result — sgscatm dominates, and the gap widens in the weak regime

| scenario | sgscatm corr | STM corr | paired *p* | speedup |
|----------|:-----------:|:--------:|:----------:|:-------:|
| strong (|lr|≥1) | **0.962 (0.010)** | 0.877 (0.089) | 8.6e-5 | **14.6×** |
| weak (|lr|≤0.3) | **0.961 (0.011)** | 0.781 (0.256) | 1.8e-3 | **8.6×** |

(mean (sd) over 25 resamples). Recovery RMSE: sgscatm 0.004 vs STM 0.008 in both scenarios.
Figures: [bes_recovery_strong_weak.pdf](figures/bes_recovery_strong_weak.pdf),
[bes_timing_strong_weak.pdf](figures/bes_timing_strong_weak.pdf); table:
[bes_strong_weak.tex](tables/bes_strong_weak.tex).

Three findings:

1. **sgscatm is signal-strength invariant and stable.** Its recovery is 0.96 in *both* regimes
   with tiny across-resample SD (~0.01) — the global spectral solution has no initialization
   variance.
2. **STM degrades and destabilizes as signal weakens.** Mean recovery falls 0.877 → 0.781, and its
   across-resample SD *explodes* 0.089 → 0.256: on many weak-regime draws its variational EM fails
   to extract the faint covariate signal from ultra-short texts, on others it succeeds. The
   advantage of sgscatm is therefore both higher accuracy and far lower run-to-run variance.
3. **sgscatm is 9–15× faster per fit** on the subsamples, and **63× faster on the full corpus**
   (1.56 s incl. closed-form SEs vs 97.6 s for STM's EM). The recovery difference is statistically
   significant in both scenarios (paired *t*, p < 0.002).

## Base application (full 18,836 docs)

sgscatm, K=6, λ=3, fit + closed-form SE in **1.56 s**. Joint Wald tests: **every** covariate is
significant, `leave` overwhelmingly so (χ²=1991 on 5 df, p≈0; age χ²=696, lr χ²=422, female χ²=371,
educ χ²=246). Leave–Remain topic-prevalence directions are substantively sensible — the
immigration/"illegal/uncontrolled mass" topics (T1, T2) are Leave-leaning, while the
climate/NHS/poverty topics (T3, T4, T6) are Remain-leaning. Object saved to `output/bes_base.rds`.

## Interpretation

The BES MII corpus is an adversarial setting for topic models (1–2 words per document), and it is
exactly where the difference shows: STM's per-document variational inference has little to work with
per document and must lean on the prevalence prior, so it is noisy and unstable when the covariate
signal is weak; sgscatm aggregates the covariate–word covariance across the whole corpus in one
eigendecomposition, so it recovers the same effect at 0.96 correlation regardless of regime, with
valid closed-form inference, orders of magnitude faster. This is consistent with — and complements
— the microbiome result: the method's edge is **fast, stable, coefficient-level compositional
regression with closed-form inference**, most visible precisely where per-document methods struggle.

## Reproduction

```
Rscript replication/application/bes_strong_weak.R      # 25-resample strong/weak experiment
Rscript replication/application/bes_base_and_figures.R # full-data fit + figures + tables
```
