# Empirical verification report — CoDa-favorable regime + CRC microbiome application

**Package:** `sgscatm` (ILR-Spectral-GSCA structural topic model)
**Run date:** 2026-07-06 / 07 · **R:** 4.5.1 · loaded via `devtools::load_all(".")`
**Pre-registration:** [output/preregistration.md](preregistration.md) (written before any run; gates fixed there)
**Scripts:** `replication/microbiome/` (runnable end-to-end, numbered 00–08)

This run answers the reviewer objection that the delocalization condition (Assumption 5)
restricts the operating regime. It (i) shows near-nominal closed-form inference in the favorable
regime via a microbiome-calibrated simulation and traces the delocalization boundary, and (ii)
applies the method to a real colorectal-cancer (CRC) gut-metagenomic cohort. Every pre-registered
gate is reported below, **including the ones that did not pass as stated**.

---

## Executive summary

| Gate | What it tests | Verdict | Headline number |
|------|----------------|---------|-----------------|
| **G1** | Delocalization on real data | **Partial** | `r_fit` ≈ 6.5 (genus) / 7.2 (species), stable across all prevalence filters; predicted genus≪species gap **not** observed |
| **G2** | Coverage in favorable regime | **Pass** | coverage 0.94–0.97 across 23 cells; SE/SD ≈ 1.0; never pinned at 1.000 |
| **G3** | Delocalization crossover | **Pass (bias axis)** | recovery RMSE rises ×150 (0.005→0.79) as `r_true` grows 7k→644k; SE stays calibrated |
| **G4** | SE calibration on real data *(critical)* | **Pass** | analytical/bootstrap SE ratio median **1.15**, 81% within ±25%; legacy SE **9.3× too wide** |
| **G5** | Concordance on real data | **Pass (directional)** | **6/6** known CRC genera load CRC-ward (sign-test p=0.016); agrees with PERMANOVA on disease term |
| **G6** | Speed on real data | **Pass** | sgscatm 0.02 s vs PERMANOVA 0.76 s (38×), STM 1.01 s (50×), ALDEx2 3.45 s (172×) |

**Bottom line.** The corrected closed-form standard error (`sgscatm_vcov`, `R/vcov.R`) is
well-calibrated both in simulation and on real CRC data — the central de-risking result. The
previously-diagnosed over-wide analytical SE (`ilr_se_analytical`) reproduces here as a 9.3× (real
data) to ~290× (simulation) inflation, confirming it is the scale-collapse failure the prior round
found. The method sits exactly where the paper positions it: a fast, coefficient-level compositional
regression with valid closed-form inference between PERMANOVA (omnibus) and ALDEx2/ANCOM-BC
(per-taxon), and it does **not** replace per-taxon testing (see G5).

### Deviations from the plan (logged)

1. **Analytical SE implementation.** The plan named `ilr_se_analytical()`. That closed form is the
   miscalibrated one the prior round flagged; the repository already contains the corrected
   `sgscatm_vcov()` (`R/vcov.R`, standardized scale + eigenvalue-fluctuation term + `1/M` factor +
   native Procrustes/identification support). We use `sgscatm_vcov()` as the primary SE and report
   `ilr_se_analytical()` alongside as the failure baseline. **No package `R/` code was modified.**
2. **Real dataset.** `curatedMetagenomicData` **cannot be installed** on this system: on
   Bioconductor 3.22 its load chain (via `mia`) references `rbiom::unifrac`, which the installed
   `rbiom` 3.1.0 no longer exports (`object 'unifrac' is not exported from 'namespace:rbiom'`). This
   also collaterally blocks **ANCOM-BC** (needs `mia`). We substitute the **Zeller 2014** French CRC
   metagenomic cohort bundled in **SIAMCAT** — one of the same cohorts cMD pools, and the canonical
   dataset that established *Fusobacterium/Peptostreptococcus* CRC enrichment. Consequence:
   single-cohort (no study/country meta-adjustment), and ANCOM-BC is reported as blocked. This is
   disclosed, not silently substituted (the pre-registration reserved this contingency).
3. **ALDEx2 mode.** `aldex.glm` hits an ALDEx2/R-4.5 S4-class bug ("condition has length > 1"); we
   run ALDEx2 in its canonical two-group (CRC vs control) mode instead. This drops covariate
   adjustment for ALDEx2 but keeps it a fair competitor for timing (G6) and taxa concordance (G5).

---

## Phase 1 — Microbiome-calibrated simulation

Reuses `replication/simulation/sim_dgp.R` + `sim_utils.R`. Dimensions swept:
N(taxa)∈{200,500}, M(samples)∈{500,2000}, K∈{5,8}, P=4, NegBin library sizes (mean 2e4), peaked
taxon profiles (`alpha_beta=0.05`). 23 cells, 50–60 replicates each.

**Estimand handling (why the metrics are scale-free).** `B_z` is identified only up to
rotation/sign, and the estimator standardizes ILR scores to unit variance, so `B_z` is **not** on
the raw generative `Bz0` scale (a naive comparison gives RMSE ≈ √E[Bz0²] because the estimate is
~0 on that scale). We therefore validate the SE against the **empirical sampling SD across
replicates** and report **coverage of the generalized-Procrustes across-replicate mean** — both
scale-free — and report recovery of `Bz0` separately.

### G2 — Coverage (PASS)

Across all 23 cells: coverage ∈ [0.944, 0.966], median SE/analytical-vs-empirical-SD ratio ∈
[0.96, 1.11], delocalization ratio `r_fit` ∈ [4.2, 7.5]. Coverage is **never** pinned at 1.000, so
this is genuine calibration, not the over-wide failure signature. A 60-replicate probe at the
favorable cell (M=2000,N=200,K=5,b_max=0.25) gave median SE/SD = 0.89 and coverage 0.926.

Figures: [phase1_coverage_vs_bmax.pdf](figures/phase1_coverage_vs_bmax.pdf),
[phase1_calibration_vs_bmax.pdf](figures/phase1_calibration_vs_bmax.pdf).
Tables: [phase1_bmax_sweep.tex](tables/phase1_bmax_sweep.tex),
[phase1_dim_robustness.tex](tables/phase1_dim_robustness.tex),
[phase1_all_cells.csv](tables/phase1_all_cells.csv).

A document-bootstrap cross-check inside the simulation (B=80) put the analytical SE within **2%**
of the bootstrap SE (ratio 1.02 at b_max=0.25; 0.999 at b_max=1.5), pre-validating the G4 machinery.

### G3 — Delocalization crossover (PASS on the bias axis)

The main sweep measured `r` on the *standardized* eigenvector scores, which is scale-invariant and
stayed flat (~6) — so a dedicated stress test (`03_phase1b_crossover.R`) drives covariate strength
with small residual noise and measures delocalization on the **natural** ILR scale:

| b_max | `r_true` (natural) | max‖z‖ | coverage | SE/SD | recovery RMSE |
|-------|-------|--------|----------|-------|---------------|
| 0.5 | 7,137 | 3.76 | 0.954 | 1.02 | 0.005 |
| 1.0 | 25,359 | 7.09 | 0.954 | 1.01 | 0.067 |
| 1.5 | 64,727 | 11.33 | 0.927 | 0.93 | 0.167 |
| 2.0 | 111,429 | 14.87 | 0.960 | 1.05 | 0.180 |
| 3.0 | 318,105 | 25.12 | 0.951 | 0.97 | 0.411 |
| 4.0 | 376,514 | 27.36 | 0.963 | 1.04 | 0.629 |
| 5.0 | 643,539 | 35.66 | 0.951 | 1.01 | 0.787 |

**Finding (precise).** As compositions leave the centroid, recovery of the true generative `Bz0`
degrades sharply — RMSE rises ~150× — which is the Proposition 16 crossover on the **bias**
axis. Simultaneously, **coverage and SE/SD stay near-nominal**: the closed-form SE remains a valid
estimate of the estimator's own sampling variability even deep in the delocalized regime. The
crossover therefore manifests as growing *bias* of the point estimate for the natural-scale
coefficient, not as *variance* miscalibration. This is a stronger statement about the inference than
the gate anticipated and is reported as such. Figure:
[phase1b_crossover.pdf](figures/phase1b_crossover.pdf).

---

## Phase 2 — CRC application (Zeller 2014 via SIAMCAT)

136 complete-case samples (5 dropped for missing BMI), 720 genera / 1742 species; genus table at
prevalence ≥10% has **N=108** genera (sparsity 0.45). Covariates: study_condition (CRC=1), age
(std), BMI (std), gender (M=1). All six known CRC-enriched genera retained.

### G1 — Delocalization on real data (PARTIAL)

Fit-based ratio `r = M·max‖z_i‖²/(K−1)` (as pre-registered), K=5:

| level | prev | N | `r_fit` |
|-------|------|---|---------|
| genus | 10% | 108 | 6.55 |
| genus | 20% | 87 | 6.54 |
| genus | 30% | 74 | 6.54 |
| species | 10% | 284 | 7.20 |
| species | 20% | 217 | 7.19 |
| genus (unfiltered) | — | 720 | 6.55 |
| species (unfiltered) | — | 1742 | 7.21 |

**Finding.** The favorable-regime half of the prediction **holds robustly**: `r_fit` is small
(order 6–7) and essentially constant across taxonomic rank, prevalence filter, and even the
unfiltered tables. The predicted *degradation at species* did **not** materialize — species is only
~10% higher than genus. The leading covariance directions stay delocalized because the disease
signal is spread across many samples, so no single document localizes the eigenvector. This
strengthens rather than weakens the method's applicability to real compositions. (The alternative
CLR-extremity measure `r_clr` is strongly pseudocount/N-dependent and is reported in
[phase2_delocalization.tex](tables/phase2_delocalization.tex) only as a sensitivity note; it is not
a reliable delocalization diagnostic.)

### G4 — SE calibration on real data (PASS — critical gate)

Document bootstrap (B=300), each replicate Procrustes + sign-aligned to the point estimate, vs the
analytical `sgscatm_vcov` SE:

- **median analytical/bootstrap ratio = 1.147**; **81%** of the 16 entries within ±25%.
- Per-covariate ratios 1.0–1.4 (analytical SE marginally conservative).
- **Legacy `ilr_se_analytical` median SE = 1.086 vs corrected 0.116 → 9.3× too wide** — the
  scale-collapse failure, reproduced on real data.

Passes both pre-registered thresholds (median in [0.8,1.25]; ≥60% within ±25%). Figure:
[phase2_se_calibration.pdf](figures/phase2_se_calibration.pdf), table:
[phase2_se_calibration.tex](tables/phase2_se_calibration.tex).

### G5 — Concordance on real data (PASS — directional)

**(a) Covariate significance vs PERMANOVA** (`vegan::adonis2`, Aitchison distance, marginal):

| covariate | sgscatm joint-Wald p | PERMANOVA p | agree? |
|-----------|------|------|--------|
| study_condition (CRC) | **0.0065** | **0.001** | ✔ both significant |
| age | 0.185 | 0.230 | ✔ both n.s. |
| BMI | 0.206 | 0.011 | ✗ disagree (PERMANOVA sig) |
| gender | 0.117 | 0.149 | ✔ both n.s. |

Agreement on the disease term and on 3/4 covariates.

**(b) Known CRC taxa.** The disease (CRC) direction in score space, mapped to genus loadings via
`Psi`: **all 6/6 known CRC-enriched genera** (*Fusobacterium, Peptostreptococcus, Parvimonas,
Gemella, Porphyromonas, Solobacterium*) load in the CRC-ward (positive) direction — sign-test
**p = 0.016**. ALDEx2 (two-group) independently flags Peptostreptococcus, Fusobacterium, Parvimonas,
Porphyromonas among its top genera; sgscatm's top-25 by signal-to-noise overlaps ALDEx2's top-25 by
10 genera.

**Honest caveat (consistent with the method's stated role).** The known biomarkers are
low-abundance, so although their *direction* is correctly recovered, their loading *magnitudes* are
small and they do **not** top the |loading| ranking (top-25 overlap with the known set is not
significant, hypergeometric p=0.42). sgscatm's dominant loadings are high-abundance genera
(*Bacteroides* +5.0, *Ruminococcus* −3.8, *Eubacterium* −2.3) that carry the most compositional
variance. This is exactly the pre-registered positioning: a **joint, coefficient-level** method that
orients the whole community shift, **not** a replacement for per-taxon differential-abundance
ranking. Figure: [phase2_known_taxa.pdf](figures/phase2_known_taxa.pdf), forest of `Bz`:
[phase2_Bz_forest.pdf](figures/phase2_Bz_forest.pdf).

### G6 — Speed on real data (PASS)

| method | wall-clock | speedup |
|--------|-----------|---------|
| **sgscatm** (fit + analytical SE) | **0.02 s** | — |
| PERMANOVA (adonis2, 999 perm) | 0.76 s | 38× |
| STM (variational EM, K=5) | 1.01 s | 50× |
| ALDEx2 (128 MC, two-group) | 3.45 s | 172× |
| ANCOM-BC | — | blocked (mia/rbiom conflict) |

sgscatm is 1.6–2.2 orders of magnitude faster; its optional bootstrap SE (B=300) adds 2.1 s and is
only needed for validation, since the analytical SE is calibrated (G4). Figure:
[phase2_timings.pdf](figures/phase2_timings.pdf), table:
[phase2_timings.tex](tables/phase2_timings.tex).

---

## Discussion of the non-clean gates

- **G1 partial** — The "degrades at species" sub-prediction did not hold on this cohort; the method
  stays in the favorable regime at both ranks. This is favorable evidence for the method but means
  the delocalization boundary is not exercised by real, prevalence-filtered genus/species tables.
  Phase 1b shows *where* it would break (natural-scale `r_true` in the tens of thousands, i.e.
  extreme near-degenerate compositions), which real filtered data do not reach.
- **G3 nuance** — The crossover is real but shows up as point-estimate bias, not SE miscalibration.
  If a downstream user cares about recovering the natural-scale coefficient (not just the
  standardized direction), the delocalized regime is where bias, not coverage, is the risk.
- **G5 caveat** — Directional concordance is significant; rank concordance with per-taxon tests is
  not, by construction. Reported plainly so the method is not oversold as a biomarker-ranking tool.
- **ANCOM-BC / cMD** — Blocked by a genuine upstream dependency break (`rbiom` API change on Bioc
  3.22), not by the method. Reproducing on a system with a compatible `rbiom`/`mia`/cMD stack would
  restore the pooled multi-cohort design and ANCOM-BC comparison.

## Reproduction

```
Rscript replication/microbiome/00_install_deps.R      # vegan, ALDEx2, ANCOMBC, stm, ...
Rscript replication/microbiome/00c_install_siamcat.R  # Zeller CRC data source
Rscript replication/microbiome/01_phase1_simulation.R # G2 + calibration sweep
Rscript replication/microbiome/03_phase1b_crossover.R # G3 crossover
Rscript replication/microbiome/02_phase1_figures.R
Rscript replication/microbiome/05_phase2_prep.R       # build genus/species compositions
Rscript replication/microbiome/06_phase2_fit.R        # G1, G4, sgscatm CIs
Rscript replication/microbiome/07_phase2_competitors.R# PERMANOVA/STM/ANCOMBC + G5a
Rscript replication/microbiome/07b_competitors_fix.R  # ALDEx2 two-group, STM fix, G5b
Rscript replication/microbiome/08_phase2_figures.R    # figures + tables
```
