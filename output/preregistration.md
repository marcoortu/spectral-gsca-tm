# Pre-registration — CoDa-favorable regime + microbiome application (spectral GSCA)

**Written:** 2026-07-06, before executing any simulation or data analysis.
**Package:** `sgscatm` (ILR-Spectral-GSCA structural topic model), loaded via `devtools::load_all(".")`.
**Purpose:** Answer the reviewer objection that the delocalization condition (Assumption 5)
restricts the operating regime, by (i) demonstrating near-nominal behaviour in the
favorable regime via a microbiome-calibrated simulation, and (ii) applying the method to a
real colorectal-cancer (CRC) microbiome meta-analysis, with honest reporting of every gate.

This is a confirmatory run. Predictions and pass/fail gates are fixed below **before** any
result is seen. Every gate is reported as a finding, failures included. A gate is not
retro-fitted to the data.

---

## Estimand and identification (fixed conventions)

- The target is the ILR path-coefficient matrix **B_z** (P × (K−1)), regressing ILR topic
  scores on covariates C.
- **B_z is identified only up to an orthogonal rotation / sign of the score basis**
  (spectral solution). Therefore *every* comparison of an estimate to a reference — truth in
  simulation, point estimate in bootstrap — is made **only after Procrustes alignment**
  (`procrustes_align`) of the estimate to that reference. Skipping alignment inflates any
  spread/SE spuriously; this is treated as a hard methodological requirement, not an option.
- Coverage and RMSE are computed **entry-wise on the aligned B_z**.

## SE implementation — deviation from the plan, logged

The plan's working assumptions name `ilr_se_analytical()` as the closed-form SE. A prior
verification round diagnosed that closed form as **miscalibrated** (over-wide, from a
unit-norm eigenvector-score scale collapse). The repository now contains a corrected
implementation, `sgscatm_vcov()` in `R/vcov.R` (untracked at run start), which:
(a) works on the standardized score scale `Z̃ = √M Z*` where `M⁻¹ Z̃ᵀZ̃ = I`;
(b) adds the eigenvalue-fluctuation influence term (missing before);
(c) carries the final `1/M` so SE = O(M^{−1/2}); and
(d) natively supports Procrustes/varimax rotation (`rotation=`) and projects out the
    rotational tangent (`identified=TRUE`) — exactly what aligned coverage/calibration need.

**Decision:** use `sgscatm_vcov()` as the **primary analytical SE** throughout, passing the
Procrustes rotation R so the SE is expressed in the aligned basis. `ilr_se_analytical()` is
also computed and reported side-by-side, so the reader sees whether the old closed form is
the failure signature G2 warns about. `ilr_se()` (nonparametric document bootstrap) is the
reference for G4. This deviation and its reasoning are logged here and repeated in the report.

---

## Gates (pass/fail fixed now)

**G1 — Delocalization (real data).** Delocalization ratio
`r = M · max_i ||z_i||² / (K−1)`. Predict: `r` is small (order ~1–10) and stable across
prevalence filters at **genus** level; it **degrades** (grows and/or destabilises) at
**species** level (sparser, more extreme compositions). PASS = genus `r` small and stable AND
species `r` materially larger. This is a directional prediction; we report `r` at both ranks.

**G2 — Coverage (simulation, truth known).** Analytical 95% CI for B_z reaches ≈ nominal
(target 0.90–0.97 mean entry-wise coverage) in the moderate-composition / large-library
regime. **Explicit failure signature:** coverage ≈ 1.000 on all entries = over-wide SE
(scale collapse), scored as FAIL, not success. PASS = coverage in [0.90, 0.97] and not pinned
at 1.000.

**G3 — Crossover (simulation).** As covariate strength `b_max` grows and `||z_i||` leaves the
centroid, coverage drops below nominal and RMSE(B_z) rises — a monotone degradation matching
Proposition 16. PASS = coverage decreasing and RMSE increasing across the `b_max` sweep, with
the delocalization ratio `r` rising in step.

**G4 — SE calibration (real data), CRITICAL.** Analytical SE matches the document-bootstrap
reference (≥200 resamples, each Procrustes + sign-aligned to the point estimate) within ~25%
on **most** entries (target: median |log ratio| within log(1.25); ≥60% of entries within
±25%). PASS = median analytical/bootstrap SE ratio in [0.8, 1.25] AND ≥60% of entries within
±25%. Reported both for `sgscatm_vcov` (primary) and `ilr_se_analytical` (legacy).

**G5 — Concordance (real data).** (a) Covariate directions the method flags significant agree
with PERMANOVA omnibus significance for `study_condition`. (b) Genera loading most strongly on
the disease-covarying ILR direction overlap known CRC-enriched taxa
(*Fusobacterium*, *Peptostreptococcus*, *Parvimonas*, *Gemella*, *Porphyromonas*,
*Solobacterium*) above chance (hypergeometric p < 0.05 for the overlap of the top-loading set
with the known set). PASS = both (a) and (b).

**G6 — Speed (real data).** sgscatm wall-clock is 1–3 orders of magnitude below ALDEx2 (Monte
Carlo), STM (variational EM), and PERMANOVA (permutation) on the same data. PASS = sgscatm
fit+SE time ≤ 1/10 of each competitor's time.

---

## Phase 1 — Microbiome-calibrated simulation (run first; de-risks Phase 2)

Reuse `replication/simulation/sim_dgp.R` + `sim_utils.R` (do not reimplement). Dimensions:
- N (taxa) ∈ {200, 500}; M (samples/documents) ∈ {500, 2000}; K ∈ {5, 8}; P = 4
  (2 continuous ~N(0,1), 2 binary via thresholded latent).
- Library sizes ~ NegBin(mean 1e4–1e5) via `sim_dgp_variable_length` length function.
- `alpha_beta` set to a microbiome-like sparse value (0.05) so topic-word (taxon) profiles are
  peaked, mimicking taxonomic profiles.
- Covariate-strength sweep `b_max` ∈ {0.10, 0.25, 0.50, 0.75, 1.00, 1.50} tracing the
  delocalization boundary.
- ≥ 50 replicates per cell for coverage/RMSE (analytical SE, fast). A reduced set of cells
  gets a document-bootstrap SE (B = 100) on ≥ 20 replicates for the bootstrap-vs-analytical
  ratio, to keep compute tractable — this reduction is disclosed, not hidden.

Per cell record: RMSE(B_z) after Procrustes, analytical 95% CI coverage, bootstrap/analytical
SE ratio, mean delocalization ratio `r`. Outputs: coverage-vs-`b_max` and RMSE-vs-`b_max`
plots (`output/figures/`), tables (`output/tables/`). Gates: **G2, G3**, plus the calibration
sanity check.

## Phase 2 — CRC real-data application (primary new deliverable)

Data: `curatedMetagenomicData`, pooled CRC cohorts. Select via `sampleMetadata` studies with
`study_condition ∈ {CRC, control}` and non-missing `age`, `BMI`, `gender` (Wirbel/Thomas-style
meta-analysis set). Pull relative abundance, aggregate to **genus**, prevalence-filter
(present in ≥ 10% of samples), renormalize to compositions.
Covariates: `study_condition` (binary), `age` (standardized), `BMI` (standardized), `gender`
(binary); adjust for study/country.
Steps: (1) report M, N(genera), library-size and sparsity, compute `r` at genus AND species →
**G1**; (2) fit sgscatm (K by scree + small λ grid via the existing calibrate routine), report
B_z with analytical CIs and wall-clock; (3) competitors on the same data/covariates — PERMANOVA
(`vegan::adonis2`, Aitchison distance), ALDEx2 (`aldex.glm`, 128 MC), STM
(samples=documents, genera=words), ANCOM-BC if straightforward — record wall-clock →
**G6**; (4) document-bootstrap (≥ 200) refit, Procrustes + sign-align each to the point
estimate, empirical vs analytical SE → **G4**; (5) map ILR directions to genus loadings,
compare disease-associated taxa vs PERMANOVA/ALDEx2 and known CRC taxa → **G5**.

**Contingency (disclosed now):** if the `curatedMetagenomicData`/ExperimentHub install or data
pull is infeasible in this environment (offline/blocked/timeout), that is reported as a gate
that could not be executed, with the exact failure, rather than substituted silently. Phase 1
stands on its own as the de-risking evidence.

## Phase 3 — Soil + pH robustness (stretch)

Only if Phases 1–2 are clean. Large soil 16S dataset with continuous pH covariate; repeat the
delocalization check and recover B_z(pH) with analytical CI.

---

## Principles

Pre-register before running. Report failed gates as findings. Keep bias-correction and
identification concerns separate. Minimal edits to package `R/` code (log any). Position the
method as joint, coefficient-level compositional regression with closed-form inference,
occupying the gap between PERMANOVA (omnibus) and ALDEx2/ANCOM-BC (per-taxon) — **not** as a
replacement for per-taxon differential-abundance testing.
