# Paper edits implied by the chain certification

Edits are listed for the maintainer to apply to the `.tex` (not edited here). Ordered by
importance. Numbers marked "certified" are from this session's chain runs; those marked
"pending" require the full-scale re-run.

## 1. Central reconciliation — which estimator and which SE the paper actually uses

- **Section 6 (applications) was produced by the raw spectral pilot + `sgscatm_vcov`**, not by
  the three-stage chain of Algorithm 1 + the Lemma-17 sandwich. Figure 3's axis literally reads
  "analytical SE (sgscatm_vcov)". The prose (6.1.0.3) attributing the SE to "the analytical
  sandwich of Lemma 17" is therefore **incorrect as written**. Either (a) re-run Section 6 with
  `sgscatm_chain()` + `vcov.sgscatm_chain()` and update the numbers, or (b) state plainly that
  the in-regime applications use the standardized-scale pilot SE `sgscatm_vcov`, which coincides
  with the sandwich only when both rate gates close.

## 2. Frozen-Φ vs joint refinement (Section 3.7, Theorem 16, Table 3)

- The manuscript's "frozen Φ̂, no Φ update in the baseline" does **not** match the estimator that
  produces the reported feasible numbers, which re-estimates Φ each sweep (unconstrained; the V4
  variant in `replication/feasibility/`). **Certified evidence** (`CHANGES.md`): frozen-Φ is best
  after ~1 sweep and then monotonically degrades — a direct confirmation of Proposition 21 — while
  joint refinement improves monotonically. Recommend: present the deliverable estimator as the
  joint (Prop-19 alternating) refinement, and frame frozen-Φ as the theoretical idealization whose
  finite-L optimum is inconsistent (which the paper already argues). If frozen-Φ is kept as the
  inferential object, the Lemma-17 sandwich must be applied at the **early-stopped** iterate (k≈1).

## 3. The Lemma-17 sandwich under-covers for the chain — the coverage claim needs the
##    anchor-variance-aware SE (Section 4.6, Corollary 18, Section 5.3, Section 6.1.0.3)

- **Certified**: chain + plain Lemma-17 sandwich at L=1e4, M=2000, N=500 gives coverage 0.35 and
  SE ≈ 5× too small (G2c FAIL). The plain sandwich omits the anchor/orientation-stage variance,
  which dominates when √M·δ_Φ is not small. The paper should either (a) state that calibrated
  coverage requires the split-document **jackknife-inflated** variance (already in the replication
  code, `fs_coverage_entry(var_add=(B_A−B_B)²/4)`) or a full-chain bootstrap — not the plain
  sandwich — or (b) restrict the calibrated-coverage claim to the oracle-orientation / large-M
  regime where δ_Φ is negligible, and report the feasible-chain coverage separately. As written,
  "the analytical sandwich of Lemma 17 ... coverage nominal" overstates what the plain sandwich
  delivers for the feasible chain.

## 4. Table 1 principal-angle claim (Section 5.2)

- **Certified**: at the base design (N=500, L=200) the pilot recovers only the **leading**
  score directions to high accuracy (principal cosines 0.95, 0.83) while the trailing directions
  are noise-limited (0.71, 0.62), giving a largest principal angle ≈0.9–1.25 rad, **not** the
  <1e-5 stated in Table 1. The raw-regression M^{-1/2} collapse and the refined→oracle descent DO
  reproduce. Recommend qualifying the machine-precision claim to the leading subspace / a cleaner
  (larger-L or stronger-signal) design, or reporting the per-direction principal cosines.

## 5. Data-source and tool citations (Section 6.1, the `?` placeholders)

- CRC cohort: **Zeller et al. 2014**, accessed via the `SIAMCAT` example data
  (`feat.crc.zeller` / `meta.crc.zeller`); `curatedMetagenomicData` is blocked on the local
  Bioc 3.22 by an `rbiom` API break (`unifrac` no longer exported, via the `mia` load chain).
- Tools: PERMANOVA = `vegan::adonis2` (Oksanen et al.); ALDEx2 = Fernandes et al. (run in
  two-group mode; `aldex.glm` hits an R-4.5 S4 bug); STM = Roberts et al. 2014.

## 6. Minor

- Reconcile the Fig. 4 sign-test p (text ≈0.01 vs figure 0.016).
- Fix truncated figure titles (Figs 3–4).
- Align repo URL / package name: the paper's Software section already points to
  `github.com/marcoortu/spectral-gsca-tm` (the repo moved from `irl-egsca-tm`); ensure README and
  Data/Software statements match.
- Two `R/` files (`egscatm_fit.R`, `sgscatm_fit.R`) both define `sgscatm()`; the alphabetically
  later `sgscatm_fit.R` wins. Remove the dead `egscatm_fit.R` definition to avoid an
  `R CMD check` multiple-definition note (logged; not changed this session).
