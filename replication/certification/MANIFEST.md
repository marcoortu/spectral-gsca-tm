# Certification manifest — manuscript object → script → gate → certified vs draft

Status legend: **certified** (chain run this session, reduced reps) · **pending** (scaffolded,
not run at full scale this session) · **stands** (prior pilot-based value, unchanged).

| Manuscript object | Generating script | Gate | Draft number | Chain-certified number | Status |
|---|---|---|---|---|---|
| Estimator `sgscatm_chain()` | `R/chain.R` (+ anchors/orient/refine/sandwich) | — | pilot only | 3-stage chain in package, 17/17 tests | **certified** |
| Raw-pilot collapse (Cor. 11) | `cert_sim.R` G0 | G0 | ∝M^{-1/2} | 0.060→0.042→0.025 (M=1e3,2e3,5e3) | **certified** ✓ |
| Table 1 principal angle | `cert_sim.R` G0 | G0 | <1e-5 | 0.89–1.25 rad (leading dirs only) | **certified** ✗ |
| Table 1 refined→oracle RMSE | `cert_sim.R` G0 | G0 | refined→oracle floor | chain 0.067 vs oracle 0.015 (M=5e3, ≈4.6×) | **certified** (gap>2×) |
| Start-independence | `cert_sim.R` G7 | G7 | — | 4.18e-4 ≤ 5e-4 | **certified** ✓ |
| In-regime coverage (NEW) | `cert_sim.R` G2c | G2c | (new) nominal | **cov 0.346, SE/SD 0.19** | **certified** ✗ (key) |
| Coverage mechanism | `cert_G2c_mechanism.R` | G2c | — | feasible vs oracle orientation | **certified** (see log) |
| Table 2 coverage vs √M/L | (scaffold) | G2a/b | .95/.94/.78/.61 | — | **pending** |
| Table 3 crossover vs STM | (scaffold; feasibility F4) | G3 | feas 0.058/0.126→0.0015/0.0038 | — | **pending** |
| Fig 1 SE/SD vs b_max | (scaffold) | G2a | ~1±10% | — | **pending** |
| b(z) / G0 gradient | (round-5 harness) | G8 | matches closed form | — | **pending** |
| Table 4 delocalization (CRC) | `06_phase2_fit.R` | G1 | r 6.55/7.20 | pilot value stands | **stands** |
| Table 5 Wald vs PERMANOVA (CRC) | `06/07_phase2` | G5a | 0.0065/0.001 … | needs chain re-run | **pending** (chain) |
| Table 6 speed (CRC) | `07b/08_phase2` | G6 | 38–173× | needs chain re-run | **pending** (chain) |
| Fig 2 forest (CRC) | `08_phase2_figures.R` | — | pilot+sgscatm_vcov | needs chain+sandwich | **pending** (chain) |
| Fig 3 SE vs bootstrap (CRC) | `06/08_phase2` | G4 | median 1.15 (pilot) | needs chain sandwich | **pending** (chain) |
| Fig 4 known genera (CRC) | `07b/08_phase2` | G5b | 6/6, p≈0.016 | pilot value stands | **stands** |
| Table 7 / Fig 5 (BES) | `bes_strong_weak.R` | — | 0.962/0.877, 0.961/0.781 | out-of-regime; pilot values stand | **stands** |

## Notes

- The **only manuscript numbers this session moves** are the certification (sim) gates: the raw
  collapse and start-independence reproduce; the Table-1 principal angle and the new in-regime
  coverage (G2c) do **not** reproduce as claimed and are the reportable findings.
- The CRC/BES application numbers (Tables 4–7, Figs 2–5) were produced by the **raw pilot +
  `sgscatm_vcov`**; re-running them on `sgscatm_chain()` + the (anchor-variance-aware) sandwich is
  the outstanding work, gated on fixing the coverage SE first (see `PAPER_EDITS.md` §3).
- Full deterministic reproduction (seeds `base + regime·1000 + rep`) is scripted in
  `replication/feasibility/04_run.R` (F1–F4) and `replication/certification/cert_sim.R`.

---

## Phase 1 re-certification addendum (frozen-Φ chain, `PREREG_PHASE1.md`)

| object | gate | draft (pilot) | chain-certified | status |
|--------|------|---------------|-----------------|--------|
| feasible in-regime coverage (T2-inreg) | key | (new) | **cov 0.49–0.60, SE/SD 0.13–0.16** | **FAIL** — sandwich under-covers feasible even clean-anchor/large-L; bootstrap required |
| CRC delocalization (T4/G1) | G1 | r 6.5 / 7.2 | **r ≈ 10,400 / 10,700** | MOVED — chain natural-scale scores delocalize; refinement hits sweep cap |
| CRC known genera (G5b) | G5b | 6/6, p=0.016 | **5/6, p=0.11** | MOVED — directional, not significant |
| CRC speed (T6/G6) | G6 | 0.02 s, 38–173× | **0.99 s ≈ competitors** | MOVED — 38–173× was the pilot, not the chain |
| CRC study_condition Wald (G5a/G4) | G5a | — | p 0.0097 (sandwich)/0.006 (boot), PERMANOVA 0.001 | STANDS — disease signal robust; BMI SE-sensitive |
| package | — | pilot only | B-stationarity backstop; dead `egscatm_fit.R` removed | fixed, documented |
| T1/T2/T3/BES/G8 | — | — | — | NOT RUN (decisive gate failed; budget) |

**Sections whose numbers move if the chain is adopted** (for the deferred PAPER_EDITS step): §5 Table 1
(principal cosines, not <1e-5 angle), §5.3/5.6 (feasible coverage is NOT calibrated by the sandwich —
bootstrap required; not the draft's nominal claim), §6.1 (CRC delocalization r, known-genera p, speed —
all move; the 38–173× speedup is the pilot, not the chain), §6.1.0.3 (SE attribution). **Recommendation
on the estimator/SE:** the raw pilot + `sgscatm_vcov` (original Section 6) is the stronger *calibrated*
deliverable; the frozen-Φ chain is theory-faithful but its feasible closed-form inference is not
calibrated (T2-inreg FAIL) and it degrades on the CRC data — report as the honest state of Cor. 18's
second gate, not as a solved closed-form-inference result.
