# Research Planning — Do Epigenetic Age Acceleration Estimates Agree Across Clock Algorithms?

Phase: `resource_finder` (Phase 1). This document records the direction search,
the scoring, the retained directions, and the pruned ones.

## 1. Hypothesis under test

> Different epigenetic clock algorithms produce age-acceleration estimates that
> correlate at **r < 0.5** across individuals in the same cohort, and the ranking
> of which individuals are biologically older changes in **more than 30%** of cases
> when switching clocks — indicating that epigenetic age acceleration is partly a
> property of the clock algorithm rather than a stable biological measurement.

Two distinct claims are bundled here and must be separated experimentally:

| # | Claim | Operationalisation |
|---|-------|--------------------|
| H1 | Cross-clock agreement of age acceleration is low | Pairwise Pearson/Spearman r of AA residuals; fraction of pairs with r < 0.5 |
| H2 | Individual-level rankings are unstable | Rank-discordance rate; quantile reclassification rate; Kendall τ |
| H3 (implied) | The disagreement is *algorithmic*, not just measurement noise | Disattenuate cross-clock r by each clock's technical ICC; compare to the noise ceiling |

H3 is the load-bearing claim of the title ("property of the clock algorithm rather
than a stable biological measurement"). Without it, low r is equally consistent
with "every clock is a noisy measurement of one underlying quantity".

## 2. Literature anchors that constrain the design

Evidence gathered in Phase 1 (see `literature_review.md` for full detail):

- **Crimmins et al. 2021 (HRS, N=4,018, 13 clocks, EPIC).** Correlations *among
  age-acceleration residuals* are markedly lower than among raw epigenetic ages.
  Highest AA pair: DunedinPoAm38 × GrimAgeAccel **r = 0.64**. PhenoAgeAccel ×
  GrimAgeAccel **r = 0.35**; PhenoAgeAccel × DunedinPoAm38 **r = 0.25**. First-generation
  clocks (Horvath1/Horvath2/Hannum) correlate relatively well with each other but
  weakly with GrimAge/DunedinPoAm.
- **Li et al. 2020 (SATSA, N=845, nine biological ages).** After regressing out
  chronological age, only moderate residual correlations remain — Horvath × Hannum
  DNAmAge **r = 0.35**.
- **Belsky et al. 2018 (Dunedin, N=964, 11 measures).** "Contrary to expectation, we
  found low agreement between different measures of biological aging… various
  proposed approaches… may not measure the same aspects of the aging process."
- **Jansen et al. 2021 (NESDA, ~3,000, five omics clocks).** All cross-clock
  r < 0.2 (across omics layers, so an upper bound on how low within-omics can go).
- **Higgins-Chen et al. 2022 (Nature Aging).** Technical noise alone produces
  replicate-to-replicate deviations of **up to 9 years** (median 0.9–2.4 y,
  max 4.5–8.6 y) for six prominent clocks, against an AA standard deviation of only
  3–5 years. PC-based retraining raises AA ICC to 0.97–0.99. **This sets the noise
  ceiling H3 must be tested against.**
- **Sugden et al. 2020 (Patterns).** Probe-level reliability on 450K/EPIC is highly
  heterogeneous; unreliable probes attenuate associations and generate false negatives.
- **Tong et al. 2024 (Nature Aging).** ~66–75% of Horvath-clock accuracy (90% for
  Zhang) is reproducible by a purely stochastic DNAm-change model; only 63% for
  PhenoAge — i.e. clocks differ in how much non-stochastic biology they encode.

**Implication for scope.** The qualitative answer ("clocks disagree") is already
established. The contribution available to this project is *quantitative and
mechanistic*: measure the agreement structure over a much wider clock panel
(29 clocks vs 13 in the largest prior study), in multiple independent public
cohorts, and decompose it into (a) technical-noise floor, (b) algorithm-design
factors. That framing is what the direction ranking below optimises for.

## 3. Direction search and scoring

Scored 1–5 on four axes; **Total** is the unweighted sum (max 20).

- **Ev** — evidence base / literature support that the direction is answerable
- **Rel** — relevance to the stated hypothesis
- **IG** — expected information gain beyond what the literature already reports
- **Feas** — implementation feasibility inside this workspace (data, compute, tooling)

| ID | Direction | Ev | Rel | IG | Feas | Total | Verdict |
|----|-----------|----|-----|----|------|-------|---------|
| **D1** | **Cross-clock agreement + rank stability of AA**, 29 clocks × ≥3 independent public cohorts (450K and EPIC) | 5 | 5 | 4 | 5 | **19** | **KEEP (primary)** |
| **D2** | **Noise-floor decomposition**: per-clock technical ICC from GSE55763 replicate pairs; disattenuated cross-clock r; compare observed disagreement with the reliability ceiling | 5 | 5 | 5 | 4 | **19** | **KEEP** |
| **D3** | **Structural explanation of (dis)agreement**: predict pairwise AA r from CpG-set Jaccard overlap, training target (age vs mortality vs pace), training tissue, clock generation, n-CpG | 4 | 5 | 5 | 5 | **19** | **KEEP** |
| D4 | Cross-cohort replication of the agreement matrix | 4 | 4 | 3 | 5 | 16 | **Folded into D1** (a robustness axis, not a separate question) |
| D5 | Cell-composition confounding (IEAA/EEAA-style adjustment via Houseman/Salas deconvolution) | 5 | 3 | 3 | 4 | 15 | Pruned — a covariate-adjustment sensitivity analysis, not an independent direction; biolearn deconvolution is available if D1 needs it |
| D6 | Downstream consequence: do different clocks yield different sex/smoking/BMI/disease associations? | 5 | 3 | 3 | 4 | 15 | Pruned — Crimmins 2021 already reports exactly this for 13 clocks; low marginal IG |
| D7 | Preprocessing/normalisation sensitivity (noob, BMIQ, funnorm) | 4 | 2 | 3 | 2 | 11 | Pruned — needs IDATs, not beta matrices; McEwen 2018 covers it |
| D8 | Build a consensus/latent "meta-clock" (PCA over AA across clocks) | 3 | 2 | 3 | 5 | 13 | Pruned as a direction; retained as a *descriptive* output of D1 (variance explained by PC1 quantifies how much shared signal exists) |
| D9 | Missing-CpG / imputation sensitivity | 3 | 2 | 3 | 4 | 12 | Pruned — becomes a documented confound control in D1 (fixed imputation policy) |
| D10 | Cross-tissue (buccal/saliva vs blood) agreement | 4 | 2 | 3 | 2 | 11 | Pruned — different question; covered by the 2025 cross-tissue comparison paper |
| D11 | Longitudinal within-person stability of AA | 4 | 3 | 4 | 1 | 12 | Pruned — no suitable public longitudinal 450K/EPIC series identified with repeat draws |
| D12 | Non-human / mammalian clocks | 3 | 1 | 2 | 2 | 8 | Pruned — out of scope |

### Retained directions (the direction budget of 3)

**D1 — Agreement and rank stability (primary test of H1 + H2).**
For each cohort: run all 29 clocks, define age acceleration as the residual of
clock output on chronological age (and, for rate clocks, the raw value plus its
age residual), then compute the full pairwise Pearson and Spearman matrix, the
fraction of pairs below r = 0.5, Kendall τ, the fraction of individual pairs whose
relative order flips between clocks, and top/bottom-quartile reclassification
rates. Replicate the whole matrix in ≥3 cohorts spanning 450K and EPIC.

**D2 — Noise floor (test of H3).**
GSE55763 contains 36 whole-blood samples assayed in duplicate across separate
batches (age range 37.3–74.6), the same benchmark used by Higgins-Chen et al.
Compute each clock's AA ICC on those 36 pairs, then disattenuate the D1
correlations: `r_true ≈ r_obs / sqrt(ICC_a · ICC_b)`. If disattenuated r stays
well below 1, the disagreement is *not* explained by measurement noise and the
hypothesis's "property of the algorithm" reading is supported; if it approaches 1,
the honest conclusion is that clocks are noisy measurements of a shared quantity.

**D3 — What predicts agreement.**
Using the extracted coefficient files, compute for every clock pair: Jaccard
overlap of CpG sets, cosine similarity of signed weights on shared CpGs, same/different
training outcome, same/different training tissue, generation (1st = trained on
chronological age, 2nd = trained on mortality/phenotype, 3rd = trained on pace).
Regress observed AA correlation on these design features. This turns "clocks
disagree" into "clocks disagree *because* of X".

### Re-opening the search space
Per the direction budget, the search space stays fixed unless new evidence
invalidates the ranking. Any change must be recorded in `STATE.md` with the
evidence that forced it.

## 4. Planned data (see `datasets/README.md` for download instructions)

| Role | Dataset | N | Platform | Notes |
|------|---------|---|----------|-------|
| Primary cohort | GSE40279 | 656 | 450K | Whole blood, ages 19–101, the Hannum training cohort; widest age range |
| Replication A | GSE42861 | 689 | 450K | Whole blood, RA cases + controls, age + sex + smoking |
| Replication B | GSE157131 | 946 | 450K | Peripheral blood leukocytes, age + sex + BP |
| Replication C | GSE51032 | 845 | 450K | EPIC-Italy cohort |
| EPIC platform check | GSE132203 | 795 | EPIC | Grady Trauma Project; tests platform dependence |
| Pilot / fast iteration | GSE41169 | 95 | 450K | Already cached and validated end-to-end |
| **Technical replicates (D2)** | GSE55763 | 72 (36 pairs) | 450K | Duplicate assays in separate batches; the Higgins-Chen reliability benchmark |

## 5. Known methodological traps to control

1. **AA definition.** Residual-from-regression vs simple difference (DNAmAge − age)
   give different correlations. Use residuals as primary (the field standard),
   report differences as a sensitivity analysis.
2. **Age-range dependence.** Correlations among AA measures depend on the cohort's
   age spread. Report age range per cohort and stratify.
3. **Shared chronological-age variance.** Raw epigenetic ages correlate highly
   simply because they all track age; only the residuals are informative. This is
   exactly the gap between Crimmins Fig 1A and Fig 1B.
4. **Clock-dependent imputation.** biolearn applies per-model imputation of missing
   CpGs (`sesame_450k` / `dunedin` / `averaging` gold standards). Fix one policy
   across all clocks and record it; a differing policy would itself create
   disagreement.
5. **GrimAge uses age and sex as inputs**, so its "acceleration" is not comparable
   in kind to clocks that use methylation only. Flag it separately.
6. **Non-age-scaled outputs** (DunedinPACE, EpiTOC, MiAge, DNAmTL, Zhang_10) are
   rates/risks, not ages. They belong in the matrix but must be z-scored and
   labelled, not treated as "years".
7. **biolearn mutates `geo_data.dnam` inside `GrimageModel.predict`** (appends Age,
   Female and Intercept rows). Run GrimAge on a copy or last in the sequence.

## 6. Success criteria for Phase 2

- Full clock × sample matrices for ≥3 cohorts persisted to `results/`
- Pairwise agreement matrices (Pearson, Spearman, Kendall) with bootstrap CIs
- Rank-flip and quartile-reclassification rates with CIs — the direct H2 numbers
- Per-clock technical ICC from GSE55763 and disattenuated correlations — the H3 number
- Design-feature regression explaining pairwise agreement — the D3 number
- An explicit verdict on each of H1/H2/H3, including the case where the data
  contradict the stated hypothesis

## 7. What Phase 1 already established

Three validation checks were run during resource gathering, so Phase 2 starts from a
verified pipeline rather than an untested one:

1. **Cross-engine agreement.** All 22 clocks implemented by both biolearn and pyaging
   agree at r = 1.00000 on identical input. This matters because it rules out the
   trivial alternative explanation for D1 — that apparent cross-clock disagreement is
   an implementation artefact. It also caught a genuine coefficient-file defect
   (biolearn's `DNAmTL` intercept sign; patched and re-verified).

2. **External benchmark.** Running the panel on the 36 GSE55763 replicate pairs
   reproduces Higgins-Chen et al. 2022's published reliability figures for the same
   pairs (PhenoAge 2.43 y median / 8.58 y max vs published 2.4 / 8.6; GrimAge ICC
   0.9889 vs 0.989).

3. **D2's inputs are computed.** `results/GSE55763_clock_reliability.csv` holds
   per-clock technical ICC for epigenetic age and for age acceleration. Headline:
   **AA ICC spans 0.714–0.994 and 15 of 26 clocks fall below 0.90.** PhenoAge — a
   first-tier clock — has AA ICC 0.758. Disattenuation is therefore a material
   correction, not a formality: an observed r = 0.35 between two clocks at ICC ≈ 0.75
   corrects to ≈ 0.47.

This changes one thing about the plan: **D2 should be reported alongside D1 rather
than after it**, because the raw correlation matrix is uninterpretable without the
per-clock reliability ceiling sitting next to it.
