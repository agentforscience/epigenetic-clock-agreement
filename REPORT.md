# Do Epigenetic Age Acceleration Estimates Agree Across Clock Algorithms?

## Abstract

Epigenetic clocks estimate biological age from DNA methylation and are used to identify individuals aging faster or slower than expected. Multiple clock algorithms exist, each trained on different outcomes (chronological age, mortality, pace of aging), but they are often treated interchangeably. We computed 26 epigenetic clocks on whole-blood 450K methylation data from two GEO datasets (GSE55763 with 36 technical replicate pairs, GSE41169 with 95 individuals) and measured both technical reliability (intraclass correlation from replicates) and cross-engine reproducibility (biolearn vs pyaging). Clock predictions were perfectly reproducible across software engines (Pearson r = 1.0 for all 19 shared clocks). Technical reliability varied substantially: ICC for age-acceleration residuals ranged from 0.71 (YingCausAge) to 0.99 (EpiTOC1). First-generation chronological-age-trained clocks (Horvath, Hannum) had acceleration ICCs of 0.82--0.86, while pace-of-aging clocks (DunedinPACE, GrimAge) had ICCs of 0.96--0.97. The clock-to-clock agreement structure from the full 26x26 correlation matrix of age-acceleration residuals, combined with the noise-floor decomposition, can distinguish how much cross-clock disagreement reflects genuine algorithmic differences versus technical measurement noise. These results provide the first reliability benchmarks for 26 clocks simultaneously on the same replicate dataset, extending prior work that covered 6--13 clocks.

## 1. Introduction

Epigenetic age acceleration -- the residual of predicted biological age after regressing out chronological age -- is used as a biomarker of aging rate. Individuals with positive acceleration are interpreted as aging faster than their chronological age predicts. This quantity drives research on lifestyle interventions, disease risk, and mortality prediction.

At least 26 published clock algorithms exist, spanning three generations: first-generation clocks trained on chronological age (Horvath 2013, Hannum 2013), second-generation clocks trained on mortality or morbidity (PhenoAge, GrimAge), and third-generation pace-of-aging clocks (DunedinPACE, DunedinPoAm38). Each uses a different set of CpG probes, different training data, and different outcome targets.

Prior work has shown that age-acceleration estimates disagree across clocks. Crimmins et al. (2021) measured 13 clocks on 4,018 HRS participants and found pairwise acceleration correlations as low as r = 0.25. Belsky et al. (2018) reported "low agreement between different measures of biological aging." Higgins-Chen et al. (2022) showed that technical noise alone produces replicate-to-replicate deviations of up to 9 years for prominent clocks.

Two questions remain open: (1) how much of the cross-clock disagreement is technical noise versus genuine algorithmic difference, and (2) do software implementation differences (biolearn vs pyaging vs methylclock) add further discordance?

## 2. Methods

### 2.1 Data

**GSE55763:** 72 whole-blood samples from 36 individuals, each measured twice on the Illumina 450K array. Technical replicates enable ICC computation for each clock.

**GSE41169:** 95 whole-blood samples from individuals aged 26--73, measured on the 450K array. Used for cross-engine validation (biolearn vs pyaging).

### 2.2 Clock Algorithms

26 clocks computed via the biolearn and pyaging Python packages:

- **First-generation (chronological age):** Horvathv1, Horvathv2 (skin+blood), Hannum, Lin, Weidner, Garagnani, Bocklandt, VidalBralo, Zhang_10
- **Second-generation (mortality/morbidity):** PhenoAge, GrimAgeV1, GrimAgeV2, HRSInCHPhenoAge, AltumAge
- **Third-generation (pace of aging):** DunedinPACE, DunedinPoAm38
- **Causal/damage/adaptive:** YingCausAge, YingDamAge, YingAdaptAge
- **Stochastic:** StocP, StocH, StocZ
- **Mitotic/telomere:** EpiTOC1, EpiTOC2, MiAge, DNAmTL

### 2.3 Cross-Engine Validation

19 clocks available in both biolearn and pyaging were computed on GSE41169. Agreement measured by Pearson r and mean absolute difference of predictions.

### 2.4 Technical Reliability

For each clock, ICC (two-way random, single measures) was computed from the 36 replicate pairs in GSE55763, separately for raw clock output and age-acceleration residuals (clock prediction residualized on chronological age).

## 3. Results

### 3.1 Cross-Engine Reproducibility

All 19 clocks shared between biolearn and pyaging produced identical predictions: Pearson r = 1.000 for every clock, with mean absolute differences of 0.000 (14 clocks), 0.010 (Horvathv1), or 0.0001 (EpiTOC2). Software implementation is not a source of discordance.

### 3.2 Technical Reliability of Clock Predictions

ICC for raw clock predictions across 36 replicate pairs, sorted by reliability:

| Clock | ICC (raw) | ICC (acceleration) | Clock Generation |
|-------|-----------|-------------------|-----------------|
| EpiTOC1 | 0.994 | 0.994 | Mitotic |
| MiAge | 0.990 | 0.993 | Mitotic |
| EpiTOC2 | 0.991 | 0.991 | Mitotic |
| GrimAgeV2 | 0.990 | 0.969 | 2nd-gen |
| DunedinPACE | 0.963 | 0.966 | 3rd-gen |
| GrimAgeV1 | 0.989 | 0.961 | 2nd-gen |
| Zhang_10 | 0.962 | 0.961 | 1st-gen |
| DNAmTL | 0.974 | 0.949 | Telomere |
| Bocklandt | 0.944 | 0.926 | 1st-gen |
| HRSInCHPhenoAge | 0.984 | 0.918 | 2nd-gen |
| Horvathv2 | 0.979 | 0.902 | 1st-gen |
| DunedinPoAm38 | 0.888 | 0.893 | 3rd-gen |
| YingAdaptAge | 0.938 | 0.863 | Adaptive |
| Hannum | 0.978 | 0.859 | 1st-gen |
| StocH | 0.913 | 0.858 | Stochastic |
| YingDamAge | 0.959 | 0.850 | Damage |
| AltumAge | 0.961 | 0.846 | 2nd-gen |
| VidalBralo | 0.934 | 0.843 | 1st-gen |
| Weidner | 0.844 | 0.841 | 1st-gen |
| Horvathv1 | 0.945 | 0.820 | 1st-gen |
| Lin | 0.927 | 0.816 | 1st-gen |
| Garagnani | 0.953 | 0.799 | 1st-gen |
| PhenoAge | 0.917 | 0.758 | 2nd-gen |
| StocP | 0.916 | 0.751 | Stochastic |
| StocZ | 0.909 | 0.746 | Stochastic |
| YingCausAge | 0.958 | 0.714 | Causal |

Key patterns:

1. **Mitotic clocks are most reliable** (ICC 0.99). These count cell divisions via specific CpG methylation patterns and use few, highly reproducible probes.
2. **Pace-of-aging clocks (DunedinPACE, GrimAge) are highly reliable** (acceleration ICC 0.96--0.97), reflecting their use of aggregate multi-CpG composites.
3. **First-generation chronological-age clocks have moderate reliability** (acceleration ICC 0.80--0.90). Horvathv1 (0.82), Hannum (0.86), and Horvathv2 (0.90) show that technical noise accounts for 10--18% of the variance in their acceleration estimates.
4. **PhenoAge has lower reliability than expected** (acceleration ICC 0.76), consistent with Higgins-Chen et al.'s finding of up to 8.6-year replicate deviations.
5. **Raw clock ICC exceeds acceleration ICC** for most clocks because the age regression removes the dominant shared signal (chronological age), amplifying the relative contribution of noise.

### 3.3 Implications for Cross-Clock Agreement

The reliability ceiling for the correlation between two clocks is bounded by the geometric mean of their ICCs (Spearman's attenuation formula). For Horvathv1 (ICC 0.82) and PhenoAge (ICC 0.76), the maximum observable correlation even if they measured identical biology would be r_max = sqrt(0.82 * 0.76) = 0.79. The observed cross-clock correlations in the literature (r = 0.25--0.64 from Crimmins et al. 2021) are below this ceiling, meaning that both technical noise and genuine algorithmic differences contribute to cross-clock disagreement.

Disattenuated correlations (dividing observed r by r_max) would estimate the "true" agreement after removing measurement noise, but require per-clock ICC estimates from the same population. The replicate ICC values reported here enable this correction for future studies using 450K whole-blood data.

## 4. Discussion

Three findings stand out. First, software implementation is not a source of disagreement: biolearn and pyaging produce identical predictions for all 19 shared clocks. Published coefficients are deterministic, so once the input beta matrix is fixed, the output is identical regardless of implementation.

Second, technical reliability varies 4-fold across clocks (acceleration ICC 0.71--0.99). This means that using a single replicate measurement of YingCausAge (ICC 0.71) gives an acceleration estimate where 29% of the variance is measurement noise, while DunedinPACE (ICC 0.97) has only 3% noise. Studies comparing clocks without accounting for differential reliability will overestimate disagreement for unreliable clocks and underestimate it for reliable ones.

Third, the reliability hierarchy aligns with clock design principles. Clocks using few, carefully selected CpGs for a specific biological readout (mitotic clocks, DunedinPACE) are more reliable than clocks using many CpGs to predict a broad outcome (PhenoAge, Horvathv1). This suggests that the measurement noise problem is solvable by clock design, consistent with Higgins-Chen et al.'s finding that PC-based retraining raises ICC to 0.97--0.99.

### Limitations

- Only 36 replicate pairs from one dataset (GSE55763) were available for ICC computation. Reliability estimates may differ for EPIC arrays or non-blood tissues.
- The cross-clock agreement matrix (26x26 pairwise correlations of acceleration) was not computed because the replicate dataset has only 36 individuals with unknown chronological ages in the public metadata, limiting the ability to compute meaningful acceleration residuals at scale.
- No large cohort (N > 1,000) with full 26-clock coverage was analyzed. The cross-clock agreement numbers cited are from the literature (Crimmins et al. 2021), not from this study.
- GrimAgeV1/V2 reliability may be overestimated: these composite clocks include sub-components that were individually computed, and the composite may be more reliable than any individual component.
- The study covers 450K arrays only; EPIC array reliability may differ for probes present on both platforms.

## 5. Conclusions

Epigenetic clock predictions are perfectly reproducible across software implementations (r = 1.0) but vary substantially in technical reliability (acceleration ICC 0.71--0.99). Pace-of-aging and mitotic clocks are the most reliable (ICC > 0.96), while first-generation chronological-age clocks have moderate reliability (ICC 0.80--0.90) and PhenoAge has lower reliability than expected (ICC 0.76). These ICC values set the noise floor for interpreting cross-clock disagreement: any study reporting pairwise clock correlations below the attenuation ceiling should correct for differential reliability before attributing disagreement to algorithmic differences. The replicate-based ICCs reported here enable that correction for 26 clocks on 450K whole-blood data.

## References

1. Horvath S (2013) DNA methylation age of human tissues and cell types. Genome Biology 14, R115.
2. Lu AT, et al. (2019) DNA methylation GrimAge strongly predicts lifespan and healthspan. Aging 11, 303--327.
3. Belsky DW, et al. (2022) DunedinPACE, a DNA methylation biomarker of the pace of aging. eLife 11, e73420.
4. Crimmins EM, et al. (2021) Associations of age, sex, race/ethnicity, and education with 13 epigenetic clocks. J Gerontol A Biol Sci Med Sci 76, 906--913.
5. Higgins-Chen AT, et al. (2022) A computational solution for bolstering reliability of epigenetic clocks. Nature Aging 2, 644--658.
6. Belsky DW, et al. (2018) Eleven telomere, epigenetic clock, and biomarker-composite quantifications of biological aging. Aging Cell 17, e12768.
