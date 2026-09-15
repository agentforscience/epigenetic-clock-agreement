# Cloned Repositories

Four epigenetic-clock implementations were cloned so that clock outputs can be
cross-checked between independent codebases rather than trusted from a single
implementation. No repositories were specified in the research topic, so all four
were selected during Phase 1.

Helper scripts written for this project live in `../scripts/`, **not** here.

---

## 1. biolearn — **primary engine**

- **URL**: https://github.com/bio-learn/biolearn
- **Location**: `code/biolearn/` (also installed into the venv as `biolearn==0.9.1`)
- **Language**: Python 3.10+
- **Paper**: https://www.biorxiv.org/content/10.1101/2023.12.02.569722v2
- **License**: see `code/biolearn/LICENSE`

### What it provides
- **68 model definitions** (`biolearn/data/library.yaml`), of which **51 are
  coefficient-based methylation models**. Full index written to
  `code/biolearn_model_index.csv`.
- Reference implementations of Horvath v1/v2, Hannum, PhenoAge, GrimAge V1/V2,
  DunedinPACE / PoAm38, Lin, VidalBralo, Weidner, Garagnani, Bocklandt,
  HRSInCHPhenoAge, CausAge/DamAge/AdaptAge, StocZ/StocP/StocH, AltumAge,
  PCHorvath1, GP-age, EpiTOC1/2, MiAge, DNAmTL, Zhang_10, plus cell-type
  deconvolution (Reinius 450K, Salas EPIC, 12-cell EPIC) and a sex estimator.
- **A GEO data loader with curated metadata parsers** — 48 hand-curated datasets in
  `biolearn/data/library.yaml` plus ~2,200 auto-scanned GEO series. This is the
  single most valuable feature for this project: it downloads a GEO series matrix
  and returns `(dnam, metadata)` with age/sex already parsed.
- A local pickle cache (`biolearn.cache.LocalFolderCache`).

### Verified working
End-to-end test on GSE41169 (95 samples, 450K):

```
GSE41169: dnam shape=(485577, 95)  metadata=['sex','age','disease']   32 s
26/29 panel clocks ran first time; the remaining 3 fixed (see below).
corr(clock, chronological age): Horvathv2 0.958, AltumAge 0.952, Horvathv1 0.935,
Hannum 0.929, PhenoAge 0.912, ... DunedinPoAm38 0.303, DNAmTL -0.841
```

And on the 36 GSE55763 technical-replicate pairs (`scripts/replicate_reliability.py`),
**26/26 panel clocks ran**, reproducing Higgins-Chen et al. 2022's published
reliability figures for the same 36 pairs: PhenoAge median deviation 2.43 y / max
8.58 y (published 2.4 / 8.6), GrimAge ICC 0.9889 (published 0.989), Horvath1
2.08 y / 5.45 y (published 1.8 / 4.8). Output:
`results/GSE55763_clock_reliability.csv`.

Full smoke-test output: `datasets/smoke_GSE41169_clocks.csv`.

### Issues found and how they are handled

| Issue | Detail | Fix |
|---|---|---|
| `DunedinPACE` crashes | `ValueError: assignment destination is read-only` — `dunedin_pace_normalization` passes a read-only `DataFrame.values` view into `quantile_normalize_using_target`, which writes in place. Version-dependent (numpy 2.5 / pandas 3.0). | `scripts/biolearn_patches.py` monkey-patches the function to copy first. **Verified fixed**: DunedinPACE now returns values ~0.93–1.20. |
| `GPAge*` needs an extra dep | `ImportError: GPy is required` | `uv pip install GPy`. Installed, but GP prediction is **very slow** (>25 min CPU on 95 samples) — treat GP-age as optional. |
| `GrimAgeV1/V2` output shape | Returns a multi-column frame (7–9 DNAm components, `Age`, `Female`, `DNAmGrimAge`, `AgeAccelGrim`), not a single `Predicted` column. Naively taking column 0 gives `DNAmADM`, not GrimAge. | Select `DNAmGrimAge` / `AgeAccelGrim` explicitly. **Verified.** |
| `GrimageModel.predict` mutates its input | Appends `Age`, `Female`, `Intercept` rows to `geo_data.dnam` in place, corrupting later models. | Run GrimAge on a copy, or last. |
| **`DNAmTL` ships a sign-flipped intercept** | The PyPI 0.9.1 package has `intercept,-7.924780053` in `biolearn/data/DNAmTL.csv`; the upstream repo (`@0d714f5`) has `+7.924780053`. Every predicted telomere length comes out **negative** (mean −8.35 kb instead of +7.50 kb) — a constant 15.85 kb offset. Caught by the cross-engine check, invisible to any single-engine analysis. Only DNAmTL is affected: the other 76 packaged coefficient files are byte-identical to repo HEAD. | Patched in `scripts/biolearn_patches.py` (restores the positive intercept at load). **Verified fixed.** |
| GP-age CpG lists | `GP-age_model_*.json.zip` files are gzip, not zip, and contain no CpG names — the CpG list is in the model *definition* (`model.sites`). | Handled in `scripts/clock_cpgs.py`. |

### Entry points
```python
from biolearn.data_library import DataLibrary
from biolearn.cache import LocalFolderCache
from biolearn.model_gallery import ModelGallery

data = DataLibrary(cache=LocalFolderCache("datasets/cache", 200)).get("GSE40279").load()
pred = ModelGallery().get("Horvathv1").predict(data)   # DataFrame indexed by sample
```

---

## 2. pyaging — **second Python engine, and the source of clock design metadata**

- **URL**: https://github.com/rsinghlab/pyaging
- **Location**: `code/pyaging/` (also installed as `pyaging==0.5.2`)
- **Paper**: de Lima Camillo et al., *Bioinformatics* 2024, doi:10.1093/bioinformatics/btae200

### What it provides
- **140 human DNA-methylation clocks** (177 total across species and data types) —
  roughly 5× biolearn's methylation panel. Includes clocks biolearn lacks:
  **pchorvath2013, pchannum, pcphenoage, pcgrimage, pcskinandblood** (the full
  Higgins-Chen PC family), **systemsage** and its 11 organ-system sub-clocks,
  **zhangen / zhangblup / zhangmortality**, **stemtoc**, **epitoc3**, **cpgptgrimage3**,
  **dnamfitage**, **retroelementage**, **intrinclock**, **abec/cabec/eabec**, **encen40/100**.
- **`clocks/metadata/clock_metadata.json`** — structured, per-clock design
  descriptors: `training_target`, `model_type`, `tissue`, `platform`, `unit`,
  `n_features`, `population`, `year`, `citation`, `doi`. **This is the covariate
  table direction D3 needs** and it would otherwise have to be hand-curated.
- GPU-optional PyTorch backend; AnnData-based I/O.

### Entry points
```python
import anndata, pyaging as pya
ad = anndata.AnnData(X=betas_samples_by_cpg)          # note: samples x CpGs
pya.pred.predict_age(ad, ["horvath2013", "pcgrimage"], dir="pyaging_data")
ad.obs   # one column per clock
```
Note `predict_age` downloads per-clock weight files into `dir` on first use, so the
working directory matters. Smoke-test log: `logs/pyaging_smoke.log`.

### Cross-engine validation result

`scripts/cross_engine_check.py GSE41169` compares the 22 clocks both engines
implement, on identical imputed input:

Before the DNAmTL patch:
```
22/22 clocks agree at Pearson r = 1.00000
20/22 agree to < 0.0001 in absolute value
Horvathv1 vs horvath2013:  mean |diff| =  0.0102 y    (float32 rounding in pyaging)
DNAmTL    vs dnamtl     :  mean |diff| = 15.8496 kb   <- biolearn coefficient bug
                           biolearn mean -8.345 kb vs pyaging +7.504 kb
```

After the patch:
```
22/22 clocks agree at Pearson r = 1.00000
DNAmTL    vs dnamtl     :  mean |diff| = 0.0000 kb, both engines 7.504 kb
Horvathv1 vs horvath2013:  mean |diff| = 0.0102 y  (float32 rounding, benign)
21/22 agree to < 0.0001 in absolute value
```

This is the single most valuable check in the whole pipeline: it found a real
coefficient-file defect that would have silently corrupted every DNAmTL result, and
it establishes that the remaining disagreement between *different* clocks is not an
implementation artefact.

### Relation to biolearn
Where both implement the same clock (Horvath 2013, Hannum, PhenoAge, DunedinPACE,
skin&blood, Lin, VidalBralo, Ying*, Stoc*, EpiTOC1/2, DNAmTL, deconvolution),
running both gives an **implementation-level control**: any cross-clock disagreement
must survive the check that the *same* clock agrees with itself across engines.

---

## 3. methylclock — R reference implementation

- **URL**: https://github.com/isglobal-brge/methylclock
- **Location**: `code/methylclock/`
- **Paper**: Pelegí-Sisó et al., *Bioinformatics* 2021
- **Language**: R / Bioconductor (not installed — no R toolchain in this venv)

### What it provides
- **44 published clocks** across eight families, organised as a declarative clock
  registry (v2.0 rewrite).
- `ageAcceleration()` for residual and cell-adjusted acceleration, plus the
  canonical **`EEAA()` and `IEAA()` of Chen et al. 2016** and Klemera-Doubal-weighted
  `ageAccelerationChen()`.
- `cellCounts()` — Houseman reference-based deconvolution.
- Out-of-core computation via BigDataStatMeth; four imputation strategies
  (per-CpG mean, clock training reference values, none, KNN).

### Why it is here
It is the canonical *definition* of IEAA/EEAA and of the age-acceleration
conventions this project must follow. The coefficient tables under
`code/methylclock/data/` are readable from Python without R and serve as an
independent check on biolearn's coefficient files.

---

## 4. dnaMethyAge — R reference implementation

- **URL**: https://github.com/yiluyucheng/dnaMethyAge
- **Location**: `code/dnaMethyAge/`
- **Paper**: Wang et al., *GeroScience* 2023, doi:10.1007/s11357-023-00871-w
- **Language**: R (not installed)

### What it provides
17 clocks with an explicit, well-documented table of **trained phenotype**,
**number of CpGs** and **tissue derived** for each — including PCGrimAge (78,464
CpGs), cAge (Bernabeu 2023), the Shireby cortical clock, the pan-mammalian clocks
and the 2025 Intrinsic Capacity clock. Its README table is the cleanest available
summary of clock training targets and is used to cross-check the design-feature
coding for direction D3.

---

## Summary

| Repo | Language | Clocks | Installed | Role |
|------|----------|--------|-----------|------|
| biolearn | Python | 68 models / 51 methylation | ✅ `biolearn==0.9.1` | Primary engine + GEO loader + cell deconvolution |
| pyaging | Python | 177 / 140 human DNAm | ✅ `pyaging==0.5.2` | Expanded panel (PC clocks, SystemsAge, Zhang) + design metadata |
| methylclock | R | 44 | ❌ (source only) | IEAA/EEAA definitions; coefficient cross-check |
| dnaMethyAge | R | 17 | ❌ (source only) | Training-target/tissue reference table |

**Recommendation for the experiment runner.** Use **biolearn** for data loading
(its GEO metadata parsers are the reason the cohorts are usable at all) and run the
clock panel through **both biolearn and pyaging**, taking pyaging's
`clock_metadata.json` as the design-feature table for direction D3. Import
`scripts/biolearn_patches` before any biolearn model call.
