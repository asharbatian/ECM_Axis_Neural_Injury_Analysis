# Independent GEO validation of the ECM-axis hypothesis

Goal: test, with **independent public datasets**, whether the ECM axes defined in this
project (esp. the Hyaluronan / ECM-DAMP axis) are engaged after CNS injury, and build a
reproducible, auditable evidence set. Discovery sets (GSE35338, GSE5296) are *not* reused here.

## Cohorts (all verified against GEO metadata)
| GEO | Species | Model | Contrast used | n (inj vs ctrl) |
|---|---|---|---|---|
| GSE45006 | rat | SCI, clip impact-compression, GPL1355 | D1/D3/D7/D14/D56 vs sham | 4 vs 4 each |
| GSE47681 | mouse | SCI contusion, GPL1261 (WT only) | D1/D3/D7 vs sham | 3-4 vs 4 |
| GSE45997 | rat | unilateral TBI, GPL1355 | ipsilateral vs naive (and vs contralateral) | 3 vs 3 |

**QC decision:** two GSE47681 samples have conflicting title vs. characteristics
(GSM1154530 "DAY1WT rep1" = Day 7; GSM1154541 "DAY1WT rep2" = control) and were excluded.

## Method (pre-specified before looking at results)
1. Axis gene sets: unchanged from `scripts/02_annotation_mapping.R` (HA axis split into synthesis / degradation / receptor-DAMP sub-modules).
2. Probes collapsed to genes (highest mean probe); log2 expression.
3. Gene-level variance-moderated effect size (SAM-style s0 = median SD) injured vs control.
4. Axis statistic = mean effect size; **competitive permutation test** (10,000 random gene sets, same size, matched on mean-expression decile) -> z and one-sided p.
5. **Primary pre-specified contrasts:** GSE45006 D3, GSE47681 D3, GSE45997 ipsi-vs-naive; Stouffer meta-analysis across the 3 cohorts. All other contrasts are secondary/descriptive.

Reproduce: `python3 geo_validation/run_geo_validation.py && python3 geo_validation/make_figures.py` (needs network to `ftp.ncbi.nlm.nih.gov`; deps: pandas, scipy, numpy, matplotlib, requests).

## Results
Primary meta-analysis (k = 3 cohorts, `results/meta_primary_contrasts.csv`):

| Axis | Stouffer z | BH q |
|---|---|---|
| Provisional matrix | 9.6 | <1e-4 |
| HA receptor/DAMP arm (Cd44, Tlr2, Tlr4, Cd14, Hmmr) | 8.3 | <1e-4 |
| **Hyaluronan (full axis)** | **7.2** | **<1e-4** |
| Proteases/regulators | 6.1 | <1e-4 |
| Crosslinking/fibrosis | 6.0 | <1e-4 |
| Basement membrane | 2.5 | 0.009 |
| HA degradation | 1.7 | 0.06 |
| PNN/CSPG | -1.2 | 0.89 |

HA-axis z is positive in **every one of 11 contrasts**, including late SCI (D56, GSE45006).
Figures: `results/fig1_axis_heatmap.png`, `results/fig2_HA_axis_by_contrast.png`.

## What this does and does NOT support (read before citing)
Supported
* The HA/ECM-DAMP axis is activated after SCI (two species, two labs/platforms) and, more weakly, TBI; persists up to 8 weeks in SCI.
* The signal is carried by the **DAMP-receptor arm** (Cd44, Tlr2, Tlr4, Cd14) - consistent with an "ECM-DAMP signalling" framing.

Not supported / caveats
* **HA synthesis (Has1/2/3) is not induced** in these data (z ~ 0 in GSE47681; Has2 q>0.1). Claims about increased HA production are not backed by these cohorts. HA synthesis gene coverage on GPL1355 is only 1/3 genes.
* **Not specific to HA**: Provisional matrix ranks higher than HA in the meta-analysis. Injury broadly activates ECM remodelling.
* Tlr2/Tlr4/Cd14/Cd44 are also microglia/macrophage genes; the HA-receptor signal may partly reflect immune-cell infiltration. Needs a control (immune-marker adjustment) or single-cell/cell-sorted data (e.g. GSE96055 microglia SCI time course).
* GSE45997 (TBI) is n = 3 vs 3 with unknown timepoint; effect is modest.
* Permutation null samples genes independently, whereas axis genes are co-regulated, so z-scores are likely anti-conservative; treat p-values as relative evidence, effect sizes/consistency across cohorts are more informative. n = 3-4 per group.
* Only 2 SCI + 1 TBI cohorts; no stroke or implant cohort yet.

## Candidate next datasets (verified to exist, not yet analysed)
GSE96055 (SCI microglia RNA-seq time course), GSE80174 (rat LFPI RNA-seq, 3 mo; paper highlights CD44), GSE2871 (rat TBI 4h/24h, RG-U34A - poor ECM coverage), GSE464 (rat SCI, 544 arrays), GSE32529 (mouse ischemia, 224 arrays), GSE298545 / GSE319931 (recent SCI RNA-seq).
