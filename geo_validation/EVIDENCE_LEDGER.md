# Evidence ledger: public GEO data vs. the HA-molecular-weight review

Review under test: "A Review of Hyaluronic Acid at the Intersection of Tissue Architecture, Immune Signaling
and Metabolism" (Sharbatian, Ashouri, Stieglitz).

## What GEO expression data can and cannot test
The review's central claim is about HA **molecular weight** (HMW restrains, LMW-fragments drive innate->adaptive
transition). Gene expression does not measure polymer size, abundance or flux. GEO data can therefore test only
*transcript-level* predictions about the synthesis / stabilization / fragmentation / clearance machinery and the
receptor arm (CD44, TLR2/4, CD14, HMMR). The checkpoint hypothesis itself is NOT testable here and nothing below
should be cited as direct evidence for it. The review itself says it is "not established" and calls for paired
size/receptor/outcome measurements.

## Ledger (every line pre-specified; see PREREGISTRATION_*.md and commit timestamps)
| # | Claim from the review | Data | Verdict | Strength |
|---|---|---|---|---|
| E1 | CNS = stabilization module high (HAPLN) | mouse atlas GSE10246 (46 tissues); human GSE1133 | **Supported** (mouse q=3e-4, human q=0.004) | Strong, but partly known biology: Hapln2/4 are brain-specific, so this is a consistency check, not novel support |
| E2 | CNS has little fragmentation (HYAL) | mouse atlas; human | **Supported in mouse** (q=3e-4; amygdala ranks last of 46); **not supported in human** (p=0.51) | Mixed |
| E3 | CNS "unique" stabilization-without-fragmentation | mouse; human | **Supported** (mouse q=3e-4; human q=0.003) | Moderate (driven by E1) |
| E4 | CNS has little clearance (CEMIP/TMEM2) | mouse; human | **Not supported** (p=0.11 / 0.86) | - |
| E5 | Liver/spleen lead fragmentation | mouse | **Liver yes (rank 1/46), spleen no (28/46)**; group test n.s. (p=0.18) | Partial |
| E6 | GI dominates clearance; bladder #1, adipose/lung top synthesis | mouse | **Not supported** (GI p=0.59; bladder ranks 8/46; adipose p=0.67) | Fails |
| E7 | Injury -> HA receptor/DAMP sensing up | SCI rat GSE45006, SCI mouse GSE47681, TBI rat GSE45997 (whole tissue) | **Supported** in all 10 contrasts; meta Stouffer z=7.2 (HA axis), 8.3 (receptor arm) | Moderate (n=3-4/group; see caveats) |
| E8 | HA *synthesis* rises in injury | same | **Not supported** (Has genes ~0) | - |
| E9 | E7 is cell-intrinsic, not just infiltration | sorted microglia, SCI, GSE96055 | **Partly supported**: receptor arm up in microglia (z 2.2-3.0, p 0.014-0.022 at 72h/1w/2w), driven by Cd44 (+2.1 log2), Hmmr (+2.3), Tlr2 (+1.1); **Tlr4 and Cd14 NOT up** (q>0.2) -> their bulk rise likely infiltration/other cells. Full HA axis at 72h borderline (p=0.057) | Moderate-weak (n=3 vs 5-6) |
| E10 | Specificity of HA | all injury sets | **Not specific**: Provisional matrix (Fn1/Spp1/Thbs/Tnc...) is stronger in every injury set, including microglia | - |

## Honest summary
Supported: CNS has a stabilization-dominated, fragmentation-poor transcript profile (esp. mouse); injury engages the
CD44/TLR2/HMMR receptor arm, partly microglia-intrinsically, in two species and persisting to 8 wk in SCI.
Not supported by these data: increased HA synthesis in injury; clearance, GI, bladder, adipose/lung and spleen
predictions of the tissue-map; HA-specificity (generic ECM remodelling is stronger).
**Implication for the manuscript:** the Figure-2 narrative ("bladder leads synthesis; GI dominates clearance")
did not replicate on this independent array atlas. Either the claim needs softening/qualification (platform, species,
bulk-tissue composition) or it should be re-checked in the source dataset used for Figure 2. Please do not cite E6 as support.

## Deviations / limitations (disclosed)
* Human atlas CNS list (19 region labels) was assembled at run time after viewing sample labels only (no expression);
  the mouse list was fixed in the pre-registration. Human GPL96 lacks HAS3, HAPLN3/4 (coverage in output).
* P3 (bladder rank) is reported as a rank hit/miss; its numeric "p" in the CSV is rank/N and not a formal test.
* Permutation p-values bottom out at 1e-4 (10,000 permutations); tissue-group tests treat tissues as exchangeable.
* Injury cohorts: n=3-4/group, few cohorts (2 SCI, 1 TBI, 1 sorted-microglia). Competitive null treats genes as
  independent although axis genes are co-regulated, so p-values are optimistic; rely on cross-cohort consistency.
* Atlas is microarray (GPL1261/GPL96); the review's Figure 2 uses TPM RNA-seq. mRNA != HA flux or MW.
* Atlas tissues are bulk: stromal/immune/epithelial composition differences can drive module scores.

## Reproduce
`python3 geo_validation/run_geo_validation.py` (injury cohorts) · `run_tissue_atlas.py` · `run_microglia.py` ·
`make_figures.py`, `make_atlas_figure.py`. Outputs in `geo_validation/results/`. Raw downloads are git-ignored.

## Not yet done (candidates, verified to exist)
Disease cohorts from the review's Table 3 (e.g. IPF lung GSE47460/GSE32537 for HYAL/CD44/TLR; NAFLD liver GSE49541; OA),
checkpoint-readout analysis (Treg/DC-maturation genes vs. HA-receptor arm over time), GSE80174 (rat LFPI, chronic, CD44),
stroke (GSE32529).
