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

---
# Part 2 - HA-fragmentation machinery vs innate DAMP response across diseases
Pre-registered in `PREREGISTRATION_damp_diseases.md` (commit c919b2b, before any disease-cohort expression was inspected).
Code `run_damp_diseases.py`; figure `results/fig4_damp_coupling.png`; tables `results/damp_*.csv`.
Cohorts: IPF lung GSE53845, NASH liver GSE48452, OA synovium GSE55235 (RA reported separately), carotid atheroma GSE43292 (paired),
psoriasis skin GSE13355, Alzheimer brain GSE5281 (region-demeaned), bladder cancer GSE13507 + SCI/TBI (3 cohorts) = 10 primary.

| # | Prediction | Result | Verdict |
|---|---|---|---|
| D1a | DAMP innate program up in disease | Stouffer z=5.6 (8/10 cohorts positive); **8.9 if the synovium cohort is excluded** (see caveat) | **Supported**, strongest in injury and atheroma |
| D1b | HA-receptor arm (SENSE) up | z=5.6 (8/10) | **Supported**, mostly injury cohorts |
| D1c | FRAG (HYAL1/2, CEMIP, TMEM2) up | z=2.1 (9/10 z>0), but **IPF is significantly DOWN (z=-2.5)** | **Weak / partly supported**; effect size small |
| D1d | SYNTH (HAS) | z=-1.1 (3/8 positive) | Not up - consistent with Part 1 |
| D2 | FRAG co-varies with DAMP beyond leukocyte content and random gene sets | pooled partial r = **0.32 [0.22, 0.41]**, I2=27%; Stouffer z vs expression-matched null = 3.9; leave-one-out range 0.28-0.35; human-only (7) r=0.31; but only **4/10 cohorts individually p<0.05** | **Supported (modest)** |
| D2b | SENSE~DAMP | r=0.32 [0.20,0.42], I2=41%; null-z=2.9; only 2/10 individually significant | Supported (modest) |
| D4 | Innate specificity: ADAPT (T/B lymphocyte) genes should NOT couple to DAMP | ADAPT~DAMP r=-0.01 [-0.30,0.28], null-z=-1.6 | **Supported**: coupling is specific to the innate arm |
| D4b | SYNTH negative control | r=0.20 [0.03,0.36], I2=73%, null-z=1.1 (n.s.) | Weak positive, not clean negative control |
| D5 | Cross-disease: modules with larger FRAG effect have larger DAMP effect | Spearman 0.82 (SENSE 0.93) over 10 contrasts; **human-only n=7: 0.57 (SENSE 0.82)** | Exploratory only (effect sizes across studies are not comparable; n small) |

## What this says about "HA fragments -> DAMP -> innate immunity"
Supported at the transcript level: in several unrelated diseases, expression of the HA-fragmenting enzymes and the HA-fragment
receptors rises with the innate TLR/MyD88/NF-kB program, this coupling is not explained by leukocyte content (LEUK-adjusted),
exceeds random expression-matched gene sets, and is absent for adaptive-lymphocyte genes. This matches your reading that the
fragment -> DAMP arm is an *innate* response.
Not shown: (i) that fragments exist or have low molecular weight (not measurable here); (ii) direction - TNF/IL-1b strongly
induce CEMIP/HYAL/CD44 in stromal cells, so the coupling is equally compatible with inflammation driving fragmentation
(feed-forward loop) as with fragments driving inflammation; (iii) that HA is special - other stress/NF-kB-inducible genes could behave alike;
(iv) the review's own caveat that endotoxin-free HA fragments may not activate TLRs (Dong 2016) is not addressed by expression data.

## Part-2 caveats (disclosed)
* **GSE55235 (synovium):** labels are correct (MMP1/MMP3 up in OA/RA, CXCL13 only in RA) yet DAMP is strongly LOWER in OA/RA
  (z=-8.9 / -2.7) because immediate-early/stress genes (FOS, JUN, EGR1, NFKBIA, CXCL8, IL6) are higher in the "healthy" synovium,
  probably a sampling/handling artifact. The primary result was kept as pre-specified; do not read this as evidence against DAMP in arthritis.
* DAMP score contains stress/immediate-early-type genes, so tissue-procurement artifacts can bias D1 and D2.
* IPF: FRAG is lower in disease, which does not fit "more fragmentation" (CEMIP/HYAL transcripts need gene-level follow-up).
* RA shares controls with OA, so RA is reported but excluded from the pooled estimate.
* Only partial-r and null-z summarize D2; per-cohort significance is limited (n=10-165 and a wide null in small cohorts: injury n=6-8).
* Expression is bulk tissue: cell-type composition beyond the LEUK adjustment (fibroblasts, endothelium) can still drive the coupling.
* Mouse/rat symbols mapped by upper-cased gene symbol; platform differences (Affymetrix, Agilent, Illumina) are not harmonised.

---
# Part 3 - per-gene breakdown of the receptor arm (CD44, CD14, TLR2, TLR4) - EXPLORATORY
Decomposition of the pre-registered SENSE module; the gene-level split itself was NOT pre-registered. `run_damp_receptors.py`,
`results/damp_receptor_genes.csv`, `results/damp_receptor_genes_pooled.csv`. Same 10 primary cohorts.

| Gene | Cohorts up (q<0.05) / down | Pooled partial r with DAMP given leukocytes [95% CI] | I2 | null-z |
|---|---|---|---|---|
| CD44 | 3 / **2** (down in psoriasis, bladder cancer) | **0.27 [0.12, 0.41]** | 64% | 3.5 |
| TLR2 | 4 / **0** (never down) | **0.23 [0.11, 0.35]** | 41% | 2.0 |
| TLR4 | 4 / 2 | 0.15 [0.02, 0.28] | 50% | 1.1 |
| CD14 | 5 / 1 | 0.13 [-0.04, 0.30] (CI includes 0) | 70% | 1.2 |

Reading: the module-level signal is carried mainly by **CD44 and TLR2**; TLR4 and CD14 are weak/inconsistent. This agrees with the
independent sorted-microglia result (E9: Cd44, Tlr2 up in microglia; Tlr4, Cd14 not). Caveats: CD44 and TLR2 are themselves NF-kB-inducible,
so coupling with DAMP can mean the receptor is induced by inflammation rather than that HA fragments signal through it; all four bind
many ligands (CD44: osteopontin etc.; TLR2/4: many DAMPs/PAMPs), so this is not HA-specific; heterogeneity is high (I2 41-70%);
rat SCI shows negative CD14/TLR4 couplings (n=8, wide null).
