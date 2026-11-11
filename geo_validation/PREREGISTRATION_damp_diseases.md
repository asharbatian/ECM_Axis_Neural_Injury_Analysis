# Pre-registration: HA-fragmentation machinery vs. innate DAMP response across diseases

Written and committed before any expression value of the disease cohorts below was inspected
(only sample labels / group structure were viewed). Question: do transcripts of the HA fragmentation
machinery and the HA-fragment receptor arm track the innate TLR/MyD88/NF-kB (DAMP) response across diseases,
beyond what generic leukocyte infiltration or random gene sets would produce?
Transcript data cannot show HA fragments themselves; this tests the *machinery-level coupling* only.

## Modules (human symbols; rodent data matched by upper-cased symbol)
* FRAG   (fragmentation): HYAL1, HYAL2, CEMIP, TMEM2
* SYNTH  (synthesis): HAS1, HAS2, HAS3
* SENSE  (HA-fragment receptors): CD44, TLR2, TLR4, CD14, HMMR
* DAMP   (innate TLR/MyD88/NF-kB effector program; deliberately excludes SENSE and HA-matrix genes TNFAIP6/PTX3):
  MYD88, IRAK4, TRAF6, NFKB1, NFKBIA, TNFAIP3, TNF, IL1B, IL1A, IL6, CXCL8, CCL2, CXCL1, CXCL2, CXCL3, CCL20, PTGS2, NLRP3, ICAM1, SOD2
* LEUK   (pan-leukocyte content control): PTPRC, LCP1, CORO1A, CD53, LAPTM5, TYROBP, FCER1G
* ADAPT  (adaptive lymphocyte comparator): CD3D, CD3E, CD2, LCK, CD8A, CD19, MS4A1, CD79A
Module score of a sample = mean over present genes of the gene's z-score across samples of the cohort.

## Cohorts / primary contrasts (fixed now)
IPF lung GSE53845 (IPF vs control) · NASH liver GSE48452 (NASH vs Control) · synovium GSE55235 (OA vs healthy; RA vs healthy)
· carotid atheroma GSE43292 (plaque vs paired intact) · psoriasis skin GSE13355 (involved vs normal) ·
Alzheimer brain GSE5281 (AD vs normal, region-demeaned) · bladder cancer GSE13507 (primary cancer vs normal mucosa) ·
plus the three injury cohorts already analysed (GSE45006 D3, GSE47681 D3, GSE45997 ipsi vs naive).
Cohorts whose platform lacks gene symbols or that cannot be parsed are dropped and listed, not substituted.

## Predictions and tests
* D1  DAMP and SENSE modules are up in disease vs control (competitive permutation test against 10,000
      expression-matched random gene sets of equal size, same machinery as run_geo_validation.py; one-sided). FRAG up: same test.
      SYNTH: reported without a directional prediction.
* D2  Coupling: after removing the group (disease/control), pair/region effects (residuals), the partial Spearman
      correlation of FRAG with DAMP controlling for LEUK exceeds the null distribution of 2,000 expression-matched random
      4-gene sets (same procedure); success = empirical one-sided p < 0.05. Same for SENSE vs DAMP.
      Also report the raw (non-residualized, no LEUK control) correlation for transparency.
* D3  Cross-cohort: Stouffer z of D1 per module (primary contrasts, k cohorts) and inverse-variance Fisher-z mean of the D2 partial correlations
      (random effects, DerSimonian-Laird) with I^2.
* D4  Specificity: the ADAPT module is analysed identically; the "innate" reading is supported if FRAG/SENSE couple with DAMP
      more strongly than with ADAPT. SYNTH is the negative control (no coupling predicted).
* D5  Ecological: across primary contrasts, Spearman correlation between the module effect of FRAG and of DAMP.
Multiple testing: BH within each of D1/D2 across modules x cohorts. Negative results will be reported as such.
