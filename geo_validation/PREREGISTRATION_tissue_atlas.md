# Pre-registration: tissue-atlas test of the HA-module claims (Figure 2 of the review)

Written and committed BEFORE any expression value from GSE10246 / GSE1133 was inspected
(only sample labels were viewed). The git commit timestamp is the evidence of this.

## Data
* Primary: GSE10246 (mouse GNF GeneAtlas V3, GPL1261). Cell lines, sorted immune cells, progenitors,
  embryonic/stem lines, LPS/IgE-stimulated and osteoblast-differentiation samples are EXCLUDED
  (tissues only). Replicates averaged per tissue.
* Secondary (human, CNS/liver/lung only - atlas lacks bladder, spleen, GI): GSE1133 GPL96.

## Modules (as defined in the review text; mouse symbols)
* Synthesis: Has1, Has2, Has3
* Stabilization: Hapln1, Hapln2, Hapln3, Hapln4
* Fragmentation: Hyal1, Hyal2, Hyal4
* Clearance: Cemip, Tmem2, Spam1
Module score of a tissue = mean over genes of the gene's z-score across tissues. Genes absent from the
platform are dropped and reported.

## Tissue groups (fixed now)
* CNS: amygdala, cerebellum, cerebral_cortex, cerebral_cortex_prefrontal, dorsal_striatum, hippocampus,
  hypothalamus, nucleus_accumbens, olfactory_bulb, spinal_cord
* GI: stomach, intestine_small, intestine_large
* Fragmentation organs: spleen, liver
* Adipose: adipose_white, adipose_brown
* Single: bladder, lung

## Predictions (directional) and test
Statistic: mean percentile rank of the group (1 = highest) among all tissues for the module;
p = fraction of 100,000 random same-size tissue groups with an equal-or-more-extreme mean rank.
* P1  CNS is HIGH for Stabilization.
* P2a CNS is LOW for Fragmentation.   P2b CNS is LOW for Clearance.
* P3  Bladder has the top Synthesis score (rank 1 of N), secondary: rank <= 3.
* P4  Spleen+liver are HIGH for Fragmentation.
* P5  GI group is HIGH for Clearance.
* P6  (secondary) Adipose and lung rank 2nd-3rd for Synthesis, i.e. within the top 5.
* P7  (the "CNS-unique" claim) CNS rank on (Stabilization - Fragmentation) is HIGH.
Multiple testing: BH across P1-P7 (primary: P1-P5, P7). Failure of a prediction will be reported as a failure.

## Interpretation limits (committed in advance)
Array data from a different platform/species than the review's TPM figure; mRNA is a proxy for pathway
capacity, not HA flux or molecular weight. Agreement supports the transcript-level claim only.
