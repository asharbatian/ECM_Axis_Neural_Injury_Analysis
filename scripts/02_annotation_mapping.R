# =============================================================================
# HA Axis Validation Study — Script 02: Ortholog Mapping and Gene Annotation
# -----------------------------------------------------------------------------
# Purpose: Map mouse ECM axis genes to rat orthologs and annotate Clariom S probes.
# Inputs:  Mouse gene sets (manuscript Table 1 or equivalent in-script definitions)
# Outputs: Rat ortholog gene sets and probe-to-gene mapping (see docs/02_annotation_mapping.md)
#
# Author: Dr.-Ing Kevin Joseph
# Group Leader - Laboratory of NeuroEngineering
# Department of Neurosurgery
# Medical Center - University of Freiburg
# =============================================================================

library(tidyverse)
library(biomaRt)
library(clariomsrattranscriptcluster.db)
library(limma)

source("scripts/config.R")

# -----------------------------------------------------------------------------
# 1. Define ECM Axes (Mouse Gene Symbols from Manuscript)
# -----------------------------------------------------------------------------

cat("Defining ECM axis gene sets (mouse)...\n")

ecm_axes_mouse <- list(
  
  Hyaluronan = c(
    # Synthesis
    "Has1", "Has2", "Has3",
    # Degradation
    "Hyal1", "Hyal2", "Hyal3", "Cemip", "Tmem2", "Spam1",
    # Receptors
    "Cd44", "Hmmr", "Tlr2", "Tlr4", "Cd14"
  ),
  
  Provisional_Matrix = c(
    "Fn1", "Eda", "Tnc", "Itgav", "Itgb1", "Itgb3",
    "Itga5", "Spp1", "Thbs1", "Thbs2"
  ),
  
  PNN_CSPG = c(
    "Acan", "Vcan", "Ncan", "Bcan", "Cspg4", "Cspg5",
    "Tnr", "Hapln1", "Hapln4", "Ptprz1"
  ),
  
  Basement_Membrane = c(
    "Lama1", "Lama2", "Lama4", "Lama5", "Lamb1", "Lamb2", "Lamc1",
    "Col4a1", "Col4a2", "Col4a3", "Nid1", "Nid2", "Hspg2", "Agrn"
  ),
  
  Proteases_Regulators = c(
    "Mmp2", "Mmp3", "Mmp9", "Mmp12", "Mmp14",
    "Adamts1", "Adamts4", "Adamts5", "Adamts9",
    "Timp1", "Timp2", "Timp3"
  ),
  
  Crosslinking_Fibrosis = c(
    "Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4",
    "Tgm1", "Tgm2",
    "Col1a1", "Col1a2", "Col3a1", "Col5a1", "Col6a1",
    "Acta2", "Tagln"
  )
)

# Print summary
cat("\nMouse ECM axis sizes:\n")
for (axis in names(ecm_axes_mouse)) {
  cat(sprintf("  %s: %d genes\n", axis, length(ecm_axes_mouse[[axis]])))
}

# -----------------------------------------------------------------------------
# 2. Map Mouse to Rat Orthologs
# -----------------------------------------------------------------------------

cat("\nMapping mouse genes to rat orthologs...\n")

# Get all unique mouse genes
all_mouse_genes <- unique(unlist(ecm_axes_mouse))

# For mouse/rat, gene symbols are typically identical or very similar
# Use direct symbol matching (robust, doesn't depend on Ensembl server)
cat("  Using direct symbol matching (mouse/rat symbols are typically identical)...\n")

orthologs <- data.frame(
  mouse_symbol = all_mouse_genes,
  rat_symbol = all_mouse_genes  # Same symbol for most orthologs
)

cat(sprintf("  Mapped %d genes\n", nrow(orthologs)))

# Create rat gene sets
ecm_axes_rat <- lapply(ecm_axes_mouse, function(mouse_genes) {
  mapped <- orthologs %>%
    filter(mouse_symbol %in% mouse_genes) %>%
    pull(rat_symbol) %>%
    unique()
  return(mapped)
})

# Print mapping summary
cat("\nOrtholog mapping coverage:\n")
mapping_summary <- data.frame(
  axis = names(ecm_axes_mouse),
  mouse_genes = sapply(ecm_axes_mouse, length),
  rat_genes = sapply(ecm_axes_rat, length)
)
mapping_summary$coverage <- round(mapping_summary$rat_genes / mapping_summary$mouse_genes * 100, 1)
print(mapping_summary)

# Identify unmapped genes
unmapped <- setdiff(all_mouse_genes, orthologs$mouse_symbol)
if (length(unmapped) > 0) {
  cat("\nUnmapped mouse genes:\n")
  print(unmapped)
}

# -----------------------------------------------------------------------------
# 3. Annotate Clariom S Probes
# -----------------------------------------------------------------------------

cat("\nAnnotating Clariom S probes...\n")

# Load expression matrix to get probe IDs
expr_matrix <- readRDS(file.path(RESULTS_DIR, "expr_matrix_normalized.rds"))
probe_ids <- rownames(expr_matrix)

# Get probe annotations from Bioconductor annotation package
probe_annotation <- AnnotationDbi::select(
  clariomsrattranscriptcluster.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME", "ENTREZID"),
  keytype = "PROBEID"
)

# Remove duplicates (keep first)
probe_annotation <- probe_annotation %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(PROBEID, .keep_all = TRUE)

cat(sprintf("  %d probes with gene symbols (of %d total)\n",
            nrow(probe_annotation), length(probe_ids)))

# -----------------------------------------------------------------------------
# 4. Create Gene-Level Expression Matrix
# -----------------------------------------------------------------------------

cat("\nCreating gene-level expression matrix...\n")

# Map probes to genes
probe_to_gene <- setNames(probe_annotation$SYMBOL, probe_annotation$PROBEID)

# Filter expression matrix to annotated probes
expr_annotated <- expr_matrix[names(probe_to_gene), ]

# Collapse to gene level (use avereps from limma)
gene_symbols <- probe_to_gene[rownames(expr_annotated)]
expr_gene <- avereps(expr_annotated, ID = gene_symbols)

cat(sprintf("  Gene-level matrix: %d genes x %d samples\n",
            nrow(expr_gene), ncol(expr_gene)))

# Check axis gene coverage in expression data
cat("\nECM axis gene coverage in expression data:\n")
expr_genes <- rownames(expr_gene)
for (axis in names(ecm_axes_rat)) {
  axis_genes <- ecm_axes_rat[[axis]]
  found <- sum(axis_genes %in% expr_genes)
  cat(sprintf("  %s: %d of %d genes found (%.1f%%)\n",
              axis, found, length(axis_genes), found/length(axis_genes)*100))
}

# -----------------------------------------------------------------------------
# 5. Save Outputs
# -----------------------------------------------------------------------------

cat("\nSaving outputs...\n")

ref_dir <- file.path(RESULTS_DIR, "reference")
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)

# Save gene sets
saveRDS(ecm_axes_mouse, file.path(ref_dir, "ecm_axes_mouse.rds"))
saveRDS(ecm_axes_rat, file.path(ref_dir, "ecm_axes_rat.rds"))

# Save ortholog mapping
write_csv(orthologs, file.path(ref_dir, "mouse_rat_orthologs.csv"))
write_csv(mapping_summary, file.path(ref_dir, "ortholog_mapping_summary.csv"))

# Save probe annotations
write_csv(probe_annotation, file.path(ref_dir, "probe_annotations.csv"))

# Save gene-level expression
saveRDS(expr_gene, file.path(RESULTS_DIR, "expr_gene_level.rds"))

cat("\nScript 02 complete!\n")
cat(sprintf("  Rat ECM axes: %s\n", file.path(ref_dir, "ecm_axes_rat.rds")))
cat(sprintf("  Gene expression: %s\n", file.path(RESULTS_DIR, "expr_gene_level.rds")))
