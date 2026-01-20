# =============================================================================
# HA Axis Validation Study — Script 07: Pathway Enrichment
# -----------------------------------------------------------------------------
# Purpose: Test enrichment of HA- and ECM-related pathways and validate literature
#          signatures supporting HA as a mechanistic driver.
# Inputs:  DEG results and ECM axis gene sets (from prior scripts)
# Outputs: GSEA/enrichment result tables and plots (see docs/07_pathway_enrichment.md)
#
# Author: Dr.-Ing Kevin Joseph
# Group Leader - Laboratory of NeuroEngineering
# Department of Neurosurgery
# Medical Center - University of Freiburg
# =============================================================================

library(tidyverse)
library(fgsea)
library(clusterProfiler)
library(org.Rn.eg.db)

# Ensure dplyr functions take precedence
select <- dplyr::select
filter <- dplyr::filter

source("scripts/config.R")
source("scripts/theme_publication.R")

ha_dir <- file.path(RESULTS_DIR, "ha_analysis")
dir.create(ha_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------

cat("=== Script 07: Pathway Enrichment ===\n\n")
cat("1. Loading data...\n")

# Load DEG results
deg_results <- readRDS(file.path(RESULTS_DIR, "deg/deg_results_list.rds"))

# Load ECM axes
ecm_axes_rat <- readRDS(file.path(RESULTS_DIR, "reference/ecm_axes_rat.rds"))

cat("  Timepoints with DEG results:", paste(names(deg_results), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 2. Define Gene Sets for GSEA
# -----------------------------------------------------------------------------

cat("\n2. Defining gene sets for enrichment analysis...\n")

# Literature implant signatures
implant_signatures <- list(
  # Huff et al. - chronic implant markers
  Huff_Chronic = c("Hmox1", "Il6", "Ccl2", "Cd14", "Cd44", "Socs3", 
                   "Tlr2", "Tlr4", "Gfap", "Vim"),
  
  # Joseph et al. - acute MMP response
  Joseph_Acute = c("Mmp12", "Mmp9", "Mmp14", "Cd14", "Cd44", 
                   "Lcn2", "Gfap", "Vim", "Timp1"),
  
  # HA-specific signaling
  HA_DAMP_Signaling = c("Cd44", "Tlr2", "Tlr4", "Cd14", "Myd88", "Nfkb1", 
                        "Il1b", "Il6", "Tnf", "Ccl2"),
  
  # HA synthesis/degradation
  HA_Metabolism = c("Has1", "Has2", "Has3", "Hyal1", "Hyal2", "Cemip", 
                    "Tmem2", "Cd44", "Hmmr")
)

# Add ECM axes as gene sets
all_gene_sets <- c(implant_signatures, ecm_axes_rat)

# Convert to uppercase for matching
all_gene_sets <- lapply(all_gene_sets, toupper)

cat("  Gene sets defined:\n")
for (gs in names(all_gene_sets)) {
  cat(sprintf("    %s: %d genes\n", gs, length(all_gene_sets[[gs]])))
}

# -----------------------------------------------------------------------------
# 3. Run GSEA for Each Timepoint
# -----------------------------------------------------------------------------

cat("\n3. Running Gene Set Enrichment Analysis...\n")

gsea_results <- list()

for (tp in names(deg_results)) {
  cat(sprintf("  Processing %s...\n", tp))
  
  # Get ranked gene list (by t-statistic)
  deg_df <- deg_results[[tp]]
  
  # Ensure we have the t-statistic
  if (!"t" %in% colnames(deg_df)) {
    # Use logFC * -log10(p) as proxy
    deg_df$t <- deg_df$logFC * -log10(deg_df$P.Value + 1e-300)
  }
  
  # Create named vector, uppercase
  ranked_genes <- setNames(deg_df$t, toupper(deg_df$gene))
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  ranked_genes <- ranked_genes[!is.na(ranked_genes)]
  
  # Run fgsea
  fgsea_res <- fgsea(
    pathways = all_gene_sets,
    stats = ranked_genes,
    minSize = 3,
    maxSize = 500,
    eps = 0
  )
  
  fgsea_res$timepoint <- tp
  gsea_results[[tp]] <- fgsea_res
}

# Combine results
gsea_all <- bind_rows(gsea_results)

cat("\n  GSEA complete.\n")

# -----------------------------------------------------------------------------
# 4. Analyze HA-Related Pathway Enrichment
# -----------------------------------------------------------------------------

cat("\n4. Analyzing HA-related pathway enrichment...\n")

# Focus on HA-related pathways
ha_pathways <- c("HA_DAMP_Signaling", "HA_Metabolism", "Hyaluronan")

ha_enrichment <- gsea_all %>%
  filter(pathway %in% toupper(ha_pathways) | pathway %in% ha_pathways) %>%
  arrange(timepoint, pval)

cat("\n  HA Pathway Enrichment:\n")
print(ha_enrichment %>% select(pathway, timepoint, NES, pval, padj))

# Check if HA pathways are enriched in acute phase
acute_ha <- ha_enrichment %>% filter(timepoint %in% c("Week0", "Week1"))
if (nrow(acute_ha) > 0 && any(acute_ha$padj < 0.1)) {
  cat("\n  --> HA pathways ENRICHED in acute phase (supports manuscript hypothesis)\n")
} else {
  cat("\n  --> HA pathway enrichment not significant in acute phase\n")
}

# -----------------------------------------------------------------------------
# 5. Compare Literature Signatures
# -----------------------------------------------------------------------------

cat("\n5. Comparing with literature implant signatures...\n")

lit_enrichment <- gsea_all %>%
  filter(pathway %in% c("HUFF_CHRONIC", "JOSEPH_ACUTE", "Huff_Chronic", "Joseph_Acute"))

if (nrow(lit_enrichment) > 0) {
  cat("\n  Literature Signature Enrichment:\n")
  print(lit_enrichment %>% select(pathway, timepoint, NES, padj) %>% arrange(timepoint))
  
  # Summarize
  sig_huff <- lit_enrichment %>% 
    filter(grepl("HUFF|Huff", pathway), padj < 0.05) %>% 
    nrow()
  sig_joseph <- lit_enrichment %>% 
    filter(grepl("JOSEPH|Joseph", pathway), padj < 0.05) %>% 
    nrow()
  
  cat(sprintf("\n  Huff signature significant at %d timepoints\n", sig_huff))
  cat(sprintf("  Joseph signature significant at %d timepoints\n", sig_joseph))
}

# -----------------------------------------------------------------------------
# 6. GO Enrichment for Peak Response
# -----------------------------------------------------------------------------

cat("\n6. Running GO enrichment for peak response timepoint...\n")

# Find peak DEG timepoint
deg_counts <- sapply(deg_results, function(x) sum(x$adj.P.Val < 0.05 & abs(x$logFC) > 0.5))
peak_tp <- names(which.max(deg_counts))

cat(sprintf("  Peak response at %s (%d significant DEGs)\n", peak_tp, max(deg_counts)))

# Get significant DEGs
peak_deg <- deg_results[[peak_tp]] %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 0.5)

# Convert to Entrez IDs
gene_symbols <- peak_deg$gene
gene_map <- tryCatch({
  bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)
}, error = function(e) {
  cat("  Warning: Could not map gene symbols to Entrez IDs\n")
  data.frame(SYMBOL = character(), ENTREZID = character())
})

if (nrow(gene_map) > 10) {
  # GO Biological Process enrichment
  go_bp <- enrichGO(
    gene = gene_map$ENTREZID,
    OrgDb = org.Rn.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
    go_results <- as.data.frame(go_bp)
    
    # Filter for HA-related terms
    ha_terms <- c("hyaluronan", "glycosaminoglycan", "toll-like", "pattern recognition",
                  "innate immune", "inflammatory", "extracellular matrix", "wound")
    
    ha_go <- go_results %>%
      filter(str_detect(tolower(Description), paste(ha_terms, collapse = "|")))
    
    if (nrow(ha_go) > 0) {
      cat("\n  HA-Related GO Terms Enriched:\n")
      print(ha_go %>% select(Description, GeneRatio, p.adjust) %>% head(10))
    }
    
    write_csv(go_results, file.path(ha_dir, sprintf("go_enrichment_%s.csv", peak_tp)))
  }
}

# -----------------------------------------------------------------------------
# 7. Visualizations
# -----------------------------------------------------------------------------

cat("\n7. Generating visualizations...\n")

# 7a. GSEA NES heatmap across timepoints
gsea_matrix <- gsea_all %>%
  select(pathway, timepoint, NES) %>%
  pivot_wider(names_from = timepoint, values_from = NES) %>%
  column_to_rownames("pathway")

# Reorder columns by timepoint
tp_order <- intersect(TIMEPOINTS, colnames(gsea_matrix))
gsea_matrix <- gsea_matrix[, tp_order, drop = FALSE]

# Filter to interesting pathways
interesting <- c("HA_DAMP_Signaling", "HA_Metabolism", "Hyaluronan", 
                 "HUFF_CHRONIC", "JOSEPH_ACUTE",
                 "Provisional_Matrix", "CSPG_PNN", "Basement_Membrane",
                 "Proteases_Regulators", "Crosslinking_Fibrosis")
interesting <- intersect(interesting, rownames(gsea_matrix))

if (length(interesting) > 0) {
  gsea_plot_data <- gsea_matrix[interesting, , drop = FALSE]
  
  p_gsea_heatmap <- gsea_plot_data %>%
    rownames_to_column("pathway") %>%
    pivot_longer(-pathway, names_to = "timepoint", values_to = "NES") %>%
    mutate(timepoint = factor(timepoint, levels = TIMEPOINTS)) %>%
    ggplot(aes(x = timepoint, y = pathway, fill = NES)) +
    geom_tile(color = "grey80", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f", NES)), size = 3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                         midpoint = 0, na.value = "grey90") +
    labs(
      title = "Gene Set Enrichment Across Timepoints",
      subtitle = "NES: Normalized Enrichment Score (positive = upregulated in implant)",
      x = "", y = ""
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  save_publication_figure(p_gsea_heatmap, 
                          file.path(FIGURES_DIR, "individual", "07_gsea_heatmap"),
                          width = 9, height = 7)
}

# 7b. HA pathway trajectory
ha_traj <- gsea_all %>%
  filter(pathway %in% c("HA_DAMP_Signaling", "HA_Metabolism", "Hyaluronan", "HYALURONAN")) %>%
  mutate(
    timepoint = factor(timepoint, levels = TIMEPOINTS),
    day = case_when(
      timepoint == "Week0" ~ 0,
      timepoint == "Week1" ~ 7,
      timepoint == "Week2" ~ 14,
      timepoint == "Week4" ~ 28,
      timepoint == "Week18" ~ 126
    ),
    significant = padj < 0.05
  )

if (nrow(ha_traj) > 0) {
  p_ha_traj <- ggplot(ha_traj, aes(x = day, y = NES, color = pathway)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = significant), size = 3) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
    scale_x_continuous(breaks = c(0, 7, 14, 28, 126), trans = "pseudo_log") +
    scale_color_brewer(palette = "Set2") +
    labs(
      title = "HA Pathway Enrichment Trajectory",
      subtitle = "Filled points = FDR < 0.05",
      x = "Days Post-Implantation",
      y = "Normalized Enrichment Score"
    ) +
    theme_publication() +
    theme(legend.position = "bottom")
  
  save_publication_figure(p_ha_traj, 
                          file.path(FIGURES_DIR, "individual", "07_ha_pathway_trajectory"),
                          width = 8, height = 5)
}

# -----------------------------------------------------------------------------
# 8. Save Results
# -----------------------------------------------------------------------------

cat("\n8. Saving results...\n")

upstream_results <- list(
  gsea_all = gsea_all,
  ha_enrichment = ha_enrichment,
  lit_enrichment = lit_enrichment,
  peak_timepoint = peak_tp,
  gene_sets_used = names(all_gene_sets)
)

saveRDS(upstream_results, file.path(ha_dir, "upstream_regulator_results.rds"))
write_csv(gsea_all, file.path(ha_dir, "gsea_all_results.csv"))

cat("\n=== Script 07 Complete ===\n")

# Print key summary
cat("\nKEY FINDINGS:\n")
sig_ha <- gsea_all %>% 
  filter(grepl("HA|Hyaluronan", pathway, ignore.case = TRUE), padj < 0.05)
if (nrow(sig_ha) > 0) {
  cat("  HA-related pathways significantly enriched at:\n")
  print(sig_ha %>% select(pathway, timepoint, NES, padj))
} else {
  cat("  No significant HA pathway enrichment detected\n")
}
