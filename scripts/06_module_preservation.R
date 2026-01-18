# =============================================================================
# Script 06: Module Preservation Analysis
# =============================================================================
#
# Purpose: Test preservation of BI/SCI WGCNA modules in implant data
#
# Input:  Gene expression, module definitions from manuscript
# Output: Module preservation statistics, eigengene correlations
#
# =============================================================================

library(tidyverse)
library(WGCNA)
library(RColorBrewer)

source("analysis/config.R")
source("analysis/R/theme_publication.R")

# Enable multi-threading for WGCNA
allowWGCNAThreads()

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------

cat("Loading data...\n")

expr_gene <- readRDS(file.path(RESULTS_DIR, "expr_gene_level.rds"))
metadata <- read_csv(file.path(RESULTS_DIR, "sample_metadata.csv"), show_col_types = FALSE)

# Filter to timepoint samples
metadata <- metadata %>% filter(timepoint != "Control")
expr_gene <- expr_gene[, metadata$sample_id]

# -----------------------------------------------------------------------------
# 2. Define Reference Modules (from Manuscript)
# -----------------------------------------------------------------------------

cat("\nDefining reference modules from manuscript...\n")

# These are approximate module gene signatures based on manuscript GO enrichments
# In practice, these would come from the original WGCNA analysis

bi_modules <- list(
  
  # M1: Neuronal/metabolic - downregulated acutely
  BI_M1_neuronal = c(
    "Slc1a2", "Slc1a3", "Gria2", "Grin1", "Grin2a", "Grin2b",
    "Camk2a", "Camk2b", "Syn1", "Syp", "Snap25", "Vamp2",
    "Atp1a1", "Atp1a2", "Atp1a3", "Kcnj10", "Kcnq2"
  ),
  
  # M2: Acute inflammatory - TLR/NF-kB
  BI_M2_inflammatory = c(
    "Tlr2", "Tlr4", "Cd14", "Myd88", "Nfkb1", "Rela",
    "Il1b", "Il6", "Tnf", "Ccl2", "Cxcl1", "Cxcl2",
    "Ptgs2", "Nos2", "Hmox1"
  ),
  
  # M3: ECM remodeling
  BI_M3_ecm = c(
    "Fn1", "Tnc", "Vcan", "Has2", "Cd44", "Spp1",
    "Mmp9", "Mmp12", "Timp1", "Serpine1",
    "Tgfb1", "Col1a1", "Col3a1"
  ),
  
  # M4: Protein synthesis/ER stress
  BI_M4_protein = c(
    "Hspa5", "Hsp90b1", "Calr", "Canx", "Pdia3", "Pdia4",
    "Ddit3", "Atf4", "Xbp1", "Eif2ak3",
    "Rpl3", "Rpl4", "Rps3", "Rps6"
  ),
  
  # M5: Vesicle/lysosome
  BI_M5_vesicle = c(
    "Lamp1", "Lamp2", "Ctsd", "Ctsl", "Ctsb",
    "Lc3b", "Becn1", "Atg5", "Sqstm1",
    "Rab5a", "Rab7a", "Vps35"
  ),
  
  # M6: Metabolic recovery
  BI_M6_metabolic = c(
    "Ppargc1a", "Tfam", "Nrf1", "Cox4i1", "Cox5a",
    "Atp5a1", "Atp5b", "Ndufv1", "Sdhb",
    "Slc2a1", "Hk1", "Pfkm", "Ldha"
  )
)

# SCI modules (similar structure)
sci_modules <- list(
  
  SCI_M2_inflammatory = c(
    "Tlr2", "Tlr4", "Cd14", "Il1b", "Il6", "Tnf",
    "Ccl2", "Ccl3", "Cxcl1", "Cxcl10",
    "Nfkb1", "Stat3", "Socs3"
  ),
  
  SCI_M4_ecm = c(
    "Fn1", "Tnc", "Vcan", "Has2", "Cd44",
    "Col1a1", "Col3a1", "Col4a1",
    "Mmp2", "Mmp9", "Timp1"
  ),
  
  SCI_M5_cytokine = c(
    "Il1b", "Il6", "Tnf", "Ifng", "Il10",
    "Tgfb1", "Csf1", "Csf2",
    "Ccl2", "Ccl5", "Cxcl12"
  ),
  
  SCI_M6_repair = c(
    "Vim", "Gfap", "Nes", "Sox2", "Sox9",
    "Bdnf", "Ntf3", "Gdnf", "Cntf",
    "Igf1", "Vegfa", "Fgf2"
  )
)

# Combine all modules
all_modules <- c(bi_modules, sci_modules)

cat(sprintf("Defined %d reference modules\n", length(all_modules)))

# -----------------------------------------------------------------------------
# 3. Calculate Module Eigengenes
# -----------------------------------------------------------------------------

cat("\nCalculating module eigengenes in implant data...\n")

# Function to calculate module eigengene (PC1)
calc_eigengene <- function(expr, module_genes) {
  # Filter to genes present in data
  genes_present <- intersect(module_genes, rownames(expr))
  
  if (length(genes_present) < 3) {
    return(rep(NA, ncol(expr)))
  }
  
  module_expr <- expr[genes_present, ]
  
  # Scale and compute PC1
  scaled_expr <- t(scale(t(module_expr)))
  pca <- prcomp(t(scaled_expr), scale = FALSE)
  
  # Return PC1 as eigengene
  return(pca$x[, 1])
}

# Calculate eigengenes for all modules
module_eigengenes <- sapply(all_modules, function(genes) {
  calc_eigengene(expr_gene, genes)
})

# Report coverage
cat("\nModule gene coverage in expression data:\n")
for (mod in names(all_modules)) {
  genes <- all_modules[[mod]]
  found <- sum(genes %in% rownames(expr_gene))
  cat(sprintf("  %s: %d/%d genes (%.0f%%)\n", mod, found, length(genes), found/length(genes)*100))
}

# -----------------------------------------------------------------------------
# 4. Module Activity Analysis
# -----------------------------------------------------------------------------

cat("\nAnalyzing module activity...\n")

# Prepare data
me_df <- as.data.frame(module_eigengenes)
me_df$sample_id <- colnames(expr_gene)
me_df <- me_df %>%
  left_join(metadata, by = "sample_id")

# Test module activity: Implant vs Control at each timepoint
library(limma)

timepoints <- c("Week0", "Week1", "Week2", "Week4", "Week18")
module_stats <- data.frame()

for (mod in names(all_modules)) {
  for (tp in timepoints) {
    tp_data <- me_df %>% filter(timepoint == tp)
    
    if (nrow(tp_data) < 4) next
    
    # Get matched pairs (same replicate number)
    implant_data <- tp_data %>% filter(condition == "Implant") %>% arrange(replicate)
    control_data <- tp_data %>% filter(condition == "Control") %>% arrange(replicate)
    
    # Find common replicates
    common_reps <- intersect(implant_data$replicate, control_data$replicate)
    
    if (length(common_reps) < 2) next
    
    implant <- implant_data[[mod]][implant_data$replicate %in% common_reps]
    control <- control_data[[mod]][control_data$replicate %in% common_reps]
    
    if (length(implant) != length(control)) next
    if (any(is.na(implant)) || any(is.na(control))) next
    
    # Paired t-test
    test <- t.test(implant, control, paired = TRUE)
    
    row <- data.frame(
      module = mod,
      timepoint = tp,
      day = TIMEPOINT_MAP[tp],
      mean_implant = mean(implant),
      mean_control = mean(control),
      delta = mean(implant) - mean(control),
      p_value = test$p.value
    )
    module_stats <- rbind(module_stats, row)
  }
}

# Adjust p-values
module_stats$adj_p_value <- p.adjust(module_stats$p_value, method = "BH")

# Save results
preserve_dir <- file.path(RESULTS_DIR, "preservation")
dir.create(preserve_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(module_stats, file.path(preserve_dir, "module_activity_stats.csv"))
saveRDS(module_eigengenes, file.path(preserve_dir, "module_eigengenes.rds"))

# -----------------------------------------------------------------------------
# 5. Preservation Summary
# -----------------------------------------------------------------------------

cat("\nModule preservation summary:\n")

# Count significant activations per module
module_summary <- module_stats %>%
  group_by(module) %>%
  summarize(
    n_timepoints_sig = sum(adj_p_value < 0.05),
    mean_delta = mean(delta),
    max_delta = max(abs(delta)),
    peak_timepoint = timepoint[which.max(abs(delta))]
  ) %>%
  arrange(desc(n_timepoints_sig))

print(module_summary)

write_csv(module_summary, file.path(preserve_dir, "module_preservation_summary.csv"))

# -----------------------------------------------------------------------------
# 6. Visualizations
# -----------------------------------------------------------------------------

cat("\nGenerating visualizations...\n")

# 6a. Module activity heatmap
delta_matrix <- module_stats %>%
  select(module, timepoint, delta) %>%
  pivot_wider(names_from = timepoint, values_from = delta) %>%
  column_to_rownames("module")

library(pheatmap)

pdf(file.path(FIGURES_DIR, "06_module_activity_heatmap.pdf"), width = 8, height = 8)
pheatmap(
  as.matrix(delta_matrix),
  scale = "none",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  border_color = "grey80",
  main = "Module Activity (Implant - Control)",
  angle_col = 45
)
dev.off()

# 6b. Key module trajectories
key_modules <- c("BI_M2_inflammatory", "BI_M3_ecm", "BI_M6_metabolic", "SCI_M4_ecm")

me_long <- me_df %>%
  pivot_longer(cols = all_of(names(all_modules)), names_to = "module", values_to = "eigengene") %>%
  filter(module %in% key_modules)

p_modules <- ggplot(me_long, aes(x = TIMEPOINT_MAP[timepoint], y = eigengene, 
                                   color = condition, group = condition)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
  stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.8, linewidth = 0.8) +
  facet_wrap(~module, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = c(0, 7, 14, 28, 126), trans = "pseudo_log") +
  scale_color_manual(values = colors_condition) +
  labs(
    title = "A",
    subtitle = "Reference Module Activity in Implant Data",
    x = "Days Post-Implantation",
    y = "Module Eigengene"
  ) +
  theme_publication(base_size = 10) +
  theme(legend.position = "bottom")

save_publication_figure(p_modules, file.path(FIGURES_DIR, "06_key_module_trajectories"), width = 9, height = 7)

# 6c. Preservation bar plot
p_preserve <- ggplot(module_summary, aes(x = reorder(module, n_timepoints_sig), 
                                          y = n_timepoints_sig, fill = mean_delta > 0)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "#1A9850"),
                    labels = c("TRUE" = "Activated", "FALSE" = "Suppressed"),
                    name = "") +
  labs(
    title = "B",
    subtitle = "Module Preservation (timepoints with adj.p < 0.05)",
    x = "",
    y = "Number of Significant Timepoints"
  ) +
  theme_publication(base_size = 10) +
  theme(legend.position = "bottom")

save_publication_figure(p_preserve, file.path(FIGURES_DIR, "06_module_preservation_bar"), type = "barplot")

cat("\nScript 06 complete!\n")
cat(sprintf("  Preservation results: %s\n", preserve_dir))
