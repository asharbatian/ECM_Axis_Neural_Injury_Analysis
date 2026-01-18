# =============================================================================
# Script 03: Differential Expression Analysis
# =============================================================================
# Performs paired DEG analysis (Implant vs Control) at each timepoint using limma
# =============================================================================

cat("=== Script 03: Differential Expression Analysis ===\n")

suppressPackageStartupMessages({
  library(limma)
  library(tidyverse)
})

source("analysis/config.R")
source("analysis/R/theme_publication.R")

# =============================================================================
# 1. Load Data
# =============================================================================
cat("\n1. Loading data...\n")

expr_gene <- readRDS(file.path(RESULTS_DIR, "expr_gene_level.rds"))
metadata <- read.csv(file.path(RESULTS_DIR, "sample_metadata.csv"))

# Exclude baseline samples (no implant comparison possible)
metadata <- metadata %>%
  filter(timepoint != "Baseline") %>%
  mutate(
    timepoint = factor(timepoint, levels = TIMEPOINTS),
    condition = factor(condition, levels = c("Control", "Implant")),
    group = interaction(condition, timepoint, drop = TRUE)
  )

expr_gene <- expr_gene[, metadata$sample_id]

cat("  Expression matrix:", nrow(expr_gene), "genes x", ncol(expr_gene), "samples\n")
cat("  Timepoints:", paste(levels(metadata$timepoint), collapse = ", "), "\n")

# =============================================================================
# 2. Build Design Matrix for Paired Analysis
# =============================================================================
cat("\n2. Setting up paired design matrix...\n")

# Design: ~0 + group + animal_id (blocking on animal for paired design)
design <- model.matrix(~ 0 + group + animal_id, data = metadata)
colnames(design) <- make.names(colnames(design))

# Simplify group names
colnames(design) <- gsub("group", "", colnames(design))

cat("  Design matrix:", nrow(design), "samples x", ncol(design), "coefficients\n")

# =============================================================================
# 3. Fit Linear Model
# =============================================================================
cat("\n3. Fitting linear model...\n")

fit <- lmFit(expr_gene, design)

# Create contrasts: Implant vs Control at each timepoint
contrast_names <- paste0("Implant.", TIMEPOINTS, " - Control.", TIMEPOINTS)
names(contrast_names) <- TIMEPOINTS

# Build contrast matrix dynamically
contrast_list <- list()
for (tp in TIMEPOINTS) {
  implant_col <- paste0("Implant.", tp)
  control_col <- paste0("Control.", tp)
  if (implant_col %in% colnames(design) && control_col %in% colnames(design)) {
    contrast_list[[tp]] <- paste0(implant_col, " - ", control_col)
  }
}

contrast_matrix <- makeContrasts(contrasts = contrast_list, levels = design)
colnames(contrast_matrix) <- names(contrast_list)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

cat("  Contrasts tested:", paste(colnames(contrast_matrix), collapse = ", "), "\n")

# =============================================================================
# 4. Extract Results
# =============================================================================
cat("\n4. Extracting differential expression results...\n")

deg_results <- list()
deg_summary <- data.frame()

for (tp in colnames(contrast_matrix)) {
  res <- topTable(fit2, coef = tp, number = Inf, sort.by = "P") %>%
    rownames_to_column("gene") %>%
    mutate(
      timepoint = tp,
      direction = case_when(
        adj.P.Val < DE_PADJ_THRESHOLD & logFC > DE_LFC_THRESHOLD ~ "Up",
        adj.P.Val < DE_PADJ_THRESHOLD & logFC < -DE_LFC_THRESHOLD ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  deg_results[[tp]] <- res
  
  # Summary stats
  n_up <- sum(res$direction == "Up")
  n_down <- sum(res$direction == "Down")
  
  deg_summary <- bind_rows(deg_summary, data.frame(
    timepoint = tp,
    n_up = n_up,
    n_down = n_down,
    n_total = n_up + n_down
  ))
  
  # Save individual results
  write.csv(res, file.path(RESULTS_DIR, "deg", paste0("deg_", tp, ".csv")), row.names = FALSE)
  
  cat(sprintf("  %s: %d up, %d down (total: %d DEGs)\n", tp, n_up, n_down, n_up + n_down))
}

# =============================================================================
# 5. Create Volcano Plots
# =============================================================================
cat("\n5. Creating volcano plots...\n")

for (tp in names(deg_results)) {
  res <- deg_results[[tp]]
  
  p_volcano <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = colors_direction) +
    geom_hline(yintercept = -log10(DE_PADJ_THRESHOLD), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-DE_LFC_THRESHOLD, DE_LFC_THRESHOLD), linetype = "dashed", color = "grey50") +
    labs(
      title = paste("Implant vs Control:", tp),
      subtitle = paste(deg_summary$n_total[deg_summary$timepoint == tp], "DEGs"),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_publication() +
    theme(legend.position = "bottom")
  
  save_publication_figure(p_volcano, 
                          file.path(FIGURES_DIR, "individual", paste0("03_volcano_", tp)),
                          width = 6, height = 6)
}

# Combined volcano facet plot
all_deg <- bind_rows(deg_results)
all_deg$timepoint <- factor(all_deg$timepoint, levels = TIMEPOINTS)

p_volcano_all <- ggplot(all_deg, aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(values = colors_direction) +
  facet_wrap(~ timepoint, nrow = 1) +
  geom_hline(yintercept = -log10(DE_PADJ_THRESHOLD), linetype = "dashed", color = "grey50") +
  labs(title = "Differential Expression Across Timepoints", x = "Log2 FC", y = "-Log10 Adj. P") +
  theme_publication() +
  theme(legend.position = "bottom")

save_publication_figure(p_volcano_all, 
                        file.path(FIGURES_DIR, "supplementary", "03_volcano_all_timepoints"),
                        width = 14, height = 4)

# =============================================================================
# 6. DEG Trajectory Summary Plot
# =============================================================================
cat("\n6. Creating DEG trajectory plot...\n")

deg_summary$timepoint <- factor(deg_summary$timepoint, levels = TIMEPOINTS)

p_deg_trajectory <- deg_summary %>%
  pivot_longer(cols = c(n_up, n_down), names_to = "direction", values_to = "count") %>%
  mutate(
    count = ifelse(direction == "n_down", -count, count),
    direction = ifelse(direction == "n_up", "Upregulated", "Downregulated")
  ) %>%
  ggplot(aes(x = timepoint, y = count, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Upregulated" = colors_direction["Up"], 
                               "Downregulated" = colors_direction["Down"])) +
  labs(
    title = "DEG Count Trajectory",
    subtitle = "Implant vs Control at each timepoint",
    x = "Timepoint", y = "Number of DEGs", fill = ""
  ) +
  theme_publication() +
  theme(legend.position = "bottom")

save_publication_figure(p_deg_trajectory, 
                        file.path(FIGURES_DIR, "individual", "03_deg_trajectory"),
                        type = "barplot")  # 5 x 7 - taller than wide

# =============================================================================
# 7. Save Results
# =============================================================================
cat("\n7. Saving results...\n")

saveRDS(deg_results, file.path(RESULTS_DIR, "deg", "deg_results_list.rds"))
saveRDS(fit2, file.path(RESULTS_DIR, "deg", "limma_fit.rds"))
write.csv(deg_summary, file.path(RESULTS_DIR, "deg", "deg_summary.csv"), row.names = FALSE)

cat("\n=== Script 03 Complete ===\n")
cat("  Peak DEG response:", deg_summary$timepoint[which.max(deg_summary$n_total)], 
    "with", max(deg_summary$n_total), "DEGs\n")
