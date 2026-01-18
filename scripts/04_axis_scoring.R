# =============================================================================
# Script 04: ECM Axis Scoring (GSVA/ssGSEA)
# =============================================================================
#
# Purpose: Calculate per-sample enrichment scores for each ECM axis
#
# Input:  Gene-level expression, ECM axis gene sets
# Output: ssGSEA scores, axis activation statistics
#
# =============================================================================

library(tidyverse)
library(GSVA)
library(limma)
library(RColorBrewer)

source("analysis/config.R")
source("analysis/R/theme_publication.R")

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------

cat("Loading data...\n")

expr_gene <- readRDS(file.path(RESULTS_DIR, "expr_gene_level.rds"))
metadata <- read_csv(file.path(RESULTS_DIR, "sample_metadata.csv"), show_col_types = FALSE)
ecm_axes_rat <- readRDS(file.path(RESULTS_DIR, "reference/ecm_axes_rat.rds"))

# Filter to timepoint samples only
metadata <- metadata %>% filter(timepoint != "Control")
expr_gene <- expr_gene[, metadata$sample_id]

cat(sprintf("  %d genes x %d samples\n", nrow(expr_gene), ncol(expr_gene)))

# -----------------------------------------------------------------------------
# 2. Run ssGSEA
# -----------------------------------------------------------------------------

cat("\nRunning ssGSEA...\n")

# Filter gene sets to genes present in data
expr_genes <- rownames(expr_gene)
gene_sets_filtered <- lapply(ecm_axes_rat, function(genes) {
  intersect(genes, expr_genes)
})

# Report coverage
for (axis in names(gene_sets_filtered)) {
  cat(sprintf("  %s: %d genes\n", axis, length(gene_sets_filtered[[axis]])))
}

# Run ssGSEA
gsva_param <- ssgseaParam(
  exprData = as.matrix(expr_gene),
  geneSets = gene_sets_filtered,
  minSize = 3,
  maxSize = 500
)

ssgsea_scores <- gsva(gsva_param)

cat(sprintf("\nssGSEA complete: %d axes x %d samples\n", 
            nrow(ssgsea_scores), ncol(ssgsea_scores)))

# -----------------------------------------------------------------------------
# 3. Statistical Testing
# -----------------------------------------------------------------------------

cat("\nTesting axis activation (Implant vs Control)...\n")

# Same design as DEG analysis
metadata$group <- factor(paste(metadata$condition, metadata$timepoint, sep = "_"))
design <- model.matrix(~ 0 + group, data = metadata)
colnames(design) <- levels(metadata$group)

metadata$pair_id <- paste(metadata$timepoint, metadata$replicate, sep = "_")

# Contrasts
timepoints <- c("Week0", "Week1", "Week2", "Week4", "Week18")
contrast_formulas <- sapply(timepoints, function(tp) {
  sprintf("Implant_%s - Control_%s", tp, tp)
})
contrasts <- makeContrasts(contrasts = contrast_formulas, levels = design)
colnames(contrasts) <- timepoints

# Fit model
corfit <- duplicateCorrelation(ssgsea_scores, design, block = metadata$pair_id)
fit <- lmFit(ssgsea_scores, design, block = metadata$pair_id, correlation = corfit$consensus)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# Extract results
axis_stats <- data.frame()
for (tp in timepoints) {
  results <- topTable(fit2, coef = tp, number = Inf)
  results$axis <- rownames(results)
  results$timepoint <- tp
  results$day <- TIMEPOINT_MAP[tp]
  axis_stats <- rbind(axis_stats, results)
}

# Save
axis_dir <- file.path(RESULTS_DIR, "axis_scoring")
dir.create(axis_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(axis_stats, file.path(axis_dir, "axis_activation_stats.csv"))
saveRDS(ssgsea_scores, file.path(axis_dir, "ssgsea_scores.rds"))

# Print summary
cat("\nAxis activation summary (logFC Implant vs Control):\n")
axis_summary <- axis_stats %>%
  select(axis, timepoint, logFC, adj.P.Val) %>%
  mutate(sig = ifelse(adj.P.Val < 0.05, "*", "")) %>%
  pivot_wider(names_from = timepoint, values_from = c(logFC, sig))
print(axis_summary)

# -----------------------------------------------------------------------------
# 4. Visualizations
# -----------------------------------------------------------------------------

cat("\nGenerating visualizations...\n")

# 4a. Prepare data for plotting
score_long <- as.data.frame(t(ssgsea_scores)) %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "axis", values_to = "score") %>%
  left_join(metadata, by = "sample_id")

score_long$day <- TIMEPOINT_MAP[score_long$timepoint]

# 4b. Trajectory plot (all axes)
p_trajectory <- ggplot(score_long, aes(x = day, y = score, color = condition, group = condition)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
  stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.8, linewidth = 0.8) +
  facet_wrap(~axis, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = c(0, 7, 14, 28, 126), trans = "pseudo_log") +
  scale_color_manual(values = colors_condition) +
  labs(
    title = "A",
    subtitle = "ECM Axis Enrichment Trajectories (Mean ± SE)",
    x = "Days Post-Implantation",
    y = "ssGSEA Score"
  ) +
  theme_publication(base_size = 10) +
  theme(legend.position = "bottom")

save_publication_figure(p_trajectory, file.path(FIGURES_DIR, "04_axis_trajectories"), width = 10, height = 8)

# 4c. Delta plot (Implant - Control difference)
score_wide <- score_long %>%
  select(axis, condition, timepoint, day, replicate, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  mutate(delta = Implant - Control)

p_delta <- ggplot(score_wide, aes(x = day, y = delta, color = axis, group = axis)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#757575", linewidth = 0.8) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
  stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.8, linewidth = 0.8) +
  scale_x_continuous(breaks = c(0, 7, 14, 28, 126), trans = "pseudo_log") +
  scale_color_manual(values = colors_axes, name = "ECM Axis") +
  labs(
    title = "B",
    subtitle = "Implant-Induced ECM Axis Perturbation (Δ Score)",
    x = "Days Post-Implantation",
    y = expression(Delta~ssGSEA~Score)
  ) +
  theme_publication(base_size = 11) +
  theme(legend.position = "right")

save_publication_figure(p_delta, file.path(FIGURES_DIR, "04_axis_delta"), width = 9, height = 6)

# 4d. Heatmap of axis scores
library(pheatmap)

# Calculate mean scores per group
mean_scores <- score_long %>%
  group_by(axis, condition, timepoint) %>%
  summarize(mean_score = mean(score), .groups = "drop") %>%
  unite("group", condition, timepoint, sep = "_") %>%
  pivot_wider(names_from = group, values_from = mean_score) %>%
  column_to_rownames("axis")

# Order columns by timepoint
col_order <- c(
  paste0("Control_", timepoints),
  paste0("Implant_", timepoints)
)
mean_scores <- mean_scores[, col_order]

# Annotation
annotation_col <- data.frame(
  Condition = rep(c("Control", "Implant"), each = 5),
  Timepoint = rep(timepoints, 2),
  row.names = col_order
)

ann_colors <- list(
  Condition = colors_condition,
  Timepoint = colors_timepoint
)

pdf(file.path(FIGURES_DIR, "04_axis_heatmap.pdf"), width = 10, height = 5)
pheatmap(
  mean_scores,
  scale = "row",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  border_color = "grey80",
  main = "ECM Axis Enrichment (Z-scaled)"
)
dev.off()

# 4e. Axis ranking at acute timepoint
acute_ranking <- axis_stats %>%
  filter(timepoint == "Week0") %>%
  arrange(desc(logFC)) %>%
  mutate(rank = row_number())

cat("\nAxis ranking at Week0 (acute):\n")
print(select(acute_ranking, rank, axis, logFC, adj.P.Val))

p_rank <- ggplot(acute_ranking, aes(x = reorder(axis, logFC), y = logFC, fill = logFC > 0)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "#1A9850"), guide = "none") +
  labs(
    title = "C",
    subtitle = "ECM Axis Activation at Week 0 (Acute Phase)",
    x = "",
    y = expression(log[2]~Fold~Change~(ssGSEA~Score))
  ) +
  theme_publication(base_size = 11)

save_publication_figure(p_rank, file.path(FIGURES_DIR, "04_axis_ranking_acute"), type = "barplot")

cat("\nScript 04 complete!\n")
cat(sprintf("  ssGSEA scores: %s\n", file.path(axis_dir, "ssgsea_scores.rds")))
cat(sprintf("  Activation stats: %s\n", file.path(axis_dir, "axis_activation_stats.csv")))
