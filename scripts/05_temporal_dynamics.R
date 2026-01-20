# =============================================================================
# HA Axis Validation Study — Script 05: Temporal Dynamics and Chronic Analysis
# -----------------------------------------------------------------------------
# Purpose: Analyze temporal patterns of ECM axis engagement and classify persistent
#          vs resolving vs late-emerging signatures across time.
# Inputs:  DEG results and axis scores from prior scripts
# Outputs: Temporal classification tables and summary metrics (see docs/05_temporal_dynamics.md)
#
# Author: Dr.-Ing Kevin Joseph
# Group Leader - Laboratory of NeuroEngineering
# Department of Neurosurgery
# Medical Center - University of Freiburg
# =============================================================================

library(tidyverse)
library(RColorBrewer)

source("scripts/config.R")
source("scripts/theme_publication.R")

# Utility colors for plots
colors_utility <- c(
  significant = "#1B9E77",
  nonsignificant = "grey80"
)

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------

cat("Loading data...\n")

deg_dir <- file.path(RESULTS_DIR, "deg")
axis_dir <- file.path(RESULTS_DIR, "axis_scoring")

# Load DEGs from all timepoints
timepoints <- c("Week0", "Week1", "Week2", "Week4", "Week18")

deg_list <- lapply(timepoints, function(tp) {
  read_csv(file.path(deg_dir, sprintf("deg_significant_%s.csv", tp)), show_col_types = FALSE)
})
names(deg_list) <- timepoints

# Load axis scores and stats
ssgsea_scores <- readRDS(file.path(axis_dir, "ssgsea_scores.rds"))
axis_stats <- read_csv(file.path(axis_dir, "axis_activation_stats.csv"), show_col_types = FALSE)
ecm_axes_rat <- readRDS(file.path(RESULTS_DIR, "reference/ecm_axes_rat.rds"))

cat("DEG counts by timepoint:\n")
print(sapply(deg_list, nrow))

# -----------------------------------------------------------------------------
# 2. Classify Genes by Temporal Pattern
# -----------------------------------------------------------------------------

cat("\nClassifying genes by temporal pattern...\n")

# Define phases
early_phase <- c("Week0", "Week1")      # Days 0-7
mid_phase <- c("Week2", "Week4")        # Days 14-28
late_phase <- c("Week18")               # Day 126

# Get unique DEGs in each phase
early_genes <- unique(unlist(lapply(early_phase, function(tp) deg_list[[tp]]$gene)))
mid_genes <- unique(unlist(lapply(mid_phase, function(tp) deg_list[[tp]]$gene)))
late_genes <- unique(unlist(lapply(late_phase, function(tp) deg_list[[tp]]$gene)))

all_degs <- unique(c(early_genes, mid_genes, late_genes))

# Classify each gene
classify_gene <- function(gene) {
  in_early <- gene %in% early_genes
  in_mid <- gene %in% mid_genes
  in_late <- gene %in% late_genes
  
  if (in_early && in_mid && in_late) return("Persistent")
  if (in_early && !in_late) return("Resolving")
  if (!in_early && in_late) return("Late-emerging")
  if (in_early && in_late && !in_mid) return("Biphasic")
  if (!in_early && in_mid && !in_late) return("Transient-mid")
  if (in_early && in_mid && !in_late) return("Early-to-mid")
  if (!in_early && in_mid && in_late) return("Mid-to-late")
  return("Other")
}

gene_classification <- data.frame(
  gene = all_degs,
  pattern = sapply(all_degs, classify_gene)
)

cat("\nGene temporal patterns:\n")
print(table(gene_classification$pattern))

# Save classification
temporal_dir <- file.path(RESULTS_DIR, "temporal")
dir.create(temporal_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(gene_classification, file.path(temporal_dir, "gene_temporal_classification.csv"))

# -----------------------------------------------------------------------------
# 3. Axis Engagement by Temporal Pattern
# -----------------------------------------------------------------------------

cat("\nAnalyzing axis engagement by temporal pattern...\n")

axis_temporal <- data.frame()

for (axis in names(ecm_axes_rat)) {
  axis_genes <- ecm_axes_rat[[axis]]
  
  for (pattern in unique(gene_classification$pattern)) {
    pattern_genes <- gene_classification$gene[gene_classification$pattern == pattern]
    overlap <- intersect(axis_genes, pattern_genes)
    
    row <- data.frame(
      axis = axis,
      pattern = pattern,
      gene_count = length(overlap),
      genes = paste(overlap, collapse = ";")
    )
    axis_temporal <- rbind(axis_temporal, row)
  }
}

# Pivot for easier viewing
axis_pattern_matrix <- axis_temporal %>%
  select(axis, pattern, gene_count) %>%
  pivot_wider(names_from = pattern, values_from = gene_count, values_fill = 0)

cat("\nAxis genes by temporal pattern:\n")
print(axis_pattern_matrix)

write_csv(axis_temporal, file.path(temporal_dir, "axis_temporal_patterns.csv"))
write_csv(axis_pattern_matrix, file.path(temporal_dir, "axis_pattern_matrix.csv"))

# -----------------------------------------------------------------------------
# 4. Key Validation Questions
# -----------------------------------------------------------------------------

cat("\n" , rep("=", 60), "\n", sep = "")
cat("KEY VALIDATION QUESTIONS\n")
cat(rep("=", 60), "\n")

# Q6: Does HA axis engagement resolve by Week 18?
ha_early <- intersect(ecm_axes_rat$Hyaluronan, early_genes)
ha_late <- intersect(ecm_axes_rat$Hyaluronan, late_genes)
ha_resolving <- setdiff(ha_early, late_genes)
ha_persistent <- intersect(ha_early, late_genes)

cat("\nQ6: Does HA-axis engagement RESOLVE by Week 18?\n")
cat(sprintf("  HA genes DE at early phase: %d\n", length(ha_early)))
cat(sprintf("  HA genes DE at late phase: %d\n", length(ha_late)))
cat(sprintf("  HA genes that RESOLVE: %d (%s)\n", length(ha_resolving), paste(ha_resolving, collapse = ", ")))
cat(sprintf("  HA genes that PERSIST: %d (%s)\n", length(ha_persistent), paste(ha_persistent, collapse = ", ")))

if (length(ha_persistent) > length(ha_resolving)) {
  cat("  --> WARNING: HA axis engagement PERSISTS - potential chronic inflammation\n")
} else {
  cat("  --> HA axis shows resolution pattern - favorable biocompatibility\n")
}

# Q7: Does Crosslinking/Fibrosis show late emergence?
fibrosis_axis <- ecm_axes_rat$Crosslinking_Fibrosis
fibrosis_early <- intersect(fibrosis_axis, early_genes)
fibrosis_late <- intersect(fibrosis_axis, late_genes)
fibrosis_late_only <- setdiff(fibrosis_late, early_genes)

cat("\nQ7: Does Crosslinking/Fibrosis axis show late emergence?\n")
cat(sprintf("  Fibrosis genes DE at early phase: %d\n", length(fibrosis_early)))
cat(sprintf("  Fibrosis genes DE at late phase: %d\n", length(fibrosis_late)))
cat(sprintf("  Fibrosis genes LATE-EMERGING only: %d (%s)\n", 
            length(fibrosis_late_only), paste(fibrosis_late_only, collapse = ", ")))

# Q8: HA:Fibrosis ratio at Week 18
ha_score_w18 <- axis_stats %>% filter(axis == "Hyaluronan", timepoint == "Week18")
fibrosis_score_w18 <- axis_stats %>% filter(axis == "Crosslinking_Fibrosis", timepoint == "Week18")

cat("\nQ8: HA:Fibrosis ratio at Week 18?\n")
cat(sprintf("  HA axis logFC at Week18: %.3f (p=%.4f)\n", ha_score_w18$logFC, ha_score_w18$adj.P.Val))
cat(sprintf("  Fibrosis axis logFC at Week18: %.3f (p=%.4f)\n", fibrosis_score_w18$logFC, fibrosis_score_w18$adj.P.Val))

if (ha_score_w18$logFC != 0 && fibrosis_score_w18$logFC != 0) {
  ratio <- ha_score_w18$logFC / fibrosis_score_w18$logFC
  cat(sprintf("  HA:Fibrosis ratio: %.2f\n", ratio))
}

# Q11: Significant implant-vs-control difference at Week 18?
cat("\nQ11: Any axis significantly different at Week 18?\n")
w18_stats <- axis_stats %>% 
  filter(timepoint == "Week18") %>%
  arrange(adj.P.Val)
print(select(w18_stats, axis, logFC, P.Value, adj.P.Val))

sig_w18 <- w18_stats %>% filter(adj.P.Val < 0.05)
if (nrow(sig_w18) > 0) {
  cat(sprintf("  --> %d axes still significantly perturbed at Week 18\n", nrow(sig_w18)))
} else {
  cat("  --> No axes significantly perturbed - tissue returning to homeostasis\n")
}

# -----------------------------------------------------------------------------
# 5. Visualizations
# -----------------------------------------------------------------------------

cat("\nGenerating visualizations...\n")

# 5a. Gene classification bar chart (more publication-appropriate than pie)
pattern_counts <- as.data.frame(table(gene_classification$pattern))
colnames(pattern_counts) <- c("Pattern", "Count")
pattern_counts$Percent <- round(pattern_counts$Count / sum(pattern_counts$Count) * 100, 1)

p_patterns <- ggplot(pattern_counts, aes(x = reorder(Pattern, -Count), y = Count)) +
  geom_col(fill = colors_utility["significant"], width = 0.6, color = "black", linewidth = 0.4) +
  geom_text(aes(label = Count), vjust = -0.5, size = 3.5) +
  labs(
    title = "A",
    subtitle = "DEG Temporal Patterns",
    x = "",
    y = "Number of Genes"
  ) +
  theme_publication(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_publication_figure(p_patterns, file.path(FIGURES_DIR, "05_gene_temporal_patterns"), type = "barplot")

# 5b. Axis temporal pattern heatmap
pattern_order <- c("Persistent", "Resolving", "Late-emerging", "Biphasic", 
                   "Early-to-mid", "Transient-mid", "Mid-to-late", "Other")
axis_pattern_long <- axis_temporal %>%
  filter(pattern %in% pattern_order) %>%
  mutate(pattern = factor(pattern, levels = pattern_order))

p_axis_pattern <- ggplot(axis_pattern_long, aes(x = pattern, y = axis, fill = gene_count)) +
  geom_tile(color = "grey80", linewidth = 0.5) +
  geom_text(aes(label = gene_count), size = 3) +
  scale_fill_gradient(low = "white", high = "#B2182B", name = "Gene\nCount") +
  labs(
    title = "B",
    subtitle = "ECM Axis Gene Distribution by Temporal Pattern",
    x = "Temporal Pattern",
    y = ""
  ) +
  theme_publication(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_publication_figure(p_axis_pattern, file.path(FIGURES_DIR, "05_axis_temporal_heatmap"), type = "heatmap")

# 5c. DEG count trajectory
deg_counts <- data.frame(
  timepoint = factor(timepoints, levels = timepoints),
  day = TIMEPOINT_MAP[timepoints],
  count = sapply(deg_list, nrow)
)

p_deg_trajectory <- ggplot(deg_counts, aes(x = day, y = count)) +
  geom_area(fill = "#DEEBF7", alpha = 0.8) +
  geom_line(linewidth = 0.8, color = "#3182BD") +
  geom_point(size = 3, color = "#3182BD") +
  geom_text(aes(label = count), vjust = -1, size = 3.5) +
  scale_x_continuous(breaks = c(0, 7, 14, 28, 126), trans = "pseudo_log") +
  labs(
    title = "C",
    subtitle = expression(DEG~Count~Over~Time~("|"*log[2]*FC*"|">=1*","~adj.p<0.05*")")),
    x = "Days Post-Implantation",
    y = "Number of DEGs"
  ) +
  theme_publication(base_size = 11) +
  expand_limits(y = 0)

save_publication_figure(p_deg_trajectory, file.path(FIGURES_DIR, "05_deg_trajectory"), width = 7, height = 5)

# -----------------------------------------------------------------------------
# 5b. Calculate Biocompatibility Score
# -----------------------------------------------------------------------------

cat("\nCalculating Biocompatibility Score...\n")

# Biocompatibility Score = HA_Resolution / (Fibrosis_Persistence + 0.1)
# Higher score = better biocompatibility

# Get axis stats at early (peak) vs late
early_tp <- "Week1"  # Peak response
late_tp <- "Week18"  # Chronic

ha_early <- axis_stats %>% filter(axis == "Hyaluronan", timepoint == early_tp)
ha_late <- axis_stats %>% filter(axis == "Hyaluronan", timepoint == late_tp)
fib_early <- axis_stats %>% filter(axis == "Crosslinking_Fibrosis", timepoint == early_tp)
fib_late <- axis_stats %>% filter(axis == "Crosslinking_Fibrosis", timepoint == late_tp)

if (nrow(ha_early) > 0 && nrow(ha_late) > 0) {
  ha_resolution <- ha_early$logFC - ha_late$logFC
  
  fib_persistence <- ifelse(
    nrow(fib_early) > 0 && fib_early$logFC != 0,
    fib_late$logFC / fib_early$logFC,
    ifelse(nrow(fib_late) > 0, fib_late$logFC, 0)
  )
  
  biocompatibility_score <- ha_resolution / (abs(fib_persistence) + 0.1)
  
  cat(sprintf("\n  BIOCOMPATIBILITY SCORE: %.2f\n", biocompatibility_score))
  cat(sprintf("    HA Resolution (early-late): %.2f\n", ha_resolution))
  cat(sprintf("    Fibrosis Persistence (late/early): %.2f\n", fib_persistence))
  
  if (biocompatibility_score > 1) {
    cat("    --> FAVORABLE: HA resolves faster than fibrosis accumulates\n")
  } else if (biocompatibility_score > 0) {
    cat("    --> MODERATE: Some resolution occurring\n")
  } else {
    cat("    --> UNFAVORABLE: Persistent inflammation or increasing fibrosis\n")
  }
} else {
  biocompatibility_score <- NA
  ha_resolution <- NA
  fib_persistence <- NA
  cat("  Could not calculate biocompatibility score (missing axis data)\n")
}

# -----------------------------------------------------------------------------
# 6. Save Summary
# -----------------------------------------------------------------------------

# Create validation summary
validation_summary <- list(
  biocompatibility_score = biocompatibility_score,
  ha_resolution = ha_resolution,
  fibrosis_persistence = fib_persistence,
  ha_axis = list(
    early_genes = ha_early,
    late_genes = ha_late,
    resolving = ha_resolving,
    persistent = ha_persistent,
    resolution_ratio = length(ha_resolving) / max(length(ha_early), 1)
  ),
  fibrosis_axis = list(
    early_genes = fibrosis_early,
    late_genes = fibrosis_late,
    late_emerging = fibrosis_late_only
  ),
  week18_stats = w18_stats,
  gene_patterns = table(gene_classification$pattern)
)

saveRDS(validation_summary, file.path(temporal_dir, "validation_summary.rds"))

cat("\nScript 05 complete!\n")
cat(sprintf("  Temporal results: %s\n", temporal_dir))
