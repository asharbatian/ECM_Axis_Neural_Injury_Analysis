# =============================================================================
# HA Axis Validation Study — Script 01: Data Processing and Quality Control
# -----------------------------------------------------------------------------
# Purpose: Load CEL files, normalize with RMA, perform QC, and create metadata.
# Inputs:  CEL files in Data/arrays/Controls/ and Data/arrays/Implants/
# Outputs: Normalized expression matrix and sample metadata (see docs/01_data_processing.md)
#
# Author: Dr.-Ing Kevin Joseph
# Group Leader - Laboratory of NeuroEngineering
# Department of Neurosurgery
# Medical Center - University of Freiburg
# =============================================================================

library(tidyverse)
library(oligo)
library(limma)

source("scripts/config.R")
source("scripts/theme_publication.R")

# -----------------------------------------------------------------------------
# 1. Load CEL files
# -----------------------------------------------------------------------------

cat("Loading CEL files...\n")

# Get file paths
controls_dir <- file.path(DATA_DIR, "Controls")
implants_dir <- file.path(DATA_DIR, "Implants")

control_files <- list.files(controls_dir, pattern = "\\.CEL$", full.names = TRUE)
implant_files <- list.files(implants_dir, pattern = "\\.CEL$", full.names = TRUE)

cat(sprintf("  Found %d control files\n", length(control_files)))
cat(sprintf("  Found %d implant files\n", length(implant_files)))

# Rename files to include condition prefix before loading
# Create temp symlinks or just process each folder separately and combine after RMA

# Process Controls
cat("  Processing control samples...\n")
raw_controls <- read.celfiles(control_files)
eset_controls <- rma(raw_controls)
expr_controls <- exprs(eset_controls)
colnames(expr_controls) <- paste0("Control_", colnames(expr_controls))

# Process Implants
cat("  Processing implant samples...\n")
raw_implants <- read.celfiles(implant_files)
eset_implants <- rma(raw_implants)
expr_implants <- exprs(eset_implants)
colnames(expr_implants) <- paste0("Implant_", colnames(expr_implants))

# -----------------------------------------------------------------------------
# 2. Combine Expression Matrices
# -----------------------------------------------------------------------------

cat("Combining normalized expression matrices...\n")

# Combine the two matrices (same probes, different samples)
expr_matrix <- cbind(expr_controls, expr_implants)

cat(sprintf("  Expression matrix: %d probes x %d samples\n", 
            nrow(expr_matrix), ncol(expr_matrix)))

# -----------------------------------------------------------------------------
# 3. Create Sample Metadata
# -----------------------------------------------------------------------------

cat("Creating sample metadata...\n")

# Parse sample names (now have prefix: Control_Week1_1.CEL or Implant_Week1_1.CEL)
sample_names <- colnames(expr_matrix)

# Extract condition from prefix
conditions <- ifelse(grepl("^Control_", sample_names), "Control", "Implant")

# Extract timepoint and replicate from filename
parse_sample <- function(name) {
  # Remove condition prefix and extension
  # Format is: Control_Week1_1.CEL or Implant_Control_1.CEL
  base <- gsub("^(Control_|Implant_)", "", name)
  base <- gsub("\\.CEL$", "", base)
  
  # Parse timepoint and replicate
  if (grepl("^Control_", base)) {
    # This is a baseline control sample (e.g., "Control_1")
    timepoint <- "Baseline"
    replicate <- as.numeric(gsub("Control_", "", base))
  } else {
    # This is a timepoint sample (e.g., "Week1_3")
    parts <- strsplit(base, "_")[[1]]
    timepoint <- parts[1]
    replicate <- as.numeric(parts[2])
  }
  
  return(c(timepoint = timepoint, replicate = replicate))
}

parsed <- t(sapply(sample_names, parse_sample))

metadata <- data.frame(
  sample_id = sample_names,
  condition = conditions,
  timepoint = parsed[, "timepoint"],
  replicate = as.numeric(parsed[, "replicate"]),
  stringsAsFactors = FALSE
)

# Add day mapping
metadata$day <- TIMEPOINT_MAP[metadata$timepoint]
metadata$day[metadata$timepoint == "Control"] <- NA

## NOTE (unpaired analysis):
## Animal identity is not tracked across implant/control due to randomization.
## The `replicate` number is treated only as a sample index within a timepoint,
## not a biological pairing variable.

# Clean up column names
colnames(expr_matrix) <- metadata$sample_id

cat("Sample counts by condition and timepoint:\n")
print(table(metadata$condition, metadata$timepoint))

# -----------------------------------------------------------------------------
# 4. Quality Control
# -----------------------------------------------------------------------------

cat("\nPerforming quality control...\n")

# 4a. Sample correlations
cor_matrix <- cor(expr_matrix, method = "spearman")
mean_cors <- rowMeans(cor_matrix)

# Flag low-correlation samples
low_cor_threshold <- 0.85
low_cor_samples <- names(mean_cors[mean_cors < low_cor_threshold])

if (length(low_cor_samples) > 0) {
  cat(sprintf("  WARNING: %d samples with mean correlation < %.2f:\n", 
              length(low_cor_samples), low_cor_threshold))
  print(low_cor_samples)
} else {
  cat("  All samples pass correlation QC\n")
}

# 4b. PCA
pca <- prcomp(t(expr_matrix), scale = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],
  metadata
)

var_explained <- summary(pca)$importance[2, 1:3] * 100

# 4c. PCA plot
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = timepoint, shape = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = colors_timepoint) +
  labs(
    title = "A",
    subtitle = "PCA of Normalized Expression",
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_publication(base_size = 11) +
  theme(legend.position = "right")

save_publication_figure(p_pca, file.path(FIGURES_DIR, "01_pca_all_samples"), width = 8, height = 6)

# 4d. Sample correlation heatmap
library(pheatmap)
library(RColorBrewer)

annotation_col <- as.data.frame(metadata)
rownames(annotation_col) <- annotation_col$sample_id
annotation_col <- annotation_col[, c("condition", "timepoint")]

ann_colors <- list(
  condition = colors_condition,
  timepoint = c(colors_timepoint, "Baseline" = "#FFFFFF")
)

pdf(file.path(FIGURES_DIR, "01_correlation_heatmap.pdf"), width = 12, height = 10)
pheatmap(
  cor_matrix,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  border_color = "grey80",
  main = "Sample Correlation Matrix (Spearman)"
)
dev.off()

# -----------------------------------------------------------------------------
# 5. Save Outputs
# -----------------------------------------------------------------------------

cat("\nSaving outputs...\n")

# Create results subdirectory
qc_dir <- file.path(RESULTS_DIR, "qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# Save expression matrix
saveRDS(expr_matrix, file.path(RESULTS_DIR, "expr_matrix_normalized.rds"))

# Save metadata
write_csv(metadata, file.path(RESULTS_DIR, "sample_metadata.csv"))

# Save QC metrics
qc_metrics <- data.frame(
  sample_id = names(mean_cors),
  mean_correlation = mean_cors,
  pass_qc = mean_cors >= low_cor_threshold
)
write_csv(qc_metrics, file.path(qc_dir, "sample_correlations.csv"))

# Save PCA results
write_csv(pca_df, file.path(qc_dir, "pca_coordinates.csv"))

cat("\nScript 01 complete!\n")
cat(sprintf("  Expression matrix: %s\n", file.path(RESULTS_DIR, "expr_matrix_normalized.rds")))
cat(sprintf("  Sample metadata: %s\n", file.path(RESULTS_DIR, "sample_metadata.csv")))
