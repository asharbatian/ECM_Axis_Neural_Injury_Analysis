# =============================================================================
# HA Axis Validation Study — Project Configuration
# -----------------------------------------------------------------------------
# This file centralizes file paths and global parameters used across scripts.
#
# All scripts are intended to be run from the repository root so relative paths
# resolve correctly (e.g., `Data/`, `results/`, `figures/`).
# =============================================================================

# ----------------------------
# Directory layout (relative)
# ----------------------------
DATA_DIR <- file.path("Data", "arrays")
RESULTS_DIR <- "results"
FIGURES_DIR <- "figures"

# Ensure output directories exist
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Study parameters
# ----------------------------
# Timepoint naming follows the CEL file naming convention used in this project.
TIMEPOINTS <- c("Week0", "Week1", "Week2", "Week4", "Week18")

# Map timepoint label -> day post-implant (for plotting / summaries)
# Week0 corresponds to Day 0 (4h).
TIMEPOINT_MAP <- c(
  Baseline = NA_real_,
  Week0 = 0,
  Week1 = 7,
  Week2 = 14,
  Week4 = 28,
  Week18 = 126
)

# Differential expression thresholds used across scripts
DE_PADJ_THRESHOLD <- 0.05
DE_LFC_THRESHOLD <- 1.0

# Plotting palettes (used consistently across scripts)
colors_condition <- c(
  Control = "#2166AC",
  Implant = "#B2182B"
)

colors_timepoint <- c(
  Baseline = "#999999",
  Week0 = "#1B9E77",
  Week1 = "#D95F02",
  Week2 = "#7570B3",
  Week4 = "#E7298A",
  Week18 = "#66A61E"
)
