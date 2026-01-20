# =============================================================================
# HA Axis Validation Study — Publication Plot Theme Helpers
# -----------------------------------------------------------------------------
# Minimal ggplot2 theme + helper to save figures consistently.
#
# These helpers are intentionally lightweight so the analysis can run in a fresh
# environment without extra bespoke dependencies.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

theme_publication <- function(base_size = 11, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

# Save as PDF + PNG (300 dpi) using a filename *stem* without extension.
# Accepts and ignores extra arguments (e.g., legacy `type=`) for compatibility.
save_publication_figure <- function(plot, filepath_stem, width = 8, height = 6, dpi = 300, ...) {
  dir.create(dirname(filepath_stem), recursive = TRUE, showWarnings = FALSE)

  ggsave(paste0(filepath_stem, ".pdf"), plot = plot, width = width, height = height, units = "in")
  ggsave(paste0(filepath_stem, ".png"), plot = plot, width = width, height = height, units = "in", dpi = dpi)
}

