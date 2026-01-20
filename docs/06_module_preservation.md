# scripts/06_module_preservation.R

## Purpose

Tests whether WGCNA co-expression modules from brain injury (BI) and spinal cord injury (SCI) datasets show differential activity in the neural implant dataset.

## Input

| File | Source |
|------|--------|
| `results/expr_gene_level.rds` | Script 02 |
| `results/sample_metadata.csv` | Script 01 |

## Output

### Data Files
| File | Description |
|------|-------------|
| `results/preservation/module_eigengenes.rds` | Module eigengene values per sample |
| `results/preservation/module_activity_stats.csv` | Statistical comparison per module per timepoint |
| `results/preservation/module_preservation_summary.csv` | Summary of preservation per module |

### Figures
| File | Description |
|------|-------------|
| `figures/06_key_module_trajectories.pdf` | Module eigengene trajectories |
| `figures/06_module_preservation_bar.pdf` | Bar plot of significant timepoints per module |

## Interpretation (for readers)

- **What this step answers**: “Do known injury-related gene programs (modules) also show activity changes after implantation?”
- **Open first**:
  - `results/preservation/module_preservation_summary.csv` (which modules are consistently affected)
  - `figures/06_key_module_trajectories.pdf` (how module activity evolves over time)
- **How to interpret it**:
  - A “module” is a **gene program**; its “eigengene” is a single summary value representing overall module activity per sample.
  - If inflammatory/ECM modules change in implant vs control, it supports shared biology with BI/SCI responses.
- **Important note**: these modules are manuscript-derived summaries; interpret them as “program-level concordance,” not a definitive causal mechanism.

## Methods

### Module Eigengene Calculation

First principal component of module genes:

```r
calculate_eigengene <- function(expr_matrix, module_genes) {
  module_expr <- expr_matrix[module_genes, ]
  pca <- prcomp(t(module_expr), scale. = TRUE)
  eigengene <- pca$x[, 1]
  
  # Ensure positive correlation with mean expression
  if (cor(eigengene, colMeans(module_expr)) < 0) {
    eigengene <- -eigengene
  }
  return(eigengene)
}
```

### Statistical Testing

Module eigengenes are compared between conditions using an **unpaired** test within each timepoint (implant vs control; independent samples).

### Preservation Criterion

Module is "preserved" if significant (p < 0.05) at ≥2 timepoints.

## Reference Modules

| Module | Source | N Genes |
|--------|--------|---------|
| BI_M1 - BI_M6 | Brain injury dataset | 50-200 each |
| SCI_M1 - SCI_M6 | Spinal cord injury dataset | 50-200 each |

## Output Columns

### module_preservation_summary.csv
| Column | Description |
|--------|-------------|
| module | Module ID (e.g., BI_M2) |
| n_genes | Number of genes in module |
| n_timepoints_sig | Timepoints with p < 0.05 |
| mean_delta | Mean eigengene difference |
| preserved | TRUE if n_timepoints_sig ≥ 2 |

## Dependencies

```r
library(tidyverse)
library(limma)
library(ggplot2)
```
