# scripts/05_temporal_dynamics.R

## Purpose

Classifies differentially expressed genes into temporal patterns based on when they are significant across the time course.

## Input

| File | Source |
|------|--------|
| `results/deg/deg_results_list.rds` | Script 03 |
| `results/sample_metadata.csv` | Script 01 |
| `results/reference/ecm_axes_rat.rds` | Script 02 |

## Output

### Data Files
| File | Description |
|------|-------------|
| `results/temporal/gene_temporal_classification.csv` | Gene-level pattern assignments |
| `results/temporal/axis_temporal_patterns.csv` | Axis-level pattern counts |
| `results/temporal/axis_pattern_matrix.csv` | Matrix of axis × pattern counts |

### Figures
| File | Description |
|------|-------------|
| `figures/05_gene_temporal_patterns.pdf` | Bar plot of pattern distribution |
| `figures/05_axis_temporal_heatmap.pdf` | Heatmap of axis × pattern |

## Interpretation (for readers)

- **What this step answers**: “Do implantation-induced changes resolve or persist into the chronic timepoint?”
- **Open first**:
  - `results/temporal/gene_temporal_classification.csv` (gene-level patterns)
  - `figures/05_gene_temporal_patterns.pdf` (overall pattern distribution)
- **How to read the patterns**:
  - **Resolving**: early changes that are not present late (often interpreted as recovery/adaptation)
  - **Persistent**: present early and still present late (often interpreted as chronic remodeling)
  - **Late-emerging**: minimal early but appears late (often interpreted as delayed/chronic processes)
- **Paper-friendly takeaway**: this is the simplest place to quantify “acute vs chronic” biology.

## Methods

### Classification Criteria

Genes are classified based on whether they appear as significant DEGs in early, mid, and/or late phases.

| Pattern | Definition |
|---------|------------|
| Persistent | Significant in early, mid, and late phases |
| Resolving | Significant in early phase, not significant in late phase |
| Late-emerging | Not significant in early phase, significant in late phase |
| Biphasic | Significant in early and late phases, not in mid phase |
| Transient-mid | Significant only in mid phase |
| Early-to-mid | Significant in early and mid phases, not late |
| Mid-to-late | Significant in mid and late phases, not early |
| Other | Any remaining pattern |

### Code

```r
# Define phases (example)
early_phase <- c("Week0", "Week1")
mid_phase   <- c("Week2", "Week4")
late_phase  <- c("Week18")

# Collect significant genes by phase (deg_list contains per-timepoint DEG tables)
early_genes <- unique(unlist(lapply(early_phase, function(tp) deg_list[[tp]]$gene)))
mid_genes   <- unique(unlist(lapply(mid_phase,   function(tp) deg_list[[tp]]$gene)))
late_genes  <- unique(unlist(lapply(late_phase,  function(tp) deg_list[[tp]]$gene)))

classify_gene <- function(gene) {
  in_early <- gene %in% early_genes
  in_mid  <- gene %in% mid_genes
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
```

## Output Columns

### gene_temporal_classification.csv
| Column | Description |
|--------|-------------|
| gene | Gene symbol |
| n_sig_timepoints | Number of significant timepoints |
| early_significant | TRUE if DEG at Day 0-14 |
| late_significant | TRUE if DEG at Day 126 |
| pattern | Classification result |
| max_logFC | Maximum absolute logFC |
| peak_timepoint | Timepoint with max |logFC| |

## Dependencies

```r
library(tidyverse)
library(ggplot2)
library(pheatmap)
```
