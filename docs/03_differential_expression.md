# 03_differential_expression.R

## Purpose

Performs differential expression analysis between implant and control samples at each timepoint using the `limma` framework.

## Input

### Required Files
| File | Source |
|------|--------|
| `results/expr_gene_level.rds` | Script 02 |
| `results/sample_metadata.csv` | Script 01 |
| `results/reference/ecm_axes_rat.rds` | Script 02 |

## Output

### Data Files
| File | Description |
|------|-------------|
| `results/deg/deg_results_list.rds` | List of full DEG results per timepoint |
| `results/deg/limma_fit.rds` | Fitted limma model object |
| `results/deg/deg_all_Week*.csv` | All genes with statistics per timepoint |
| `results/deg/deg_significant_Week*.csv` | Significant DEGs only per timepoint |
| `results/deg/deg_summary.csv` | DEG counts per timepoint |
| `results/deg/axis_enrichment.csv` | ECM axis enrichment in DEGs |

### Figures
| File | Description |
|------|-------------|
| `figures/individual/03_volcano_*.pdf` | Volcano plots per timepoint |
| `figures/supplementary/03_volcano_all_timepoints.pdf` | Combined volcano plot |
| `figures/individual/03_deg_trajectory.pdf` | DEG count over time |

## Methods

### Statistical Model

```r
# Design matrix
design <- model.matrix(~ 0 + group, data = metadata)
# where group = condition:timepoint interaction

# Contrasts: Implant - Control at each timepoint
contrasts <- makeContrasts(
  Week0 = Implant.Week0 - Control.Week0,
  Week1 = Implant.Week1 - Control.Week1,
  ...
)

# Account for paired design
corfit <- duplicateCorrelation(expr, design, block = metadata$animal_id)
fit <- lmFit(expr, design, block = metadata$animal_id, 
             correlation = corfit$consensus)
```

### Significance Criteria

| Criterion | Threshold |
|-----------|-----------|
| Adjusted p-value | < 0.05 (Benjamini-Hochberg FDR) |
| Log2 fold change | ≥ 1.0 (absolute value) |

### ECM Axis Enrichment

For each timepoint, enrichment of ECM axis genes in DEGs is calculated:
- **Overlap**: Number of axis genes that are DEGs
- **Fold enrichment**: Observed / Expected overlap
- **P-value**: Fisher's exact test

## Dependencies

```r
library(limma)
library(tidyverse)
library(ggplot2)
```

## Key Output Columns

### deg_all_Week*.csv
| Column | Description |
|--------|-------------|
| gene | Gene symbol |
| logFC | Log2 fold change (Implant - Control) |
| AveExpr | Average expression across all samples |
| t | Moderated t-statistic |
| P.Value | Raw p-value |
| adj.P.Val | BH-adjusted p-value |
| B | Log-odds of differential expression |

### axis_enrichment.csv
| Column | Description |
|--------|-------------|
| timepoint | Week0, Week1, etc. |
| axis | ECM axis name |
| overlap | N axis genes that are DEGs |
| overlap_genes | Gene symbols |
| fold_enrichment | Enrichment ratio |
| p_value | Fisher's exact test p-value |

