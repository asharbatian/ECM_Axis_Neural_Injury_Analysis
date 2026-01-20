# scripts/04_axis_scoring.R

## Purpose

Calculates single-sample enrichment scores (ssGSEA) for each ECM axis to quantify pathway-level activation across samples and timepoints.

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
| `results/axis_scoring/ssgsea_scores.rds` | Raw ssGSEA score matrix (axes × samples) |
| `results/axis_scoring/axis_activation_stats.csv` | Statistical comparison per axis per timepoint |
| `results/axis_scoring/method_comparison.csv` | ssGSEA vs z-score correlation |
| `results/axis_scoring/zscore_scores.rds` | Alternative z-score method scores |

### Figures
| File | Description |
|------|-------------|
| `figures/04_axis_trajectories.pdf` | ECM axis scores over time by condition |
| `figures/04_axis_delta.pdf` | Delta scores (Implant - Control) over time |
| `figures/04_axis_ranking_acute.pdf` | Axis ranking at peak response |
| `figures/supplementary/04_method_comparison.pdf` | ssGSEA vs z-score comparison |

## Interpretation (for readers)

- **What this step answers**: “At the pathway/axis level, which ECM programs are most activated by implantation?”
- **Open first**:
  - `figures/04_axis_trajectories.pdf` (axis scores over time)
  - `results/axis_scoring/axis_activation_stats.csv` (statistics supporting differences)
- **How to interpret key columns**:
  - **`logFC`** > 0: higher axis score in **implant** vs control at that timepoint
  - **`adj.P.Val` (FDR)** < 0.05: statistically supported axis activation
- **Biological intuition**: axes summarize coordinated behavior of multiple genes, which can be more robust and interpretable than single-gene DEGs.

## Methods

### ssGSEA Algorithm

Single-sample Gene Set Enrichment Analysis implemented via GSVA package:

```r
library(GSVA)

gsva_params <- ssgseaParam(
  exprData = expr_gene,
  geneSets = ecm_axes,
  minSize = 3,
  maxSize = 500
)

ssgsea_scores <- gsva(gsva_params)
```

**Parameters**:
- `minSize = 3` — Minimum gene set size
- `maxSize = 500` — Maximum gene set size

### Statistical Testing

Enrichment scores compared between conditions using an unpaired limma design (implant vs control, independent samples):

```r
fit <- lmFit(ssgsea_scores, design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# Extract results per timepoint
axis_stats <- topTable(fit2, coef = "Week1", number = Inf)
```

### Sensitivity Analysis

Alternative method (mean z-score) for validation:
```r
# Z-score normalize genes
z_expr <- t(scale(t(expr_gene)))

# Calculate mean z-score per gene set
zscore_scores <- sapply(gene_sets, function(genes) {
  colMeans(z_expr[genes, ], na.rm = TRUE)
})
```

## ECM Axes Tested

| Axis | N Genes |
|------|---------|
| Hyaluronan | 12 |
| Provisional Matrix | 8 |
| PNN-CSPG | 10 |
| Basement Membrane | 8 |
| Proteases/Regulators | 12 |
| Crosslinking/Fibrosis | 13 |

## Dependencies

```r
library(GSVA)
library(limma)
library(tidyverse)
library(ggplot2)
```

## Output Columns

### axis_activation_stats.csv
| Column | Description |
|--------|-------------|
| axis | ECM axis name |
| timepoint | Week0, Week1, etc. |
| day | Numeric day (0, 7, 14, 28, 126) |
| logFC | Log2 fold change in ssGSEA score |
| t | Moderated t-statistic |
| P.Value | Raw p-value |
| adj.P.Val | BH-adjusted p-value |

