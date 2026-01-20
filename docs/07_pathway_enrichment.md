# scripts/07_pathway_enrichment.R

## Purpose

Performs fast gene set enrichment analysis (fGSEA) on ranked gene lists to test enrichment of ECM axes and literature-derived signatures.

## Input

| File | Source |
| --- | --- |
| `results/deg/deg_results_list.rds` | Script 03 |
| `results/reference/ecm_axes_rat.rds` | Script 02 |

## Output

### Data Files

| File | Description |
| --- | --- |
| `results/ha_analysis/gsea_all_results.csv` | All pathway enrichment results |
| `results/ha_analysis/go_enrichment_*.csv` | GO enrichment for peak timepoint |
| `results/ha_analysis/upstream_regulator_results.rds` | Full results object |

### Figures

| File | Description |
| --- | --- |
| `figures/individual/07_gsea_heatmap.pdf` | NES heatmap across timepoints |
| `figures/individual/07_ha_pathway_trajectory.pdf` | HA pathway NES over time |

## Interpretation (for readers)

- **What this step answers**: “Are whole pathways/gene sets (ECM axes, HA-related signaling, literature signatures) enriched after implantation?”
- **Open first**:
  - `results/ha_analysis/gsea_all_results.csv` (all enrichment results)
  - `figures/individual/07_gsea_heatmap.pdf` (quick overview across time)
- **How to interpret key columns**:
  - **NES** > 0: gene set is enriched among genes higher in **implant** (relative to control)
  - **NES** < 0: enriched among genes higher in **control**
  - **padj (FDR)** < 0.05: statistically supported enrichment after multiple testing correction
- **Biological value**: enrichment helps connect many modest gene changes into a coherent biological mechanism/signature.

## Methods

### Gene Ranking

Genes ranked by signed significance:

```r
rank_stat <- sign(logFC) * -log10(P.Value)
ranked_genes <- sort(rank_stat, decreasing = TRUE)
```

### fGSEA Parameters

```r
library(fgsea)

fgsea_results <- fgsea(
  pathways = gene_sets,
  stats = ranked_genes,
  minSize = 5,
  maxSize = 500,
  nPermSimple = 10000
)
```

### Gene Sets Tested

| Gene Set | Source | N Genes |
| --- | --- | --- |
| Hyaluronan | Manuscript Table 1 | 12 |
| Provisional_Matrix | Manuscript Table 1 | 8 |
| PNN_CSPG | Manuscript Table 1 | 10 |
| Basement_Membrane | Manuscript Table 1 | 8 |
| Proteases_Regulators | Manuscript Table 1 | 12 |
| Crosslinking_Fibrosis | Manuscript Table 1 | 13 |
| HA_DAMP_Signaling | CD44, TLR2/4, CD14, NF-κB targets | 10 |
| HA_Metabolism | HAS1-3, HYAL1-3 | 7 |
| Huff_Chronic | Huff et al. chronic markers | 10 |
| Joseph_Acute | Joseph et al. acute markers | 9 |

## Output Columns

### gsea_all_results.csv

| Column | Description |
| --- | --- |
| pathway | Gene set name |
| pval | Nominal p-value |
| padj | BH-adjusted p-value |
| log2err | Log2 error estimate |
| ES | Enrichment score |
| NES | Normalized enrichment score |
| size | Number of genes in set |
| leadingEdge | Core enriched genes |
| timepoint | Timepoint analyzed |

## Dependencies

```r
library(fgsea)
library(clusterProfiler)
library(org.Rn.eg.db)
library(tidyverse)
library(ggplot2)
```
