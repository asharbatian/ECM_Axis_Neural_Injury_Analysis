# scripts/02_annotation_mapping.R

## Purpose

Maps Affymetrix probe IDs to gene symbols and creates gene-level expression matrices. Also translates mouse ECM axis gene sets to rat orthologs.

## Input

### Required Files
| File | Source |
|------|--------|
| `results/expr_matrix_normalized.rds` | Script 01 |

### Dependencies
- `clariomsrattranscriptcluster.db` — Clariom S Rat annotation package
- `org.Rn.eg.db` — Rat gene annotation

## Output

### Data Files
| File | Description |
|------|-------------|
| `results/expr_gene_level.rds` | Gene-level expression matrix (genes × samples) |
| `results/reference/ecm_axes_mouse.rds` | Original mouse ECM axis gene sets |
| `results/reference/ecm_axes_rat.rds` | Rat ortholog ECM axis gene sets |
| `results/reference/mouse_rat_orthologs.csv` | Mouse-to-rat gene symbol mapping |
| `results/reference/ortholog_mapping_summary.csv` | Mapping success rates per axis |
| `results/reference/probe_annotations.csv` | Full probe-to-gene annotation table |

## Interpretation (for readers)

- **What this step answers**: “Are we measuring genes (not probes), and do we have the ECM axis gene lists in the same species as the data?”
- **Open first**: `results/reference/ortholog_mapping_summary.csv`
- **How to read it**:
  - High mapping success means the axis gene lists are well represented on the rat array.
  - If an axis has many unmapped genes, axis-level scores for that axis should be interpreted more cautiously.
- **Why this matters biologically**: all downstream biology (DEGs, axis scoring, enrichment) relies on accurate **gene-level** measurements and correctly defined **axis gene sets**.

## Methods

### Probe-to-Gene Mapping

1. **Annotation Source**: `clariomsrattranscriptcluster.db`
2. **ID Type**: Transcript cluster IDs → Gene symbols
3. **Multi-mapping Resolution**: Select probe with highest mean expression per gene
4. **Unmapped Probes**: Excluded from gene-level matrix

### ECM Axis Gene Sets

Six ECM axes defined from manuscript Table 1:

| Axis | Mouse Genes | Rat Mapping |
|------|-------------|-------------|
| Hyaluronan | Has1, Has2, Has3, Hyal1-3, Cd44, Hmmr, Tlr2, Tlr4, Cd14 | Direct symbol match |
| Provisional Matrix | Fn1, Tnc, Spp1, Itga5, Itgb1, Itgav, Itgb3, Thbs1 | Direct symbol match |
| PNN-CSPG | Acan, Bcan, Vcan, Tnc, Tnr, Ptprz1, Ncan, Hapln1 | Direct symbol match |
| Basement Membrane | Hspg2, Lama4, Lamb1, Lamc1, Col4a1, Col4a2 | Direct symbol match |
| Proteases/Regulators | Adamts1/4/5, Mmp2/9/12/14, Timp1-3 | Direct symbol match |
| Crosslinking/Fibrosis | Lox, Loxl1/2, Col1a1/2, Col3a1, Tgm2 | Direct symbol match |

### Ortholog Strategy

For rat arrays, gene symbols are largely conserved:
1. **Primary**: Direct symbol match (case-insensitive)
2. **Fallback**: biomaRt ortholog query (if direct match fails)
3. **Validation**: Check presence in expression matrix

## Key Code

```r
# Probe annotation
probe_annotation <- AnnotationDbi::select(
  clariomsrattranscriptcluster.db,
  keys = rownames(expr_matrix),
  columns = c("SYMBOL", "GENENAME", "ENTREZID"),
  keytype = "PROBEID"
)

# Collapse to gene level (max expression probe per gene)
expr_gene <- expr_matrix %>%
  group_by(gene_symbol) %>%
  summarize(across(everything(), ~ .[which.max(rowMeans(.))])) 
```

## Dependencies

```r
library(clariomsrattranscriptcluster.db)
library(org.Rn.eg.db)
library(AnnotationDbi)
library(tidyverse)
```

