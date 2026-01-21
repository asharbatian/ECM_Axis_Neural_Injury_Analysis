# ECM Axis in Neural Injury and Implant Response

Analysis scripts for the study: **"ECM-DAMP Signaling in Neural Injury and Implant Tissue Response"**

This repository contains both **discovery** (BI/SCI) and **validation** (neural implant) analysis pipelines.

---

## Contributors

| Analysis | Author | Folder |
|----------|--------|--------|
| Discovery (BI/SCI) | Ali Sharbatian | `discovery/` |
| Validation (Neural Implant) | Kevin Joseph | `scripts/` |

---

## Repository Structure
```
.
├── discovery/                   # Discovery Analysis (BI/SCI)
│   ├── BI/                      # Brain Injury scripts
│   │   ├── 01_bi_fetch_preprocess.R
│   │   ├── 02_bi_wgcna.R
│   │   ├── 03_bi_diffexpr_gsea.R
│   │   └── 04_bi_visualizations.R
│   └── SCI/                     # Spinal Cord Injury scripts
│       ├── 01_sci_fetch_preprocess.R
│       ├── 02_sci_wgcna.R
│       ├── 03_sci_diffexpr_gsea.R
│       └── 04_sci_visualizations.R
├── scripts/                     # Validation Analysis (Neural Implant)
│   ├── 01_data_processing.R
│   ├── 02_annotation_mapping.R
│   ├── 03_differential_expression.R
│   ├── 04_axis_scoring.R
│   ├── 05_temporal_dynamics.R
│   ├── 06_module_preservation.R
│   ├── 07_pathway_enrichment.R
│   ├── config.R
│   └── theme_publication.R
├── docs/                        # Script documentation
├── results/                     # Generated outputs [created by scripts]
└── figures/                     # Generated plots [created by scripts]
```

---

## Part 1: Discovery Analysis (BI/SCI)

WGCNA and pathway analysis to identify ECM axes in neural injury models.

### Data Sources
| Dataset | Model | Platform |
|---------|-------|----------|
| GSE35338 | Brain Injury (MCAO) | Affymetrix Mouse 430 2.0 |
| GSE5296 | Spinal Cord Injury | Affymetrix Mouse 430 2.0 |

### Software (R 4.4.2)
```
WGCNA 1.73, limma 3.62.1, GEOquery 2.74.0, fgsea 1.32.2,
clusterProfiler 4.14.4, org.Mm.eg.db 3.20.0
```

### Running Discovery Analysis
```bash
# Brain Injury
Rscript discovery/BI/01_bi_fetch_preprocess.R
Rscript discovery/BI/02_bi_wgcna.R
Rscript discovery/BI/03_bi_diffexpr_gsea.R
Rscript discovery/BI/04_bi_visualizations.R

# Spinal Cord Injury
Rscript discovery/SCI/01_sci_fetch_preprocess.R
Rscript discovery/SCI/02_sci_wgcna.R
Rscript discovery/SCI/03_sci_diffexpr_gsea.R
Rscript discovery/SCI/04_sci_visualizations.R
```

---

## Part 2: Validation Analysis (Neural Implant)

Testing whether ECM axes identified in BI/SCI are also engaged after neural probe implantation.

### Biological Question
- **Do the ECM axes (especially the Hyaluronan axis) activate after implantation?**
- **How do those responses evolve over time** (acute → subacute → chronic)?

### Study Design
| Parameter | Value |
|-----------|-------|
| Species | Female Sprague-Dawley rats |
| Platform | Affymetrix Clariom S Rat arrays |
| Design | Paired (implant vs contralateral sham control) |
| Implant | Flexible polyimide neural probes (2mm depth, cortical) |
| Samples | n = 63 (31 implant, 32 control) |
| Timepoints | Day 0 (4h), 7, 14, 28, 126 |

### Software (R 4.3.0)
```
oligo 1.62.2, limma 3.54.0, GSVA 1.46.0, fgsea 1.24.0,
clusterProfiler, org.Rn.eg.db
```

### Running Validation Analysis
The pipeline is **linear**: each script writes outputs used by the next script. Run them in order:
```bash
cd scripts/
Rscript 01_data_processing.R      # CEL file loading, RMA normalization, QC
Rscript 02_annotation_mapping.R   # Probe-to-gene mapping, ortholog translation
Rscript 03_differential_expression.R  # Implant vs control DE at each timepoint
Rscript 04_axis_scoring.R         # ssGSEA pathway scoring
Rscript 05_temporal_dynamics.R    # Gene temporal classification
Rscript 06_module_preservation.R  # Module eigengene / preservation analysis
Rscript 07_pathway_enrichment.R   # fGSEA + literature signature validation
```

### Data Requirements (Validation)
This pipeline expects **raw Affymetrix CEL files** arranged as below:
```
Data/arrays/
├── Controls/
│   ├── Control_1.CEL ... Control_4.CEL   # Baseline
│   ├── Week0_1.CEL ... Week0_5.CEL       # Day 0 (4h)
│   ├── Week1_1.CEL ... Week1_5.CEL       # Day 7
│   ├── Week2_1.CEL ... Week2_6.CEL       # Day 14
│   ├── Week4_1.CEL ... Week4_6.CEL       # Day 28
│   └── Week18_1.CEL ... Week18_6.CEL     # Day 126
└── Implants/
    └── [same structure]
```

### "What should I look at?" (paper-oriented guide)
- **Axis-level activation across time**: `results/axis_scoring/axis_activation_stats.csv`
- **Differential expression (per timepoint)**: `results/deg/deg_significant_Week*.csv`
- **Pathway/gene set enrichment**: `results/ha_analysis/gsea_all_results.csv`
- **Temporal pattern classes**: `results/temporal/gene_temporal_classification.csv`
- **Module-level summaries**: `results/preservation/module_preservation_summary.csv`

### Key Outputs
| Folder | Contents |
|--------|----------|
| `results/deg/` | Differential expression results |
| `results/axis_scoring/` | ssGSEA scores and statistics |
| `results/ha_analysis/` | Pathway enrichment results |
| `results/temporal/` | Gene classification |
| `results/preservation/` | Module preservation |
| `figures/individual/` | Individual plots from each script |
| `figures/supplementary/` | QC, all axes, GSEA heatmaps |

---

## ECM Axes Analyzed

| Axis | Function | Key Genes |
|------|----------|-----------|
| **Hyaluronan** | HA synthesis/degradation, DAMP signaling (CD44, TLR2/4) | Has1-3, Cd44, Tlr2, Tlr4, Cd14 |
| **Provisional Matrix** | Fibronectin/tenascin scaffold, integrin engagement | Fn1, Tnc, Itga5, Itgb1 |
| **PNN-CSPG** | Perineuronal nets, chondroitin sulfate proteoglycans | Acan, Bcan, Vcan, Tnr |
| **Basement Membrane** | Laminin/collagen IV, blood-brain barrier | Hspg2, Lama4, Col4a1, Col4a2 |
| **Proteases/Regulators** | MMPs, ADAMTSs, TIMPs | Adamts4, Adamts5, Mmp9, Timp1 |
| **Crosslinking/Fibrosis** | LOX enzymes, fibrillar collagens | Lox, Loxl2, Col1a1, Tgm2 |

---

## Documentation

See `docs/` for detailed documentation of each script, including inputs, outputs, and methods.

---

## Glossary

- **DEG**: Differentially expressed gene (implant vs control)
- **FDR**: False discovery rate (multiple testing-adjusted p-value)
- **logFC**: Log2 fold change (positive = higher in implant)
- **ssGSEA / axis score**: Per-sample gene set score summarizing coordinated expression of an axis
- **NES**: Normalized enrichment score from gene set enrichment analysis (fGSEA)
- **WGCNA**: Weighted gene co-expression network analysis

---

## License

MIT License

Copyright (c) 2025 Ali Sharbatian, Kevin Joseph
