# Neural Implant Transcriptomics Validation Analysis

## Background

Neural implants trigger a foreign body response (FBR) that can compromise device performance over time. Understanding the molecular dynamics of this response is critical for developing biocompatible neural interfaces.

Previous work analyzing brain injury (BI) and spinal cord injury (SCI) transcriptomics identified six extracellular matrix (ECM) axes that define tissue remodeling phases:

| Axis | Function |
|------|----------|
| **Hyaluronan** | HA synthesis/degradation, DAMP receptor signaling (CD44, TLR2/4) |
| **Provisional Matrix** | Fibronectin/tenascin scaffold, integrin engagement |
| **PNN-CSPG** | Perineuronal nets, chondroitin sulfate proteoglycans |
| **Basement Membrane** | Laminin/collagen IV, blood-brain barrier |
| **Proteases/Regulators** | MMPs, ADAMTSs, TIMPs |
| **Crosslinking/Fibrosis** | LOX enzymes, fibrillar collagens |

This analysis validates whether these ECM axes—particularly the Hyaluronan axis—are activated in response to flexible neural probe implantation.

## Study Design

| Parameter | Value |
|-----------|-------|
| Species | Female Sprague-Dawley rats |
| Platform | Affymetrix Clariom S Rat arrays |
| Design | Paired (implant vs contralateral sham control) |
| Implant | Flexible polyimide neural probes (2mm depth, cortical) |
| Samples | n = 63 (31 implant, 32 control) |
| Timepoints | Day 0 (4h), 7, 14, 28, 126 |

The paired design (implant and control from the same animal) increases statistical power by controlling for inter-animal variability.

## Repository Structure

```
analysis/
├── README.md                    # This file
├── config.R                     # File paths and parameters
├── run_analysis.sh              # Pipeline execution script
├── scripts/
│   ├── 01_data_processing.R     # CEL file loading, RMA normalization, QC
│   ├── 02_annotation_mapping.R  # Probe-to-gene mapping, ortholog translation
│   ├── 03_differential_expression.R  # Paired limma analysis
│   ├── 04_axis_scoring.R        # ssGSEA pathway scoring
│   ├── 05_temporal_dynamics.R   # Gene temporal classification
│   ├── 06_module_preservation.R # WGCNA module eigengene analysis
│   ├── 07_pathway_enrichment.R  # fGSEA, literature signature validation
│   └── 08_figure_generation.R   # Main and supplementary figures
├── R/
│   └── theme_publication.R      # ggplot2 theme and color palettes
├── docs/                        # Script documentation
├── results/                     # Analysis outputs (.csv, .rds)
└── figures/                     # Generated plots (.pdf, .png)
```

## Running the Pipeline

### Prerequisites

```r
# Bioconductor packages
BiocManager::install(c(
  "oligo", "pd.clariom.s.rat", "clariomsrattranscriptcluster.db",
  "limma", "GSVA", "fgsea", "clusterProfiler", "org.Rn.eg.db"
))

# CRAN packages
install.packages(c("tidyverse", "patchwork", "pheatmap"))
```

### Execution

Scripts must be run in order (01 → 08) as each depends on outputs from previous scripts.

## Data Requirements

CEL files should be organized as:

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

## Outputs

### Main Figure
- `figures/main/Figure_Main.pdf` — 6-panel figure

### Supplementary Figures
- `figures/supplementary/` — QC, all axes, GSEA heatmaps

### Data Tables
- `results/deg/` — Differential expression results
- `results/axis_scoring/` — ssGSEA scores and statistics
- `results/ha_analysis/` — Pathway enrichment results
- `results/temporal/` — Gene classification
- `results/preservation/` — Module preservation

## Documentation

See `docs/` for detailed documentation of each script, including inputs, outputs, and methods.

## License

MIT License
