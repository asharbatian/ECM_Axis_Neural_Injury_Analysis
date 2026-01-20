# HA Axis Validation Study — Neural Implant Transcriptomics (Paper Companion)

## What this repository is

This repository is a **paper companion** that reproduces the core transcriptomic analyses used to test whether **extracellular matrix (ECM) remodeling “axes”**—identified previously in **brain injury (BI)** and **spinal cord injury (SCI)**—are also engaged after **flexible neural probe implantation**.

If you are here to **understand the paper**, start with:

- **Main figure (compiled)**: `figures/main/Figure_Main.pdf`
- **Supplementary figures**: `figures/supplementary/`
- **Key result tables**: `results/` (subfolders described below)

## Biological question

Neural implants trigger a foreign body response (FBR) that can compromise device performance over time. Here, we ask:

- **Do the ECM axes (especially the Hyaluronan axis) activate after implantation?**
- **How do those responses evolve over time** (acute → subacute → chronic)?

Analysis of brain injury (BI) and spinal cord injury (SCI) transcriptomics identified six extracellular matrix (ECM) axes that define tissue remodeling phases:

| Axis | Function |
|------|----------|
| **Hyaluronan** | HA synthesis/degradation, DAMP receptor signaling (CD44, TLR2/4) |
| **Provisional Matrix** | Fibronectin/tenascin scaffold, integrin engagement |
| **PNN-CSPG** | Perineuronal nets, chondroitin sulfate proteoglycans |
| **Basement Membrane** | Laminin/collagen IV, blood-brain barrier |
| **Proteases/Regulators** | MMPs, ADAMTSs, TIMPs |
| **Crosslinking/Fibrosis** | LOX enzymes, fibrillar collagens |

This pipeline quantifies those axes in implant vs contralateral control tissue across time, then summarizes findings in publication-ready figures.

## Study Design

| Parameter | Value |
|-----------|-------|
| Species | Female Sprague-Dawley rats |
| Platform | Affymetrix Clariom S Rat arrays |
| Design | Paired (implant vs contralateral sham control) |
| Implant | Flexible polyimide neural probes (2mm depth, cortical) |
| Samples | n = 63 (31 implant, 32 control) |
| Timepoints | Day 0 (4h), 7, 14, 28, 126 |

## “What should I look at?” (paper-oriented guide)

- **Axis-level activation across time**: `results/axis_scoring/axis_activation_stats.csv`
- **Differential expression (per timepoint)**: `results/deg/deg_significant_Week*.csv`
- **Pathway/gene set enrichment (including literature signatures)**: `results/ha_analysis/gsea_all_results.csv`
- **Temporal pattern classes (resolving/persistent/late, etc.)**: `results/temporal/gene_temporal_classification.csv`
- **Module-level summaries (module eigengenes / preservation-style summaries)**: `results/preservation/module_preservation_summary.csv`

## Repository Structure

```
.
├── README.md                    # This file
├── scripts/
│   ├── 01_data_processing.R     # CEL file loading, RMA normalization, QC
│   ├── 02_annotation_mapping.R  # Probe-to-gene mapping, ortholog translation
│   ├── 03_differential_expression.R  # Unpaired limma analysis
│   ├── 04_axis_scoring.R        # ssGSEA pathway scoring
│   ├── 05_temporal_dynamics.R   # Gene temporal classification
│   ├── 06_module_preservation.R # Module eigengene / preservation analysis
│   ├── 07_pathway_enrichment.R  # fGSEA + literature signature validation
│   ├── config.R                 # Shared paths/parameters used by scripts
│   └── theme_publication.R      # Plot theme + save helper used by scripts
├── docs/                        # Script documentation (inputs/outputs/methods)
├── results/                     # Generated analysis outputs (.csv, .rds) [created by scripts]
└── figures/                     # Generated plots (.pdf, .png) [created by scripts]
```

## Reproducing the analysis (01 → 08)

The pipeline is **linear**: each script writes outputs used by the next script. Run them in order:

- **`scripts/01_data_processing.R`**: normalize microarray CEL files (RMA) + QC + sample metadata
- **`scripts/02_annotation_mapping.R`**: map probes → genes; build gene-level matrix; create ECM gene sets
- **`scripts/03_differential_expression.R`**: implant vs control differential expression at each timepoint
- **`scripts/04_axis_scoring.R`**: compute axis scores per sample (ssGSEA) + stats across time
- **`scripts/05_temporal_dynamics.R`**: classify genes/axes by temporal behavior across the full time course
- **`scripts/06_module_preservation.R`**: module eigengene/activity summaries using manuscript-derived modules
- **`scripts/07_pathway_enrichment.R`**: gene set enrichment and literature signature validation

For step-by-step inputs/outputs, see `docs/` (one markdown page per script).

## Data requirements (what you need to provide)

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

## Running the Pipeline (R packages)

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

## Outputs

### Main Figure plots
- `figures/individual/` -  individual plots from each script

### Supplementary Figure plots
- `figures/supplementary/` — QC, all axes, GSEA heatmaps

### Data Tables
- `results/deg/` — Differential expression results
- `results/axis_scoring/` — ssGSEA scores and statistics
- `results/ha_analysis/` — Pathway enrichment results
- `results/temporal/` — Gene classification
- `results/preservation/` — Module preservation

## Documentation

See `docs/` for detailed documentation of each script, including inputs, outputs, and methods.

## Glossary (quick definitions)

- **DEG**: differentially expressed gene (implant vs control)
- **FDR**: false discovery rate (multiple testing-adjusted p-value)
- **logFC**: log2 fold change (positive = higher in implant)
- **ssGSEA / axis score**: per-sample gene set score summarizing coordinated expression of an axis
- **NES**: normalized enrichment score from gene set enrichment analysis (fGSEA)

## License

MIT License
