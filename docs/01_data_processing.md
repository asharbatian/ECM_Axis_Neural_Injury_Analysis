# 01_data_processing.R

## Purpose

Loads raw Affymetrix CEL files, performs RMA normalization, and conducts quality control assessment of the neural implant transcriptomics dataset.

## Input

### Required Files
- CEL files in `Data/arrays/Controls/` (sham controls)
- CEL files in `Data/arrays/Implants/` (implanted samples)

### Expected File Structure
```
Data/arrays/
├── Controls/
│   ├── Control_1.CEL - Control_4.CEL    # Baseline
│   ├── Week0_1.CEL - Week0_5.CEL        # Day 0 (4h)
│   ├── Week1_1.CEL - Week1_5.CEL        # Day 7
│   ├── Week2_1.CEL - Week2_6.CEL        # Day 14
│   ├── Week4_1.CEL - Week4_6.CEL        # Day 28
│   └── Week18_1.CEL - Week18_6.CEL      # Day 126
└── Implants/
    └── [same structure]
```

## Output

### Data Files
| File | Description |
|------|-------------|
| `results/expr_matrix_normalized.rds` | RMA-normalized probe-level expression matrix |
| `results/sample_metadata.csv` | Sample annotations (condition, timepoint, animal ID) |
| `results/qc/pca_coordinates.csv` | PCA coordinates for all samples |
| `results/qc/sample_correlations.csv` | Pairwise sample correlation matrix |

### Figures
| File | Description |
|------|-------------|
| `figures/01_pca_all_samples.pdf` | PCA plot colored by condition and timepoint |

## Methods

### Normalization
- **Algorithm**: Robust Multi-array Average (RMA)
- **Steps**: Background correction → Quantile normalization → Median polish summarization
- **Package**: `oligo` (v1.62+)

### Sample Naming Convention
Samples are renamed to avoid duplicate filenames:
- Control samples: `Control_{original_filename}`
- Implant samples: `Implant_{original_filename}`

### Metadata Extraction
Metadata is parsed from filenames:
- **Timepoint**: Extracted from filename (Week0, Week1, etc.)
- **Condition**: Derived from directory (Controls → Control, Implants → Implant)
- **Animal ID**: Constructed from timepoint and replicate number
- **Pair ID**: Links implant/control samples from same animal

### Quality Control
1. **PCA**: Principal component analysis to visualize sample clustering
2. **Correlation**: Pairwise Pearson correlation between samples
3. **Clustering**: Hierarchical clustering to identify outliers

## Key Parameters

```r
# RMA normalization (default parameters)
expr_matrix <- rma(raw_data)

# No samples excluded based on QC
# Threshold for outlier detection: correlation < 0.85 (visual inspection)
```

## Dependencies

```r
library(oligo)           # CEL file processing
library(tidyverse)       # Data manipulation
library(ggplot2)         # Visualization
```

