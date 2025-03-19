# SCELA: Single-Cell Expression Landscape of Lung Adenocarcinoma

A modular pipeline for analyzing single-cell RNA sequencing data from lung cancer samples, with a focus on cell type annotation, marker gene identification, and demographic comparisons.

## Overview

This pipeline processes the GSE131907 dataset, which contains single-cell RNA-seq data from lung cancer patients. The analysis includes:

- Data preprocessing and quality control
- Dimensionality reduction (PCA, UMAP)
- Unsupervised clustering with the Leiden algorithm
- Marker gene identification for each cluster
- Cell type annotation based on known marker genes
- Demographic comparison between tumor and normal samples
- Publication-quality visualizations

## Requirements

- Python 3.8+
- Conda environment (environment name: "cancer")

## Dependencies

- scanpy
- pandas
- numpy
- matplotlib
- scikit-learn
- anndata
- scipy
- seaborn

## Project Structure

```
SCELA/
├── main.py                  # Main entry point for the pipeline
├── data_analysis/         # Core analysis modules
│   ├── __init__.py
│   ├── config.py            # Configuration parameters
│   ├── data_loader.py       # Data loading functions
│   ├── marker_analysis.py   # Marker gene identification
│   ├── cell_annotation.py   # Cell type annotation
│   ├── demographic_analysis.py  # Demographic comparisons
│   └── visualization.py     # Visualization utilities
└── data/                    # Data directory (not included in repo)
    └── GSE131907/
        ├── GSE131907_anndata.h5ad  # Input data file
        └── analysis_results/       # Output directory
            ├── cluster_markers.csv
            └── visualizations/     # Output visualizations
```

## Usage

1. Activate the conda environment:

   ```
   conda activate <env_name>
   ```
2. Run the main analysis script:

   ```
   python main.py
   ```

## Output

The pipeline generates the following outputs in the `data/GSE131907/analysis_results/` directory:

- `cluster_markers.csv`: Top marker genes for each identified cluster
- `visualizations/`: Directory containing:
  - `cell_types_umap.pdf`: UMAP visualization colored by cell types
  - `tumor_normal_comparison.pdf`: Comparison between tumor and normal samples
  - `sample_distribution.pdf`: Distribution of cells across sample types

## Features

- Modular, well-organized codebase for easy maintenance
- Publication-ready visualizations with proper formatting
- Flexible parameters configurable in `config.py`
- Comprehensive cell type annotation based on literature-derived marker genes
- Resolution of dataframe fragmentation warnings
- Improved clustering using the Leiden algorithm with igraph flavor
