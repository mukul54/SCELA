# RNA-Seq Data Analysis Pipeline

This folder contains a modular analysis pipeline for single-cell RNA sequencing (scRNA-seq) data, specifically designed for the GSE131907 lung cancer dataset. The pipeline is structured into several Python modules, each handling a specific aspect of the analysis workflow.

## Overview

This analysis pipeline processes scRNA-seq data through several steps:
1. Data loading and quality control
2. Cell type annotation using Leiden clustering
3. Marker gene identification and visualization
4. Demographic comparison of cell distributions
5. Comprehensive visualization generation

## Module Descriptions

### `config.py`

Configuration settings for the entire analysis pipeline.

**Key Settings:**
- Data paths (DATA_DIR, RESULTS_DIR, VIS_DIR)
- Analysis parameters:
  - N_PCS (50): Number of principal components to use
  - N_NEIGHBORS (15): Number of neighbors for neighborhood graph
  - RESOLUTION (0.5): Resolution for Leiden clustering
  - N_TOP_GENES (2000): Number of highly variable genes to select
  - TOP_N_MARKERS (5): Number of top marker genes to show per cell type

### `data_loader.py`

Handles loading of preprocessed AnnData objects and integration with metadata.

**Key Functions:**
- `load_data()`: Loads the preprocessed AnnData object from the h5ad file and integrates metadata from the CSV file.

### `cell_annotation.py`

Performs cell type annotation using the Leiden clustering algorithm.

**Key Functions:**
- `annotate_cell_types(adata)`: Annotates cell types based on Leiden clustering, with the following steps:
  1. Computes Leiden clusters if not already present
  2. Uses Leiden clusters as cell type annotations
  3. Identifies marker genes for each cluster using the Wilcoxon rank-sum test
  4. Generates cluster distribution plots
  5. Creates UMAP visualizations colored by Leiden clusters

### `marker_analysis.py`

Identifies and visualizes marker genes for each cell cluster.

**Key Functions:**
- `identify_marker_genes(adata)`: Performs clustering and identifies marker genes for each cluster
  - Computes PCA if not already done
  - Finds highly variable genes
  - Builds neighborhood graph
  - Performs Leiden clustering
  - Identifies marker genes using rank_genes_groups
  
- `visualize_marker_genes(adata, flavor)`: Generates visualizations for marker genes
  - Creates heatmaps of marker gene expression across clusters
  - Generates dot plots showing marker gene expression
  - Produces violin plots for top marker genes

### `demographic_analysis.py`

Analyzes and compares cell type distributions across different demographic groups.

**Key Functions:**
- `compare_demographics(adata)`: Compares cell type composition between demographic groups
  - Male vs. Female comparison
  - Age group comparison
  - Tumor vs. Normal tissue comparison
  - Creates stacked bar charts and heatmaps to visualize differences

### `visualization.py`

Creates comprehensive visualizations for cell type analysis.

**Key Functions:**
- `generate_visualizations(adata)`: Creates various plots for analysis
  - Cell type distribution across samples
  - UMAP visualizations with different metadata overlays
  - Gene expression overlays on UMAP plots
  - Cell density plots

## Workflow

The typical workflow for using this pipeline is:

1. **Load Data**: Use `data_loader.py` to load preprocessed AnnData
2. **Annotate Cell Types**: Apply `cell_annotation.py` to perform Leiden clustering and annotate cell types
3. **Identify Markers**: Use `marker_analysis.py` to identify characteristic genes for each cluster
4. **Demographic Analysis**: Run `demographic_analysis.py` to compare cell distributions across patient groups
5. **Generate Visualizations**: Create comprehensive visualizations with `visualization.py`

## Output

The pipeline generates various outputs in the configured directories:

- **RESULTS_DIR**: Contains analysis results including tables of marker genes
- **VIS_DIR**: Contains visualizations including:
  - Leiden cluster UMAP plots
  - Cluster distribution charts
  - Marker gene heatmaps and dot plots
  - Demographic comparison visualizations

## Usage Example

```python
from data_analysis.data_loader import load_data
from data_analysis.cell_annotation import annotate_cell_types
from data_analysis.marker_analysis import identify_marker_genes, visualize_marker_genes
from data_analysis.demographic_analysis import compare_demographics
from data_analysis.visualization import generate_visualizations

# Load the data
adata = load_data()

# Annotate cell types using Leiden clustering
adata = annotate_cell_types(adata)

# Identify marker genes
adata = identify_marker_genes(adata)

# Visualize marker genes
visualize_marker_genes(adata)

# Compare demographics
compare_demographics(adata)

# Generate comprehensive visualizations
generate_visualizations(adata)
```

## Notes

- The pipeline is designed to work with preprocessed data in AnnData format
- Leiden clustering is used for unsupervised identification of cell types
- The pipeline automatically handles cases where certain computations (e.g., PCA, UMAP) have already been performed
