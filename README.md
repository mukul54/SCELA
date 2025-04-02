# SCELA: Single-Cell Expression Landscape of Lung Adenocarcinoma

A comprehensive pipeline for analyzing single-cell RNA sequencing data from lung cancer samples, with demographic stratification by age, sex, and smoking status.

## Overview

This pipeline processes the GSE131907 dataset, which contains single-cell RNA-seq data from lung adenocarcinoma patients. The analysis includes:

- Data preprocessing and quality control
- Cell type annotation based on known marker genes
- Differential expression analysis between tumor and normal samples
- Demographic comparisons (age, sex, smoking status)
- Visualization of gene expression patterns
- Identification of demographic-specific biomarkers

## Project Structure

```
SCELA/
├── gse131907_data_processing.py      # Primary data processing script
├── data_analysis/                    # Analysis scripts for visualization and exploration
│   ├── age/                          # Age-specific visualizations
│   │   └── age_tumor_analysis_post.py
│   ├── sex/                          # Sex-specific visualizations
│   │   ├── cell_specific_marker_violins.py
│   │   ├── improved_volcano_plots.py
│   │   ├── marker_violins_simple.py
│   │   ├── sex_marker_individual_plots.py
│   │   ├── sex_marker_volcano_plots.py
│   │   └── sex_specific_celltype_analysis.py
│   └── smoking/                      # Smoking-specific visualizations
│       ├── opposite_regulation_viz.py
│       ├── opposite_regulation_violin_plots.py
│       └── smoking_tumor_analysis.py
├── src/                              # Source code for core analyses
│   ├── age/                          # Age-specific differential expression analysis
│   │   ├── age_specific_degs.py
│   │   ├── age_stratified_analysis.py
│   │   ├── age_tumor_normal_degs.py
│   │   ├── basic_test_age.py
│   │   └── further_age_analysis.py
│   ├── sex/                          # Sex-specific differential expression analysis
│   │   ├── sex_specific_degs.py
│   │   └── sex_tumor_interaction.py
│   └── smoking/                      # Smoking-specific differential expression analysis
│       ├── smoking_specific_degs.py
│       └── smoking_tumor_interaction.py
├── age_analysis/                     # Additional age-related analyses and results
├── sex_analysis/                     # Additional sex-related analyses and results
├── smoking_analysis/                 # Additional smoking-related analyses and results
├── figures/                          # Generated figures and plots
└── data/                             # Data directory (not included in repo, download from hf)
    └── GSE131907/                    # Contains raw and processed data files
        ├── GSE131907_raw_UMI_matrix.txt
        ├── GSE131907_cell_annotation.txt.gz
        ├── GSE131907_metadata.csv
        ├── GSE131907_processed.h5
        └── age_celltype_annotation.h5ad.h5ad
```

## Requirements

- Python 3.10+
- Conda environment

## Environment Setup

Create and activate a new conda environment:

```bash
# Create conda environment
conda create -n scela python=3.10
conda activate scela

# Install the required packages
pip install -r requirements.txt
```

## Dependencies

The following Python packages are required:

```
# requirements.txt
scanpy==1.9.3
pandas==1.5.3
numpy==1.24.3
matplotlib==3.7.1
seaborn==0.12.2
scipy==1.10.1
anndata==0.8.0
scikit-learn==1.2.2
statsmodels==0.13.5
h5py==3.8.0
tables==3.8.0
adjustText==0.8
```

## Data Processing

To process the raw GSE131907 dataset:

```bash
python gse131907_data_processing.py
```

This script converts the raw UMI matrix to HDF5 format and creates an AnnData object with metadata.

## Demographic-Specific Analyses

### Sex-Specific Analysis

Analyze differential gene expression between males and females, stratified by cell types:

```bash
# Run sex-specific DEG analysis
python src/sex/sex_specific_degs.py --output_dir results/sex_degs

# Visualize sex-specific markers across cell types
python data_analysis/sex/sex_specific_celltype_analysis.py --base_dir results/sex_degs --output_dir figures/sex_analysis
```

### Age-Specific Analysis

Analyze differential gene expression across age groups:

```bash
# Run age-stratified analysis
python src/age/age_tumor_normal_degs.py

# Generate age-specific visualization
python data_analysis/age/age_tumor_analysis_post.py
```

### Smoking-Specific Analysis

Analyze differential gene expression between smokers and non-smokers:

```bash
# Run smoking-specific DEG analysis
python src/smoking/smoking_specific_degs.py

# Visualize smoking-specific patterns
python data_analysis/smoking/smoking_tumor_analysis.py
```

## Key Features

- Modular, well-organized codebase for demographic-stratified analyses
- Identification of sex-specific gene expression differences in cancer risk
- Age-stratified analysis of tumor vs normal tissues
- Smoking-related gene expression changes in lung adenocarcinoma
- Publication-ready visualizations with proper formatting
- Command-line interfaces for flexible parameter configuration
- Comprehensive cell type annotation based on literature-derived marker genes

## Output

Each analysis module generates the following outputs:

- Differential expression tables (`*_significant.csv`)
- Volcano plots showing significantly regulated genes
- Cell type-specific expression visualizations
- Summary statistics of demographic differences
- Interaction effect analyses between demographics and tumor/normal status
