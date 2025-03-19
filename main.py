#!/usr/bin/env python3
"""
Main entry point for GSE131907 lung cancer cell type annotation pipeline.
This modular implementation organizes the analysis into separate components
for better maintainability and debugging.
"""

import os
import sys
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import traceback

# Import analysis modules
from data_analysis.config import RESULTS_DIR, VIS_DIR
from data_analysis.data_loader import load_data
from data_analysis.marker_analysis import identify_marker_genes, visualize_marker_genes
from data_analysis.cell_annotation import annotate_cell_types
from data_analysis.demographic_analysis import compare_demographics
from data_analysis.visualization import generate_visualizations

def main():
    """
    Main function for the lung cancer cell type annotation pipeline
    """
    print("\n=== GSE131907 Lung Cancer Cell Type Annotation Pipeline ===\n")
    
    # Set pandas display options and warning handling
    pd.set_option('display.max_columns', 100)
    pd.set_option('display.max_rows', 100)
    
    # Create output directories if they don't exist
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(VIS_DIR, exist_ok=True)
    
    # Load the data
    adata = load_data()
    if adata is None:
        return
    
    try:
        # Identify marker genes for each cluster
        adata = identify_marker_genes(adata)
        
        # Make sure columns are properly typed before annotation
        for col in adata.obs.columns:
            if col.endswith('_type') and not isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
                print(f"Converting {col} to categorical type")
                adata.obs[col] = adata.obs[col].astype('category')
        
        # Annotate cell types
        adata = annotate_cell_types(adata)
        
        # Ensure cell type columns are categorical after annotation
        for col in ['cell_type', 'predicted_cell_type']:
            if col in adata.obs.columns and not isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
                print(f"Converting {col} to categorical type after annotation")
                adata.obs[col] = adata.obs[col].astype('category')
        
        # Compare demographics
        adata = compare_demographics(adata)
        
        # Generate visualizations
        generate_visualizations(adata)
        
        # Visualize marker genes
        visualize_marker_genes(adata)
        
        # Save the annotated data
        print("\nSaving annotated data...")
        adata.write(f"{RESULTS_DIR}/GSE131907_annotated.h5ad")
        
        print("\n=== Analysis Pipeline Completed Successfully! ===")
    except Exception as e:
        print(f"\nError during analysis: {e}")
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
