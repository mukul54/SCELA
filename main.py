#!/usr/bin/env python3
"""
Main entry point for GSE131907 lung cancer cell type annotation pipeline.
This modular implementation organizes the analysis into separate components
for better maintainability and debugging.
"""

import os
import sys
import argparse
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import traceback

# Import analysis modules
from data_analysis.config import (RESULTS_DIR, VIS_DIR, 
                              FORCE_RECOMPUTE_PREPROCESSING, FORCE_RECOMPUTE_CLUSTERING,
                              FORCE_RECOMPUTE_MARKERS, FORCE_RECOMPUTE_DEMOGRAPHICS,
                              FORCE_RECOMPUTE_VISUALIZATION, MARKER_GENE_RESOLUTIONS)
from data_analysis.data_loader import load_data
from data_analysis.marker_analysis import identify_marker_genes, visualize_marker_genes
from data_analysis.cell_annotation import annotate_cell_types
from data_analysis.demographic_analysis import compare_demographics
from data_analysis.visualization import generate_visualizations

def main():
    """
    Main function for the lung cancer cell type annotation pipeline
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='GSE131907 Lung Cancer Cell Type Annotation Pipeline')
    parser.add_argument('--force-preprocess', action='store_true', default=FORCE_RECOMPUTE_PREPROCESSING,
                        help='Force recomputation of preprocessing steps')
    parser.add_argument('--force-clustering', action='store_true', default=FORCE_RECOMPUTE_CLUSTERING,
                        help='Force recomputation of clustering results')
    parser.add_argument('--force-markers', action='store_true', default=FORCE_RECOMPUTE_MARKERS,
                        help='Force recomputation of marker genes')
    parser.add_argument('--force-demographics', action='store_true', default=FORCE_RECOMPUTE_DEMOGRAPHICS,
                        help='Force recomputation of demographic analysis')
    parser.add_argument('--force-visualizations', action='store_true', default=FORCE_RECOMPUTE_VISUALIZATION,
                        help='Force regeneration of visualizations')
    parser.add_argument('--force-all', action='store_true', default=False,
                        help='Force recomputation of all analysis steps')
    args = parser.parse_args()
    
    # If force-all is specified, set all other force flags to True
    if args.force_all:
        args.force_preprocess = True
        args.force_clustering = True
        args.force_markers = True
        args.force_demographics = True
        args.force_visualizations = True
    
    print("\n=== GSE131907 Lung Cancer Cell Type Annotation Pipeline ===\n")
    
    # Set pandas display options and warning handling
    pd.set_option('display.max_columns', 100)
    pd.set_option('display.max_rows', 100)
    
    # Create output directories if they don't exist
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(VIS_DIR, exist_ok=True)
    
    # Check if annotated data already exists
    annotated_file = f"{RESULTS_DIR}/GSE131907_annotated.h5ad"
    preprocessed_file = f"{RESULTS_DIR}/GSE131907_preprocessed.h5ad"
    
    # Determine whether to load existing data or reprocess
    if os.path.exists(annotated_file) and not (args.force_preprocess or args.force_clustering):
        print(f"Loading existing annotated data from {annotated_file}")
        adata = sc.read_h5ad(annotated_file)
    elif os.path.exists(preprocessed_file) and not args.force_preprocess:
        print(f"Loading existing preprocessed data from {preprocessed_file}")
        adata = sc.read_h5ad(preprocessed_file)
    else:
        # Load and preprocess the data from scratch
        print("Loading and preprocessing data from raw files...")
        adata = load_data()
        if adata is None:
            return
        
        # Save preprocessed data if it doesn't exist
        if not os.path.exists(preprocessed_file) or args.force_preprocess:
            print("Saving preprocessed data...")
            adata.write(preprocessed_file)
    
    try:
        # Only run clustering if needed
        if not os.path.exists(annotated_file) or args.force_clustering:
            print("Running clustering and cell type annotation...")
            # Make sure columns are properly typed before annotation
            for col in adata.obs.columns:
                if col.endswith('_type') and not isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
                    print(f"Converting {col} to categorical type")
                    adata.obs[col] = adata.obs[col].astype('category')
            
            # Annotate cell types (which includes clustering)
            adata = annotate_cell_types(adata)
            
            # Ensure cell type columns are categorical after annotation
            for col in ['cell_type', 'predicted_cell_type']:
                if col in adata.obs.columns and not isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
                    print(f"Converting {col} to categorical type after annotation")
                    adata.obs[col] = adata.obs[col].astype('category')
            
            # Save the annotated data
            print("\nSaving annotated data...")
            adata.write(annotated_file)
        
        # Identify marker genes if needed
        if args.force_markers or not os.path.exists(f"{RESULTS_DIR}/marker_genes_res{MARKER_GENE_RESOLUTIONS[0]}.csv"):
            print("Identifying marker genes...")
            adata = identify_marker_genes(adata)
        
        # Run demographic analysis if needed or requested
        if args.force_demographics or not os.path.exists(f"{VIS_DIR}/gender_comparison.pdf"):
            print("Performing demographic analysis...")
            adata = compare_demographics(adata)
        
        # Generate visualizations if needed or requested
        if args.force_visualizations or not os.path.exists(f"{VIS_DIR}/umap.pdf"):
            print("Generating visualizations...")
            generate_visualizations(adata)
        
        # Visualize marker genes if needed or requested
        if args.force_visualizations or not os.path.exists(f"{VIS_DIR}/marker_heatmap.pdf"):
            print("Visualizing marker genes...")
            visualize_marker_genes(adata)
        
        print("\n=== Analysis Pipeline Completed Successfully! ===")
    except Exception as e:
        print(f"\nError during analysis: {e}")
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
