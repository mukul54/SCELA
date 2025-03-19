"""
Data loading utilities for cancer single-cell RNA-seq analysis
"""

import os
import scanpy as sc
import pandas as pd
import numpy as np
from .config import DATA_DIR

def load_data():
    """
    Load AnnData object and integrate metadata
    """
    print("\n=== Loading Data ===")
    
    # Path to preprocessed anndata file
    anndata_file = os.path.join(DATA_DIR, "GSE131907_anndata.h5ad")
    
    if not os.path.exists(anndata_file):
        print(f"Error: Anndata file not found at {anndata_file}")
        print("Please run the data preprocessing script first.")
        return None
    
    # Load the preprocessed data
    print(f"Loading AnnData from {anndata_file}...")
    adata = sc.read_h5ad(anndata_file)
    
    # Load and integrate metadata
    print("Loading and integrating metadata...")
    metadata_file = os.path.join(DATA_DIR, "GSE131907_filtered_metadata.csv")
    if os.path.exists(metadata_file):
        metadata = pd.read_csv(metadata_file, index_col=0)
        
        # Merge metadata with existing obs
        for col in metadata.columns:
            if col not in adata.obs.columns:
                adata.obs[col] = metadata.loc[adata.obs.index, col].values
    
    # Print dataset dimensions
    print(f"Dataset dimensions: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Show available metadata columns
    print(f"Available metadata: {', '.join(adata.obs.columns)}")
    
    return adata
