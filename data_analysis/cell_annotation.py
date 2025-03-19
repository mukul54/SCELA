# cell_annotations.py
import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from .config import RESULTS_DIR, VIS_DIR

def annotate_cell_types(adata):
    """
    Annotate cell types based on Leiden clustering at multiple resolutions
    """
    print("\n=== Annotating Cell Types ===")
    
    # Calculate neighborhood graph if not already calculated
    if 'neighbors' not in adata.uns:
        print("Computing neighborhood graph...")
        sc.pp.neighbors(adata, n_neighbors=15)
    else:
        print("Using existing neighborhood graph")
    
    # Make sure we have a UMAP embedding
    if 'X_umap' not in adata.obsm:
        print("Computing UMAP embedding...")
        sc.tl.umap(adata)
    
    # Create visualizations directory
    os.makedirs(VIS_DIR, exist_ok=True)
    
    # List of resolutions to try
    resolutions = [0.25, 0.5, 1.0]
    leiden_keys = []
    
    # Run Leiden clustering at different resolutions
    for res in resolutions:
        leiden_key = f'leiden_{res}'
        leiden_keys.append(leiden_key)
        
        print(f"Computing Leiden clusters at resolution {res}...")
        sc.tl.leiden(adata, resolution=res, key_added=leiden_key)
        print(f"Generated {len(adata.obs[leiden_key].unique())} clusters at resolution {res}")
        
        # Make sure the column is categorical
        if not pd.api.types.is_categorical_dtype(adata.obs[leiden_key]):
            adata.obs[leiden_key] = adata.obs[leiden_key].astype('category')
            
        # Plot cluster distribution
        plt.figure(figsize=(10, 6))
        adata.obs[leiden_key].value_counts().sort_index().plot(kind='bar')
        plt.title(f'Leiden Cluster Distribution (Resolution {res})')
        plt.xlabel('Cluster ID')
        plt.ylabel('Number of Cells')
        plt.tight_layout()
        cluster_dist_file = f"{VIS_DIR}/cluster_distribution_res_{res}.pdf"
        plt.savefig(cluster_dist_file)
        plt.close()
        print(f"Saved cluster distribution plot to {cluster_dist_file}")
        
        # Create UMAP visualization for this resolution
        plt.figure(figsize=(12, 10))
        sc.pl.umap(
            adata, 
            color=leiden_key,
            legend_loc='on data',
            legend_fontweight='bold',
            legend_fontoutline=2,
            legend_fontsize=12,
            palette='tab20',
            title=f'Leiden Clusters (Resolution {res})',
            frameon=False,
            s=40,
            alpha=0.8,
            show=False
        )
        
        # Save the UMAP visualization
        cluster_umap_file = f"{VIS_DIR}/leiden_clusters_res_{res}_umap.pdf"
        plt.savefig(cluster_umap_file, bbox_inches='tight', dpi=150)
        plt.close()
        print(f"Saved Leiden clusters UMAP for resolution {res} to {cluster_umap_file}")
    
    # Use the middle resolution (0.5) as the default cell type annotation
    default_leiden = leiden_keys[1]  # This will be leiden_0.5
    
    # Set the default leiden clusters as the cell_type
    adata.obs['leiden'] = adata.obs[default_leiden]
    adata.obs['cell_type'] = adata.obs[default_leiden]
    adata.obs['cell_type_display'] = 'Cluster ' + adata.obs[default_leiden]
    
    print(f"Using Leiden clustering at resolution 0.5 as default cell type annotation.")
    print(f"Found {len(adata.obs['cell_type'].unique())} clusters")
    
    # Print statistics about the default clusters
    for cluster, count in adata.obs['cell_type'].value_counts().items():
        print(f"  Cluster {cluster}: {count} cells")
    
    # Identify marker genes for each cluster using the default resolution
    print("Identifying marker genes for default clusters...")
    try:
        # Calculate marker genes for the default clusters
        sc.tl.rank_genes_groups(adata, default_leiden, method='wilcoxon')
        print("Successfully identified marker genes for all clusters")
    except Exception as e:
        print(f"Error identifying marker genes: {e}")
    
    return adata
