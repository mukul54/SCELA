"""
Functions for marker gene identification and analysis
"""

import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Configure scanpy settings to avoid common warnings
# Use sc.settings directly instead of importing scanpy.settings
sc.settings.verbosity = 1  # Reduce verbosity 
sc.settings.set_figure_params(dpi=150)
from .config import RESULTS_DIR, VIS_DIR, N_PCS, N_NEIGHBORS, RESOLUTION, N_TOP_GENES, TOP_N_MARKERS

def identify_marker_genes(adata):
    """
    Identify marker genes for each cluster
    """
    print("\n=== Identifying Marker Genes ===")
    
    # Format axis labels for plots
    plt.rcParams['axes.titlesize'] = 14  # Larger title font
    plt.rcParams['axes.labelsize'] = 12  # Larger axis label font
    
    # Perform clustering if not already done
    if 'leiden' not in adata.obs:
        print("Performing clustering...")
        
        # Compute PCA if not already computed
        if 'X_pca' not in adata.obsm:
            print("Computing PCA...")
            # Scale the data for PCA if needed
            if 'log1p' not in adata.uns:
                print("Normalizing and log-transforming data...")
                sc.pp.normalize_total(adata, target_sum=1e4)
                sc.pp.log1p(adata)
                adata.uns['log1p'] = {'base': None}
                
            # Calculate highly variable genes if needed
            if 'highly_variable' not in adata.var.columns:
                print("Finding highly variable genes...")
                sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES, flavor='seurat')
                
            # Compute PCA
            sc.pp.pca(adata, use_highly_variable=True, svd_solver='arpack')
            print(f"Computed {adata.obsm['X_pca'].shape[1]} principal components")
            
        # Computing neighborhood graph
        print("Computing neighbor graph...")
        sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
        
        # Use igraph flavor to avoid future warnings
        print("Running Leiden clustering...")
        sc.tl.leiden(adata, resolution=RESOLUTION, flavor='igraph', n_iterations=2, directed=False)
        print(f"Found {len(adata.obs['leiden'].unique())} clusters")
    
    # Compute UMAP embedding if not already present
    if 'X_umap' not in adata.obsm:
        print("Computing UMAP embedding...")
        sc.tl.umap(adata)
    
    # Logarithmize data if not already done
    if 'log1p' not in adata.uns:
        print("Logarithmizing the data before finding marker genes...")
        sc.pp.log1p(adata)
        adata.uns['log1p'] = {'base': None}
    
    # Find markers for each cluster
    print("Finding marker genes for each cluster...")
    # Copy data to avoid performance warnings related to DataFrame fragmentation
    adata_copy = adata.copy() 
    sc.tl.rank_genes_groups(adata_copy, 'leiden', method='wilcoxon')
    # Copy results back to original adata
    adata.uns['rank_genes_groups'] = adata_copy.uns['rank_genes_groups'].copy()
    del adata_copy
    
    # Extract top markers into a dataframe
    print("Extracting top markers...")
    markers_df = sc.get.rank_genes_groups_df(adata, group=None)
    
    # Save the marker genes to file
    marker_file = os.path.join(RESULTS_DIR, "cluster_markers.csv")
    markers_df.to_csv(marker_file)
    print(f"Saved cluster markers to {marker_file}")
    
    return adata

def visualize_marker_genes(adata, flavor='igraph'):
    """
    Generate expression heatmaps for marker genes across cell types
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    flavor : str
        Leiden algorithm flavor to use (defaults to igraph to avoid warnings)
    """
    # Configure scanpy to use recommended parameters
    sc.settings.set_figure_params(dpi=150, frameon=False)
    # Use lower verbosity to reduce output noise
    verbosity_save = sc.settings.verbosity
    sc.settings.verbosity = 1
    print("\n=== Visualizing Marker Genes ===")
    
    # Set up figure parameters for better readability
    plt.rcParams['figure.figsize'] = (14, 10)  # Default figure size
    plt.rcParams['figure.dpi'] = 150          # Higher resolution figures
    
    # Cell type column to use
    cell_type_column = 'cell_type' if 'cell_type' in adata.obs.columns else 'predicted_cell_type'
    
    # Ensure cell_type_column is categorical - using isinstance to avoid deprecation warnings
    if not isinstance(adata.obs[cell_type_column].dtype, pd.CategoricalDtype):
        adata.obs[cell_type_column] = adata.obs[cell_type_column].astype('category')
        
    # Create a display-friendly version of the column name for plot titles
    cell_type_display = 'Cell Type' if cell_type_column == 'cell_type' else 'Predicted Cell Type'
    
    # Take top N genes from previously computed markers
    marker_genes_per_cluster = {}
    genes_to_show = []
    if 'rank_genes_groups' in adata.uns and adata.uns['rank_genes_groups']:
        try:
            # Create dictionary of marker genes grouped by cell type for better visualization
            for cluster in adata.obs[cell_type_column].cat.categories:
                # Find cells belonging to this cell type
                cluster_cells = adata.obs[cell_type_column] == cluster
                if not any(cluster_cells):
                    continue
                    
                # Find the most common Leiden cluster for this cell type
                leiden_counts = adata.obs.loc[cluster_cells, 'leiden'].value_counts()
                if leiden_counts.empty:
                    continue
                    
                leiden_cluster = leiden_counts.idxmax()
                
                # Get top marker genes for this leiden cluster
                top_genes = sc.get.rank_genes_groups_df(
                    adata, 
                    group=leiden_cluster,
                    key='rank_genes_groups'
                ).head(TOP_N_MARKERS)['names'].tolist()
                
                marker_genes_per_cluster[cluster] = top_genes
                genes_to_show.extend(top_genes)
            
            # De-duplicate genes
            genes_to_show = list(dict.fromkeys(genes_to_show))
            
        except Exception as e:
            print(f"Error extracting marker genes: {e}")
            # Use top genes by expression as fallback
            print("Using top genes by mean expression as fallback...")
            genes_to_show = adata.var_names[adata.var.mean_counts.argsort()[::-1][:50]].tolist()
    else:
        # Use top genes by expression as fallback
        print("No precomputed marker genes found. Using top genes by mean expression...")
        if 'highly_variable' in adata.var:
            genes_to_show = adata.var_names[adata.var.highly_variable].tolist()[:50]
        elif 'mean_counts' in adata.var:
            genes_to_show = adata.var_names[adata.var.mean_counts.argsort()[::-1][:50]].tolist()
        else:
            # Calculate gene stats if not available
            sc.pp.calculate_qc_metrics(adata, inplace=True)
            genes_to_show = adata.var_names[adata.var.mean_counts.argsort()[::-1][:50]].tolist()
    
    print(f"Creating expression heatmap for {len(genes_to_show)} genes across {len(adata.obs[cell_type_column].cat.categories)} cell types...")
    
    # Create a temporary column with formatted cell type names for plotting
    temp_col = f"{cell_type_column}_display"
    if cell_type_column in adata.obs:
        # Ensure the column is categorical before accessing .cat - using isinstance to avoid deprecation warnings
        if not isinstance(adata.obs[cell_type_column].dtype, pd.CategoricalDtype):
            adata.obs[cell_type_column] = adata.obs[cell_type_column].astype('category')
            
        # Format category names by replacing underscores with spaces and capitalizing properly
        categories = adata.obs[cell_type_column].cat.categories
        category_map = {cat: cat.replace('_', ' ').replace('-', ' ').title() for cat in categories}
        
        # Create formatted column
        adata.obs[temp_col] = adata.obs[cell_type_column].map(category_map)
        adata.obs[temp_col] = adata.obs[temp_col].astype('category')
        adata.obs[temp_col] = adata.obs[temp_col].cat.reorder_categories([category_map[cat] for cat in categories])
        plot_groupby = temp_col
    else:
        plot_groupby = cell_type_column
    
    # Format gene names if needed
    formatted_genes = [gene.replace('_', ' ') for gene in genes_to_show]
    
    # Create a copy of the data to avoid fragmentation warnings
    adata_copy = adata.copy()
    
    # Prepare for dendrogram first to avoid warning
    if f'dendrogram_{plot_groupby}' not in adata.uns:
        print(f"Computing dendrogram for {plot_groupby}...")
        sc.tl.dendrogram(adata, groupby=plot_groupby)
    
    # Adjust figure size based on number of genes and cell types
    n_cell_types = len(adata.obs[cell_type_column].cat.categories)
    fig_width = min(24, max(16, n_cell_types))  # Scale width by number of cell types
    fig_height = min(24, max(16, len(genes_to_show) * 0.25))  # Scale height by number of genes
    
    # Generate the matrixplot with direct figure access
    fig = plt.figure(figsize=(fig_width, fig_height))
    sc.pl.matrixplot(
        adata, 
        genes_to_show, 
        groupby=plot_groupby,
        dendrogram=True, 
        standard_scale='var', 
        cmap='viridis',
        colorbar_title='Scaled\nexpression',
        ax=fig.gca(),
        show=False
    )
    
    # Clean up temporary column
    if temp_col in adata.obs.columns:
        del adata.obs[temp_col]
    
    # Add better axis labels
    plt.gcf().text(0.5, 0.02, 'Cell Types', ha='center', fontsize=14)
    plt.gcf().text(0.02, 0.5, 'Marker Genes', va='center', rotation='vertical', fontsize=14)
    
    # Adjust gene labels to prevent overlap
    plt.tight_layout(pad=3.0)  # Increase padding
    
    # Get the axes and adjust y-axis label properties
    for ax in plt.gcf().get_axes():
        for tick in ax.get_yticklabels():
            tick.set_fontsize(8)  # Smaller font for y-axis labels
    
    # Save the figure manually to the correct path
    marker_heatmap_file = os.path.join(VIS_DIR, "marker_genes_heatmap.pdf")
    os.makedirs(VIS_DIR, exist_ok=True)  # Ensure directory exists
    plt.savefig(marker_heatmap_file, bbox_inches='tight', dpi=150)
    plt.close()
    print(f"Saved marker heatmap to {marker_heatmap_file}")
    
    # Create a simplified version without dendrogram for easier interpretation
    print("Creating simplified expression heatmap...")
    
    # Create a temporary column with formatted cell type names for the simplified plot
    temp_col = f"{cell_type_column}_display"
    if cell_type_column in adata.obs:
        # Ensure the column is categorical before accessing .cat - using isinstance to avoid deprecation warnings
        if not isinstance(adata.obs[cell_type_column].dtype, pd.CategoricalDtype):
            adata.obs[cell_type_column] = adata.obs[cell_type_column].astype('category')
            
        # Format category names by replacing underscores with spaces and capitalizing properly
        categories = adata.obs[cell_type_column].cat.categories
        category_map = {cat: cat.replace('_', ' ').replace('-', ' ').title() for cat in categories}
        
        # Create formatted column
        adata.obs[temp_col] = adata.obs[cell_type_column].map(category_map)
        adata.obs[temp_col] = adata.obs[temp_col].astype('category')
        adata.obs[temp_col] = adata.obs[temp_col].cat.reorder_categories([category_map[cat] for cat in categories])
        simple_plot_groupby = temp_col
    else:
        simple_plot_groupby = cell_type_column
    
    # Adjust figure size based on content - be more generous with space
    n_genes_simple = min(20, len(genes_to_show))
    genes_simple = genes_to_show[:n_genes_simple]  # Use more genes for better clarity
    n_cell_types = len(adata.obs[cell_type_column].cat.categories)
    
    # Create figure with calculated dimensions
    simple_fig_width = min(22, max(12, n_cell_types * 1.2))
    simple_fig_height = min(20, max(10, n_genes_simple * 0.5))
    plt.figure(figsize=(simple_fig_width, simple_fig_height))
    
    # Plot with specific ax to avoid tight_layout issues
    ax = plt.gca()
    sc.pl.matrixplot(
        adata, 
        genes_simple,
        groupby=simple_plot_groupby,
        standard_scale='var', 
        cmap='viridis',
        colorbar_title='Scaled\nexpression',
        title=f"Gene Expression by {cell_type_display}",
        ax=ax,
        show=False
    )
    
    # Add extra space around the plot instead of using tight_layout
    fig = plt.gcf()
    fig.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.1)
    
    # Clean up temporary column
    if temp_col in adata.obs.columns:
        del adata.obs[temp_col]
        
    simple_heatmap_file = os.path.join(VIS_DIR, "simple_expression_heatmap.pdf")
    plt.savefig(simple_heatmap_file, bbox_inches='tight', dpi=150)
    plt.close()
    print(f"Saved simplified heatmap to {simple_heatmap_file}")
    
    print(f"All visualizations saved to {VIS_DIR}")
    
    return adata
