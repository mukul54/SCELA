#!/usr/bin/env python3
import scanpy as sc
import decoupler as dc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import gc
import warnings
warnings.filterwarnings('ignore')

# Configuration
sc.settings.set_figure_params(dpi=300, frameon=False, fontsize=10)
sc.settings.verbosity = 1  # Reduce output verbosity

def create_age_groups(adata):
    """Create age groups from continuous age data"""
    print("Creating age groups...")
    bins = [40, 55, 65, 75, 85]
    labels = [
        '40-54 (Mid-Adult)', 
        '55-64 (Pre-Senior)',
        '65-74 (Young Elderly)',
        '75-84 (Senior)'
    ]
    
    # Convert age to numeric, handling any non-numeric values
    adata.obs['Age'] = pd.to_numeric(adata.obs['Age'], errors='coerce')
    
    # Create age groups as categorical
    adata.obs['Age_Group'] = pd.cut(
        adata.obs['Age'],
        bins=bins,
        labels=labels,
        right=False
    )
    
    # Add categories for ages outside the bins
    adata.obs['Age_Group'] = adata.obs['Age_Group'].cat.add_categories(['<40 (Young)', '85+ (Advanced)', 'Unknown'])
    
    # Assign categories to cells outside the bins
    adata.obs.loc[adata.obs['Age'] < 40, 'Age_Group'] = '<40 (Young)'
    adata.obs.loc[adata.obs['Age'] >= 85, 'Age_Group'] = '85+ (Advanced)'
    adata.obs.loc[adata.obs['Age'].isna(), 'Age_Group'] = 'Unknown'
    
    print(f"Age group distribution:\n{adata.obs['Age_Group'].value_counts()}")
    return adata

def preprocess_data(adata, subsample_fraction=0.1, n_top_genes=2000):
    """Memory-efficient preprocessing for large datasets"""
    print("Preprocessing data with memory optimization...")
    
    # Basic filtering - maintain sparsity
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Calculate QC metrics 
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)
    
    # Filter cells by QC metrics
    adata = adata[adata.obs.pct_counts_mt < 20, :].copy()
    
    # Subsample for memory efficiency
    if subsample_fraction < 1.0:
        print(f"Original dataset: {adata.n_obs} cells")
        sc.pp.subsample(adata, fraction=subsample_fraction, random_state=42)
        print(f"Subsampled to {adata.n_obs} cells")
    
    # Memory-efficient normalization
    print("Normalizing data...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Free memory
    gc.collect()
    
    # Select highly variable genes
    print("Finding highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='seurat', batch_key=None)
    print(f"Selected {adata.var.highly_variable.sum()} highly variable genes")
    
    # Use only HVG for dimensionality reduction
    print("Performing dimensionality reduction...")
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    
    # Compute PCA (more efficient than scaling the full matrix)
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=30, svd_solver='arpack')
    
    # Compute neighborhood graph, UMAP, and clustering
    sc.pp.neighbors(adata_hvg, n_neighbors=10, n_pcs=30)
    sc.tl.umap(adata_hvg)
    sc.tl.leiden(adata_hvg, resolution=0.8)
    
    # Transfer results back to main object
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
    adata.obs['leiden'] = adata_hvg.obs['leiden']
    
    # Free memory
    del adata_hvg
    gc.collect()
    
    return adata

def get_cell_markers(specific_cell_types=None):
    """Get cell markers from PanglaoDB"""
    print("Retrieving cell markers...")
    # Query PanglaoDB through OmniPath
    markers = dc.get_resource(
        name='PanglaoDB',
        organism='human',
        license='academic'
    )
    
    # Filter by canonical_marker and human
    markers = markers[
        markers['human'].astype(bool) & 
        markers['canonical_marker'].astype(bool) & 
        (markers['human_sensitivity'].astype(float) > 0.5)
    ]
    
    # Remove duplicated entries
    markers = markers[~markers.duplicated(['cell_type', 'genesymbol'])]
    
    # Filter for specific cell types if provided
    if specific_cell_types:
        markers = markers[markers['cell_type'].isin(specific_cell_types)]
    
    print(f"Retrieved {len(markers)} marker genes for {markers['cell_type'].nunique()} cell types")
    return markers

def annotate_cell_types(adata, markers):
    """Run ORA to identify enriched cell types"""
    print("Running cell type enrichment analysis...")
    dc.run_ora(
        mat=adata,
        net=markers,
        source='cell_type',
        target='genesymbol',
        min_n=3,
        verbose=True,
        use_raw=False
    )
    
    # Extract scores for visualization
    acts = dc.get_acts(adata, obsm_key='ora_estimate')
    
    # Fix infinite values
    acts_v = acts.X.ravel()
    max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
    acts.X[~np.isfinite(acts.X)] = max_e
    
    # Only add top cell types to original object to save memory
    top_celltypes = []
    for cluster in adata.obs['leiden'].unique():
        cells_in_cluster = acts[adata.obs['leiden'] == cluster]
        mean_scores = cells_in_cluster.X.mean(axis=0)
        top_idx = np.argsort(mean_scores)[-5:]  # Top 5 cell types
        top_celltypes.extend(acts.var_names[top_idx])
    
    # Add only top cell types to observations
    top_celltypes = list(set(top_celltypes))
    for col in top_celltypes:
        # Fix: Flatten the array before assignment
        adata.obs[f'ora_{col}'] = acts[:, col].X.flatten()
    
    return adata, acts

def predict_cluster_cell_types(adata, acts):
    """Identify top cell types per cluster"""
    print("Predicting cell types for clusters...")
    try:
        # Rank cell types by clusters
        df = dc.rank_sources_groups(
            acts, 
            groupby='leiden', 
            reference='rest', 
            method='t-test_overestim_var'
        )
        
        # Get top cell type per cluster
        annotation_dict = df.groupby('group').head(1).set_index('group')['names'].to_dict()
        
        # Assign cell types to clusters
        adata.obs['predicted_cell_type'] = [annotation_dict.get(clust, 'Unknown') 
                                           for clust in adata.obs['leiden']]
        
        # Get top 3 cell types per cluster for visualization
        n_ctypes = 3
        ctypes_dict = df.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
        
        return adata, ctypes_dict
    except Exception as e:
        print(f"Error in cell type prediction: {e}")
        # Create dummy annotation
        adata.obs['predicted_cell_type'] = 'Unknown'
        ctypes_dict = {clust: ['Unknown'] for clust in adata.obs['leiden'].unique()}
        return adata, ctypes_dict

def visualize_results(adata, acts, ctypes_dict, output_dir):
    """Generate visualizations for age and cell type analysis"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Create figure directory
    os.makedirs("figures", exist_ok=True)
    
    # 1. Age distribution
    print("Generating age distribution plot...")
    plt.figure(figsize=(10, 6))
    sns.histplot(data=adata.obs, x='Age', hue='Age_Group', bins=range(30, 90, 5))
    plt.title('Age Distribution with Clinical Grouping')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/age_distribution.pdf")
    plt.close()
    
    # 2. UMAP by age group
    print("Generating UMAP by age group...")
    sc.pl.umap(adata, color='Age_Group', title='Age Group', save=f"_age_groups.pdf")
    if os.path.exists("figures/umap_age_groups.pdf"):
        os.rename("figures/umap_age_groups.pdf", f"{output_dir}/umap_age_groups.pdf")
    
    # 3. UMAP by predicted cell type - improved labels
    print("Generating UMAP by cell type...")
    # Temporarily rename for better visualization
    adata.obs['Predicted Cell Type'] = adata.obs['predicted_cell_type']
    sc.pl.umap(adata, color='Predicted Cell Type', save="_cell_types.pdf")
    if os.path.exists("figures/umap_cell_types.pdf"):
        os.rename("figures/umap_cell_types.pdf", f"{output_dir}/umap_cell_types.pdf")
    
    # 4. UMAP by leiden clusters
    print("Generating UMAP by clusters...")
    sc.pl.umap(adata, color='leiden', title='Clusters', legend_loc='on data', save="_clusters.pdf")
    if os.path.exists("figures/umap_clusters.pdf"):
        os.rename("figures/umap_clusters.pdf", f"{output_dir}/umap_clusters.pdf")
    
    # 5. Top cell types per cluster heatmap
    try:
        print("Generating cell type heatmap...")
        # Calculate dendrogram explicitly to avoid warning
        print("Calculating cluster relationships for dendrogram...")
        sc.tl.dendrogram(acts, groupby='leiden')
        
        sc.pl.matrixplot(
            acts, 
            ctypes_dict, 
            'leiden', 
            title='Top Cell Types by Cluster',
            dendrogram=True, 
            standard_scale='var',
            colorbar_title='Z-scaled scores', 
            cmap='RdBu_r',
            save="_cluster_celltypes.pdf"
        )
        if os.path.exists("figures/matrixplot_cluster_celltypes.pdf"):
            os.rename("figures/matrixplot_cluster_celltypes.pdf", f"{output_dir}/cluster_celltypes.pdf")
    except Exception as e:
        print(f"Error generating heatmap: {e}")
    
    # 6. Cell type by age group composition - IMPROVED SIZE AND LABELS
    try:
        print("Generating cell type by age group composition...")
        plt.figure(figsize=(18, 10))  # Much larger figure
        
        # Create the crosstab and clean up the display
        ct_data = pd.crosstab(
            adata.obs['predicted_cell_type'], 
            adata.obs['Age_Group'], 
            normalize='index'
        )
        
        # Plot with better formatting
        ax = ct_data.plot.bar(stacked=True, figsize=(18, 10))
        plt.title('Cell Type Composition by Age Group', fontsize=16)
        plt.xlabel('Predicted Cell Type', fontsize=14)
        plt.ylabel('Proportion', fontsize=14)
        plt.xticks(rotation=45, ha='right', fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(title='Age Group', fontsize=12, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/celltype_by_age.pdf", bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"Error generating composition plot: {e}")
    
    print(f"Results saved to {output_dir}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Age-stratified cell type annotation')
    parser.add_argument('-i', '--input', required=True,
                      help="Path to AnnData file (.h5ad)")
    parser.add_argument('-o', '--output', required=True,
                      help="Output directory for results")
    parser.add_argument('--subsample', type=float, default=0.1,
                      help="Fraction of cells to randomly subsample (0-1)")
    parser.add_argument('--hvg', type=int, default=2000,
                      help="Number of highly variable genes to use")
    args = parser.parse_args()
    
    print(f"Loading data from {args.input}...")
    try:
        adata = sc.read_h5ad(args.input)
        print(f"Loaded data with {adata.n_obs} cells and {adata.n_vars} genes")
    except Exception as e:
        print(f"Error loading data: {e}")
        return
    
    # Create age groups
    adata = create_age_groups(adata)
    
    # Preprocess data
    adata = preprocess_data(adata, subsample_fraction=args.subsample, n_top_genes=args.hvg)
    
    # Get cell markers from PanglaoDB
    markers = get_cell_markers()
    
    # Run cell type annotation
    adata, acts = annotate_cell_types(adata, markers)
    
    # Predict cluster cell types
    adata, ctypes_dict = predict_cluster_cell_types(adata, acts)
    
    # Visualize results
    visualize_results(adata, acts, ctypes_dict, args.output)
    
    # Save results
    print(f"Saving results to {args.output}/age_celltype_annotation.h5ad")
    adata.write(f"{args.output}/age_celltype_annotation.h5ad")
    print("Analysis complete!")

if __name__ == "__main__":
    main()