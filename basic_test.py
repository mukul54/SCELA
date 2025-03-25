#!/usr/bin/env python3
import sys
import os
import scanpy as sc
import matplotlib
matplotlib.use('PDF')  # Set backend for PDF output
import matplotlib.pyplot as plt
import pandas as pd

# Configure plot settings for publication quality
sc.settings.set_figure_params(
    dpi=300,
    facecolor='white',
    vector_friendly=True,
    format='pdf',
    frameon=False,
    fontsize=10,
    dpi_save=300
)

def validate_metadata(adata):
    """Ensure required clinical metadata exists in AnnData"""
    required_columns = ['Sex', 'Age', 'Smoking', 'Stages', 'EGFR', 'Histology']
    missing = [col for col in required_columns if col not in adata.obs.columns]
    
    if missing:
        raise ValueError(f"Missing critical metadata columns: {missing}")
    print("All required metadata present:", ", ".join(required_columns))

def plot_umap_by_metadata(adata, metadata_columns, output_dir):
    """Generate individual UMAP plots for each metadata column"""
    os.makedirs(f"{output_dir}/umap", exist_ok=True)
    
    for col in metadata_columns:
        fig = plt.figure(figsize=(6, 4))
        sc.pl.umap(
            adata,
            color=col,
            frameon=False,
            title=f'UMAP - {col}',
            show=False,
            legend_loc='on data' if adata.obs[col].nunique() < 10 else 'right margin'
        )
        plt.tight_layout()
        fig.savefig(f"{output_dir}/umap/UMAP_{col.replace(' ', '_')}.pdf")
        plt.close(fig)

def plot_qc_metrics(adata, output_dir):
    """Generate QC metric visualizations"""
    os.makedirs(f"{output_dir}/qc", exist_ok=True)
    
    # Violin plots for key metrics
    metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    for metric in metrics:
        fig = plt.figure(figsize=(8, 4))
        sc.pl.violin(
            adata,
            keys=metric,
            groupby='Smoking',
            stripplot=False,
            rotation=45,
            show=False
        )
        plt.title(f'{metric} by Smoking Status')
        plt.tight_layout(pad=3)
        fig.savefig(f"{output_dir}/qc/Smoking_{metric}.pdf")
        plt.close(fig)
    
    # Scatter plot of counts vs genes
    fig = plt.figure(figsize=(5, 5))
    sc.pl.scatter(
        adata,
        x='total_counts',
        y='n_genes_by_counts',
        color='pct_counts_mt',
        show=False,
        title='Counts vs Genes'
    )
    fig.savefig(f"{output_dir}/qc/Counts_vs_Genes.pdf")
    plt.close(fig)

def plot_marker_genes(adata, output_dir):
    """Generate marker gene heatmap with clean labels"""
    os.makedirs(f"{output_dir}/markers", exist_ok=True)
    
    # Calculate dendrogram for clean ordering
    sc.tl.dendrogram(adata, groupby='leiden')
    
    # Create heatmap with adjusted parameters
    fig = plt.figure(figsize=(15, 10))  # Increased height for gene labels
    sc.pl.rank_genes_groups_heatmap(
        adata,
        groups=adata.obs['leiden'].cat.categories,
        n_genes=5,
        groupby='leiden',
        show=False,
        show_gene_labels=True,  # Force show all gene labels
        var_group_rotation=45,  # Rotate cluster labels
        dendrogram=True,
        figsize=(15, 10)
    )
    plt.subplots_adjust(bottom=0.3, right=0.8)  # Adjust spacing
    fig.savefig(f"{output_dir}/markers/Marker_Genes_Heatmap.pdf")
    plt.close(fig)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True,
                       help="Path to processed H5AD file")
    parser.add_argument('-o', '--output', required=True,
                       help="Output directory for figures")
    args = parser.parse_args()
    
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    
    # Validate metadata before visualization
    validate_metadata(adata)
    
    # Generate all visualizations
    metadata_cols = ['leiden', 'Sex', 'Age', 'Smoking', 'Stages']
    plot_umap_by_metadata(adata, metadata_cols, args.output)
    plot_qc_metrics(adata, args.output)
    plot_marker_genes(adata, args.output)

if __name__ == "__main__":
    main()