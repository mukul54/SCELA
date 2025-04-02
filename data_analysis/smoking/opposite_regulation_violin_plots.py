#!/usr/bin/env python3
"""
Generate violin plots for genes with opposite regulation patterns across smoking statuses.
This script loads the genes identified as having opposite regulation patterns between
different smoking groups and creates violin plots to visualize their expression levels.
"""

import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Define paths
base_dir = "/home/mukul.ranjan/Documents/cancer/SCELA"
output_dir = os.path.join(base_dir, "smoking_tumor_analysis_results")
os.makedirs(os.path.join(output_dir, "violin_plots"), exist_ok=True)

# Input h5ad file - adjust path if necessary
adata_file = "/home/mukul.ranjan/Documents/cancer/SCELA/age_analysis/age_celltype_annotation.h5ad"

# Load opposite regulation gene data
never_ex_file = os.path.join(output_dir, "never_vs_ex_opposite_regulation.csv")
never_cur_file = os.path.join(output_dir, "never_vs_cur_opposite_regulation.csv")
ex_cur_file = os.path.join(output_dir, "ex_vs_cur_opposite_regulation.csv")

never_ex_data = pd.read_csv(never_ex_file)
never_cur_data = pd.read_csv(never_cur_file)
ex_cur_data = pd.read_csv(ex_cur_file)

print(f"Loaded opposite regulation genes:")
print(f"Never vs Ex: {len(never_ex_data)} genes")
print(f"Never vs Current: {len(never_cur_data)} genes")
print(f"Ex vs Current: {len(ex_cur_data)} genes")

def get_top_genes(df, n=5):
    """Get top N genes by magnitude of fold change difference."""
    df = df.copy()
    
    # Calculate magnitude of difference for appropriate columns
    if 'Never_log2FC' in df.columns and 'Ex_log2FC' in df.columns:
        df['magnitude'] = abs(df['Never_log2FC'] - df['Ex_log2FC'])
    elif 'Never_log2FC' in df.columns and 'Cur_log2FC' in df.columns:
        df['magnitude'] = abs(df['Never_log2FC'] - df['Cur_log2FC'])
    elif 'Ex_log2FC' in df.columns and 'Cur_log2FC' in df.columns:
        df['magnitude'] = abs(df['Ex_log2FC'] - df['Cur_log2FC'])
    
    # Return top genes by magnitude
    return df.sort_values('magnitude', ascending=False).head(n)['Gene'].tolist()

# Get top 5 genes from each comparison
top_never_ex = get_top_genes(never_ex_data, 5)
top_never_cur = get_top_genes(never_cur_data, 5)
top_ex_cur = get_top_genes(ex_cur_data, 5)

print(f"\nTop genes with opposite regulation:")
print(f"Never vs Ex: {', '.join(top_never_ex)}")
print(f"Never vs Current: {', '.join(top_never_cur)}")
print(f"Ex vs Current: {', '.join(top_ex_cur)}")

# Load annotated data
print(f"\nLoading annotated data from {adata_file}")
adata = sc.read_h5ad(adata_file)
print(f"Loaded {adata.n_obs} cells with {adata.n_vars} genes")

# Function to create violin plots
def create_violin_plots_for_genes(genes, title, output_file, cell_type=None):
    """
    Create violin plots for specific genes, comparing expression across smoking categories.
    
    Parameters:
    -----------
    genes : list
        List of gene names to plot
    title : str
        Plot title
    output_file : str
        Output file path
    cell_type : str, optional
        If provided, filter to this cell type
    """
    # Create a copy of the dataset to avoid modifying the original
    plot_data = adata.copy()
    
    # Filter to specific cell type if requested
    if cell_type:
        print(f"Filtering to cell type: {cell_type}")
        plot_data = plot_data[plot_data.obs['predicted_cell_type'] == cell_type].copy()
        if plot_data.n_obs < 20:
            print(f"  Too few cells ({plot_data.n_obs}) for cell type {cell_type}")
            return False
    
    # Create a figure with multiple subplots for each gene
    if len(genes) > 0:
        # Calculate number of rows and columns for subplots
        n_rows = min(3, len(genes))
        n_cols = int(np.ceil(len(genes) / n_rows))
        fig_height = 4 * n_rows
        fig_width = 4 * n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))
        axes = axes.flatten() if isinstance(axes, np.ndarray) else [axes]
        
        for i, gene in enumerate(genes):
            if i < len(axes):
                ax = axes[i]
                
                # Skip if gene not in dataset
                if gene not in plot_data.var_names:
                    print(f"  Gene {gene} not found in dataset")
                    ax.text(0.5, 0.5, f"Gene {gene} not found", 
                          ha='center', va='center', fontsize=12)
                    ax.axis('off')
                    continue
                
                # Get expression data for this gene
                gene_data = []
                smoking_categories = []
                
                for smoking in sorted(plot_data.obs['Smoking'].unique()):
                    cells = plot_data[plot_data.obs['Smoking'] == smoking]
                    if cells.n_obs > 0:
                        expr = cells[:, gene].X.toarray().flatten()
                        gene_data.extend(expr)
                        smoking_categories.extend([smoking] * len(expr))
                
                # Create DataFrame for plotting
                plot_df = pd.DataFrame({
                    'Expression': gene_data,
                    'Smoking': smoking_categories
                })
                
                # Create violin plot
                sns.violinplot(x='Smoking', y='Expression', data=plot_df, ax=ax, 
                              palette={"Never": "blue", "Ex": "orange", "Cur": "red"}, 
                              order=["Never", "Ex", "Cur"])
                
                # Add individual points
                sns.stripplot(x='Smoking', y='Expression', data=plot_df, ax=ax, 
                             color='black', alpha=0.3, size=2, jitter=True,
                             order=["Never", "Ex", "Cur"])
                
                ax.set_title(gene, fontsize=12, fontweight='bold')
                ax.set_ylabel('Expression' if i % n_cols == 0 else '')
                ax.set_xlabel('')
                
                # Add statistical annotation if possible
                try:
                    # Get average expression for each smoking category
                    means = plot_df.groupby('Smoking')['Expression'].mean()
                    
                    # Add text annotations for the means
                    for j, smoking in enumerate(["Never", "Ex", "Cur"]):
                        if smoking in means:
                            ax.text(j, means[smoking] + 0.05 * ax.get_ylim()[1],
                                  f"{means[smoking]:.2f}", ha='center', fontsize=9)
                except:
                    pass  # Skip if any issues
        
        # Turn off unused subplots
        for i in range(len(genes), len(axes)):
            axes[i].axis('off')
        
        plt.suptitle(title, fontsize=16, fontweight='bold')
        plt.tight_layout(rect=[0, 0, 1, 0.96])  # Make room for suptitle
        
        # Add cell type to filename if specified
        if cell_type:
            output_file = output_file.replace('.pdf', f'_{cell_type.replace(" ", "_")}.pdf')
        
        plt.savefig(output_file)
        plt.close()
        
        print(f"  Saved violin plots to {output_file}")
        return True
    
    return False

# Create violin plots for all cells
print("\nGenerating violin plots for all cells:")

# For Never vs Ex
print("  Creating violin plots for Never vs Ex genes")
create_violin_plots_for_genes(
    top_never_ex,
    "Top Genes with Opposite Regulation: Never vs Ex Smokers",
    os.path.join(output_dir, "violin_plots", "never_ex_opposite_regulation_violin.pdf")
)

# For Never vs Current
print("  Creating violin plots for Never vs Current genes")
create_violin_plots_for_genes(
    top_never_cur,
    "Top Genes with Opposite Regulation: Never vs Current Smokers",
    os.path.join(output_dir, "violin_plots", "never_cur_opposite_regulation_violin.pdf")
)

# For Ex vs Current
print("  Creating violin plots for Ex vs Current genes")
create_violin_plots_for_genes(
    top_ex_cur,
    "Top Genes with Opposite Regulation: Ex vs Current Smokers",
    os.path.join(output_dir, "violin_plots", "ex_cur_opposite_regulation_violin.pdf")
)

# Get list of top cell types with the most DEGs
def get_top_celltypes():
    """Get cell types with the most smoking-dependent DEGs."""
    summary_file = os.path.join(base_dir, "smoking_tumor_interaction/cell_types/smoking_tumor_celltype_degs_summary.csv")
    if os.path.exists(summary_file):
        summary_df = pd.read_csv(summary_file, index_col=0)
        # Sum all DEGs across smoking statuses
        total_degs = summary_df.sum(axis=1)
        return total_degs.sort_values(ascending=False).head(5).index.tolist()
    return []

# Create violin plots for top cell types
top_celltypes = get_top_celltypes()
if top_celltypes:
    print(f"\nGenerating cell type-specific violin plots for top cell types:")
    print(f"Top cell types: {', '.join(top_celltypes)}")
    
    for cell_type in top_celltypes:
        print(f"\nProcessing cell type: {cell_type}")
        
        # For Never vs Ex
        print(f"  Creating Never vs Ex violin plots for {cell_type}")
        create_violin_plots_for_genes(
            top_never_ex,
            f"Never vs Ex: {cell_type}",
            os.path.join(output_dir, "violin_plots", "never_ex_celltype_violin.pdf"),
            cell_type=cell_type
        )
        
        # For Never vs Current
        print(f"  Creating Never vs Current violin plots for {cell_type}")
        create_violin_plots_for_genes(
            top_never_cur,
            f"Never vs Current: {cell_type}",
            os.path.join(output_dir, "violin_plots", "never_cur_celltype_violin.pdf"),
            cell_type=cell_type
        )
        
        # For Ex vs Current
        print(f"  Creating Ex vs Current violin plots for {cell_type}")
        create_violin_plots_for_genes(
            top_ex_cur,
            f"Ex vs Current: {cell_type}",
            os.path.join(output_dir, "violin_plots", "ex_cur_celltype_violin.pdf"),
            cell_type=cell_type
        )

# Also create combined plots for immune-related cell types
immune_celltypes = [ct for ct in adata.obs['predicted_cell_type'].unique() 
                   if any(x in ct.lower() for x in ['t cell', 'b cell', 'macrophage', 'dendritic', 'nk'])]

if immune_celltypes:
    print(f"\nGenerating violin plots for immune cell types:")
    
    # Combine top genes from all comparisons for immune cells
    all_top_genes = list(set(top_never_ex + top_never_cur + top_ex_cur))
    
    for cell_type in immune_celltypes:
        print(f"  Processing immune cell type: {cell_type}")
        create_violin_plots_for_genes(
            all_top_genes[:6],  # Limit to 6 to keep plots manageable
            f"Opposite Regulation Genes in {cell_type}",
            os.path.join(output_dir, "violin_plots", "immune_opposite_regulation_violin.pdf"),
            cell_type=cell_type
        )

print("\nViolin plot generation complete.")
