#!/usr/bin/env python3
"""
Further analysis of age-stratified cell type annotations.
This script loads the previously processed data and enables additional analysis.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import decoupler as dc

# Set plotting parameters
sc.settings.set_figure_params(dpi=300, frameon=False, fontsize=12)
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['figure.dpi'] = 300

def load_analysis_data(h5ad_path):
    """Load previously processed data"""
    print(f"Loading annotated data from {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"Loaded data with {adata.n_obs} cells and {adata.n_vars} genes")
    return adata

def get_age_group_stats(adata):
    """Analyze cell type composition by age group"""
    print("Analyzing age group statistics...")
    
    # Get cell counts by age group
    age_counts = adata.obs['Age_Group'].value_counts().sort_index()
    
    # Get cell type counts by age group
    ct_by_age = pd.crosstab(
        adata.obs['Age_Group'],
        adata.obs['predicted_cell_type'],
    )
    
    # Calculate percentages
    ct_by_age_pct = ct_by_age.div(ct_by_age.sum(axis=1), axis=0) * 100
    
    return age_counts, ct_by_age, ct_by_age_pct

def plot_improved_visualizations(adata, output_dir, ct_by_age_pct=None):
    """Generate improved visualizations"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Recalculate crosstabs if not provided
    if ct_by_age_pct is None:
        # Get cell type counts by age group
        ct_by_age = pd.crosstab(
            adata.obs['Age_Group'],
            adata.obs['predicted_cell_type'],
        )
        
        # Calculate percentages
        ct_by_age_pct = ct_by_age.div(ct_by_age.sum(axis=1), axis=0) * 100
        
        # Save as CSV
        ct_by_age.to_csv(f"{output_dir}/celltype_counts_by_age.csv")
        ct_by_age_pct.to_csv(f"{output_dir}/celltype_percent_by_age.csv")
    
    # 1. Improved cell type composition by age with better formatting
    plt.figure(figsize=(18, 10))
    
    # Create dataframe for plotting with nicer labels
    plot_df = adata.obs.copy()
    plot_df['Predicted Cell Type'] = plot_df['predicted_cell_type']
    
    # Create the crosstab
    ct_data = pd.crosstab(
        plot_df['Predicted Cell Type'], 
        plot_df['Age_Group'], 
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
    plt.savefig(f"{output_dir}/improved_celltype_by_age.pdf", bbox_inches='tight')
    plt.close()
    
    # 2. Cell type counts across age groups - absolute values
    cell_counts = pd.crosstab(plot_df['Age_Group'], plot_df['Predicted Cell Type'])
    plt.figure(figsize=(18, 10))
    sns.heatmap(cell_counts, annot=False, fmt='d', cmap='viridis')
    plt.title('Cell Type Counts by Age Group', fontsize=16)
    plt.ylabel("Age Group")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/celltype_counts_by_age.pdf", bbox_inches='tight')
    plt.close()
    
    # Ensure we have proper clustering info for visualization
    if 'dendrogram_leiden' not in adata.uns:
        print("Computing neighborhood graph and clustering relationships...")
        # Re-calculate UMAP if needed
        if 'X_umap' not in adata.obsm:
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            sc.tl.leiden(adata)
        # Calculate dendrogram
        sc.tl.dendrogram(adata, groupby='leiden')
    
    # 3. UMAP with improved annotations
    sc.pl.umap(
        adata, 
        color=['Age_Group', 'predicted_cell_type'], 
        title=['Age Group', 'Predicted Cell Type'],
        save="_improved.pdf"
    )
    if os.path.exists("figures/umap_improved.pdf"):
        os.rename("figures/umap_improved.pdf", f"{output_dir}/umap_improved.pdf")
    
    # 4. Age distribution by cell type
    plt.figure(figsize=(16, 10))
    sns.boxplot(data=plot_df, x='Predicted Cell Type', y='Age')
    plt.title('Age Distribution by Cell Type', fontsize=16)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/age_by_celltype.pdf", bbox_inches='tight')
    plt.close()
    
    # 5. Cell type enrichment by age group
    # Calculate fold change of cell types across age groups
    baseline = ct_by_age_pct.loc['40-54 (Mid-Adult)']  # Use mid-adult as baseline
    
    for age in ct_by_age_pct.index:
        if age != '40-54 (Mid-Adult)':
            # Calculate fold change compared to baseline
            fold_changes = ct_by_age_pct.loc[age] / baseline
            fold_changes = fold_changes.sort_values(ascending=False)
            
            # Plot top enriched cell types
            plt.figure(figsize=(14, 8))
            fold_changes.head(10).plot.bar()
            plt.axhline(y=1.0, color='r', linestyle='-', alpha=0.3)
            plt.title(f'Cell Types Enriched in {age} vs Mid-Adult', fontsize=16)
            plt.ylabel('Fold Change', fontsize=14)
            plt.xlabel('Cell Type', fontsize=14)
            plt.tight_layout()
            plt.savefig(f"{output_dir}/enrichment_{age.split()[0]}.pdf", bbox_inches='tight')
            plt.close()
    
    print(f"Improved visualizations saved to {output_dir}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Further analysis of age-stratified cell types')
    parser.add_argument('-i', '--input', default="./age_analysis/age_celltype_annotation.h5ad",
                      help="Path to annotated AnnData file (.h5ad)")
    parser.add_argument('-o', '--output', default="./age_analysis/additional_plots",
                      help="Output directory for additional results")
    args = parser.parse_args()
    
    # Load previously processed data
    adata = load_analysis_data(args.input)
    
    # Get age group statistics
    age_counts, ct_by_age, ct_by_age_pct = get_age_group_stats(adata)
    
    # Save statistics to CSV
    print("Saving analysis data to CSV files...")
    age_counts.to_frame().to_csv(f"{args.output}/age_group_counts.csv")
    ct_by_age.to_csv(f"{args.output}/celltype_counts_by_age.csv")
    ct_by_age_pct.to_csv(f"{args.output}/celltype_percent_by_age.csv")
    
    # Save cell metrics to CSV
    adata.obs.to_csv(f"{args.output}/cell_annotations.csv")
    
    # Display summary statistics
    print("\nCell counts by age group:")
    print(age_counts)
    
    print("\nTop 3 cell types in each age group:")
    for age in ct_by_age_pct.index:
        top3 = ct_by_age_pct.loc[age].sort_values(ascending=False).head(3)
        print(f"\n{age}:")
        for ct, pct in top3.items():
            print(f"  {ct}: {pct:.1f}%")
    
    # Generate improved visualizations
    plot_improved_visualizations(adata, args.output, ct_by_age_pct)
    
    print("Analysis complete!")

if __name__ == "__main__":
    main()
