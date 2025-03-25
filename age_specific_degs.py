#!/usr/bin/env python3
"""
Cell type-specific differential expression analysis across age groups.
Takes the annotated data and identifies differentially expressed genes for each cell type,
with the option to filter by tissue type (tumor or normal).
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

def classify_tissue_type(tissue_origin):
    """Classify tissue origins into tumor or normal categories"""
    if tissue_origin.startswith('n'):
        return 'Normal'
    elif tissue_origin.startswith('t') or tissue_origin.startswith('m') or tissue_origin == 'PE':
        return 'Tumor'
    else:
        return 'Unknown'

def run_celltype_de_analysis(adata, output_dir, min_cells=10, tissue_type=None):
    """
    Run differential expression analysis for each cell type across age groups.
    Uses young/mid-adult (40-54) as reference against each other age group.
    
    Parameters:
    -----------
    tissue_type: str, optional
        Filter analysis to only 'Tumor' or 'Normal' tissue. If None, analyze all cells.
    """
    os.makedirs(output_dir, exist_ok=True)
    results = {}
    
    # Add tissue type classification if not already present
    if 'Tissue_Type' not in adata.obs.columns:
        adata.obs['Tissue_Type'] = adata.obs['Tissue origins'].apply(classify_tissue_type)
        print(f"Tissue type distribution:\n{adata.obs['Tissue_Type'].value_counts()}")
    
    # Filter by tissue type if specified
    if tissue_type is not None:
        print(f"Filtering to only {tissue_type} tissue")
        adata = adata[adata.obs['Tissue_Type'] == tissue_type].copy()
        print(f"After filtering: {adata.n_obs} cells")
    
    # Get unique cell types
    cell_types = adata.obs['predicted_cell_type'].unique()
    reference_group = '40-54 (Mid-Adult)'
    
    # Initialize DataFrame to track number of DEGs per cell type and age group
    summary_df = pd.DataFrame(index=cell_types, 
                             columns=[g for g in adata.obs['Age_Group'].unique() if g != reference_group])
    
    for ct in cell_types:
        print(f"\nAnalyzing cell type: {ct}")
        # Subset to only this cell type
        ct_adata = adata[adata.obs['predicted_cell_type'] == ct].copy()
        
        # Skip if too few cells
        if ct_adata.n_obs < min_cells:
            print(f"  Skipping {ct} - too few cells ({ct_adata.n_obs})")
            continue
            
        # Create directory for cell type
        ct_dir = os.path.join(output_dir, ct.replace(" ", "_").replace("/", "_"))
        os.makedirs(ct_dir, exist_ok=True)
        
        # Perform DE analysis against reference for each age group
        results[ct] = {}
        
        for age_group in ct_adata.obs['Age_Group'].unique():
            if age_group == reference_group:
                continue
                
            # Check if we have enough cells in both groups
            n_ref = np.sum(ct_adata.obs['Age_Group'] == reference_group)
            n_age = np.sum(ct_adata.obs['Age_Group'] == age_group)
            
            if n_ref < min_cells or n_age < min_cells:
                print(f"  Skipping {age_group} - insufficient cells (ref: {n_ref}, target: {n_age})")
                summary_df.loc[ct, age_group] = 0
                continue
                
            print(f"  Comparing {age_group} vs {reference_group}")
            
            # Create grouping variable
            ct_adata.obs['compare_groups'] = 'Other'
            ct_adata.obs.loc[ct_adata.obs['Age_Group'] == reference_group, 'compare_groups'] = reference_group
            ct_adata.obs.loc[ct_adata.obs['Age_Group'] == age_group, 'compare_groups'] = age_group
            
            # Only keep cells in these two groups
            subset = ct_adata[ct_adata.obs['compare_groups'].isin([reference_group, age_group])].copy()
            
            # Run differential expression analysis
            sc.tl.rank_genes_groups(subset, 'compare_groups', reference=reference_group, method='wilcoxon')
            
            # Get results
            de_results = sc.get.rank_genes_groups_df(subset, group=age_group)
            
            # Filter for significance and fold change
            significant_degs = de_results[
                (de_results['pvals_adj'] < 0.05) & 
                (abs(de_results['logfoldchanges']) > 0.5)
            ]
            
            # Save to file
            de_results.to_csv(f"{ct_dir}/{age_group.split()[0]}_vs_reference_all.csv", index=False)
            significant_degs.to_csv(f"{ct_dir}/{age_group.split()[0]}_vs_reference_significant.csv", index=False)
            
            # Track number of DEGs
            n_degs = len(significant_degs)
            summary_df.loc[ct, age_group] = n_degs
            
            # Plot top DEGs
            if n_degs > 0:
                top_genes = list(significant_degs.sort_values('logfoldchanges', ascending=False).head(10)['names'])
                top_genes += list(significant_degs.sort_values('logfoldchanges', ascending=True).head(10)['names'])
                sc.pl.violin(subset, top_genes[:10], groupby='compare_groups', 
                            rotation=90, save=f"_{ct}_{age_group.split()[0]}_up.pdf")
                sc.pl.violin(subset, top_genes[-10:], groupby='compare_groups', 
                            rotation=90, save=f"_{ct}_{age_group.split()[0]}_down.pdf")
                
                # Move files to the right directory
                prefix = f"violin_{ct}_{age_group.split()[0]}"
                for suffix in ["_up.pdf", "_down.pdf"]:
                    if os.path.exists(f"figures/{prefix}{suffix}"):
                        os.rename(f"figures/{prefix}{suffix}", f"{ct_dir}/{prefix}{suffix}")
                
                # Create volcano plot
                plt.figure(figsize=(10, 8))
                plt.scatter(de_results['logfoldchanges'], -np.log10(de_results['pvals']), 
                         alpha=0.5, s=20, color='gray')
                
                # Highlight significant genes
                sig_up = significant_degs[significant_degs['logfoldchanges'] > 0]
                sig_down = significant_degs[significant_degs['logfoldchanges'] < 0]
                
                plt.scatter(sig_up['logfoldchanges'], -np.log10(sig_up['pvals']), 
                         alpha=0.8, s=25, color='red', label=f'Up in {age_group} ({len(sig_up)})')
                plt.scatter(sig_down['logfoldchanges'], -np.log10(sig_down['pvals']), 
                         alpha=0.8, s=25, color='blue', label=f'Up in {reference_group} ({len(sig_down)})')
                
                plt.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
                plt.axvline(0.5, color='gray', linestyle='--', alpha=0.5)
                plt.axvline(-0.5, color='gray', linestyle='--', alpha=0.5)
                
                plt.xlabel('Log2 Fold Change')
                plt.ylabel('-log10(p-value)')
                plt.title(f'{ct}: {age_group} vs {reference_group}')
                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig(f"{ct_dir}/volcano_{age_group.split()[0]}.pdf")
                plt.close()
    
    # Save summary
    summary_df.to_csv(f"{output_dir}/degs_per_celltype_summary.csv")
    
    # Create a heatmap of DEG counts
    plt.figure(figsize=(12, 10))
    summary_df = summary_df.fillna(0).astype(int)
    summary_df = summary_df.loc[summary_df.sum(axis=1) > 0]
    
    if len(summary_df) > 0:
        import seaborn as sns
        sns.heatmap(summary_df, annot=True, fmt="d", cmap="YlOrRd")
        plt.title("Number of Differentially Expressed Genes by Cell Type and Age Group")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/degs_per_celltype_heatmap.pdf")
    
    return summary_df

def main():
    parser = argparse.ArgumentParser(description='Cell type-specific differential expression analysis')
    parser.add_argument('-i', '--input', default="./age_analysis/age_celltype_annotation.h5ad",
                       help="Path to annotated AnnData file (.h5ad)")
    parser.add_argument('-o', '--output', default="./age_analysis/diff_expression",
                       help="Output directory for results")
    parser.add_argument('-m', '--min_cells', type=int, default=20,
                       help="Minimum cells required for each group")
    parser.add_argument('-t', '--tissue_type', choices=['Tumor', 'Normal'],
                       help="Filter analysis to only tumor or normal tissue")
    args = parser.parse_args()

    # Load the annotated data
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells with {adata.n_vars} genes")
    
    # Run differential expression analysis
    # Adjust output directory based on tissue type
    if args.tissue_type:
        output_dir = os.path.join(args.output, args.tissue_type.lower())
    else:
        output_dir = args.output
        
    summary = run_celltype_de_analysis(adata, output_dir, args.min_cells, args.tissue_type)
    
    print("\nAnalysis complete!")
    if len(summary) > 0:
        print("\nSummary of differential expression analysis:")
        print(summary)
    else:
        print("No significant differentially expressed genes found.")

if __name__ == "__main__":
    main()
