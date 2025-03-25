#!/usr/bin/env python3
"""
Tumor vs. Normal differential expression analysis across age groups.
Identifies genes that are differentially expressed between tumor and normal tissues,
stratified by age group and cell type.
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

def analyze_tumor_vs_normal(adata, output_dir, min_cells=20):
    """Perform differential expression between tumor and normal tissues"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Add a simpler tissue type column (Tumor vs Normal)
    adata.obs['Tissue_Type'] = adata.obs['Tissue origins'].apply(classify_tissue_type)
    print(f"Tissue type distribution:\n{adata.obs['Tissue_Type'].value_counts()}")
    
    # Save tissue type information
    tissue_counts = pd.crosstab(
        adata.obs['Tissue origins'], 
        adata.obs['Tissue_Type']
    )
    tissue_counts.to_csv(f"{output_dir}/tissue_classification.csv")
    
    # Initialize results tracking
    results = {}
    
    # 1. Global tumor vs normal comparison (all age groups)
    print("\nAnalyzing overall tumor vs normal differences...")
    sc.tl.rank_genes_groups(adata, 'Tissue_Type', reference='Normal', method='wilcoxon')
    
    # Save results
    de_results = sc.get.rank_genes_groups_df(adata, group='Tumor')
    significant_degs = de_results[
        (de_results['pvals_adj'] < 0.05) & 
        (abs(de_results['logfoldchanges']) > 0.5)
    ]
    
    de_results.to_csv(f"{output_dir}/global_tumor_vs_normal_all.csv", index=False)
    significant_degs.to_csv(f"{output_dir}/global_tumor_vs_normal_significant.csv", index=False)
    
    # Create global volcano plot
    if len(significant_degs) > 0:
        plt.figure(figsize=(12, 10))
        plt.scatter(de_results['logfoldchanges'], -np.log10(de_results['pvals']), 
                 alpha=0.5, s=20, color='gray')
        
        # Highlight significant genes
        sig_up = significant_degs[significant_degs['logfoldchanges'] > 0]
        sig_down = significant_degs[significant_degs['logfoldchanges'] < 0]
        
        plt.scatter(sig_up['logfoldchanges'], -np.log10(sig_up['pvals']), 
                 alpha=0.8, s=25, color='red', label=f'Up in Tumor ({len(sig_up)})')
        plt.scatter(sig_down['logfoldchanges'], -np.log10(sig_down['pvals']), 
                 alpha=0.8, s=25, color='blue', label=f'Up in Normal ({len(sig_down)})')
        
        # Add labels for top genes
        top_up_genes = sig_up.sort_values('logfoldchanges', ascending=False).head(10)
        top_down_genes = sig_down.sort_values('logfoldchanges', ascending=True).head(10)
        
        for _, gene in top_up_genes.iterrows():
            plt.text(gene['logfoldchanges'], -np.log10(gene['pvals']), 
                    gene['names'], fontsize=9, ha='center', va='bottom')
            
        for _, gene in top_down_genes.iterrows():
            plt.text(gene['logfoldchanges'], -np.log10(gene['pvals']), 
                    gene['names'], fontsize=9, ha='center', va='bottom')
        
        plt.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
        plt.axvline(0.5, color='gray', linestyle='--', alpha=0.5)
        plt.axvline(-0.5, color='gray', linestyle='--', alpha=0.5)
        
        plt.xlabel('Log2 Fold Change', fontsize=12)
        plt.ylabel('-log10(p-value)', fontsize=12)
        plt.title('Global Tumor vs Normal Comparison', fontsize=14)
        plt.legend(loc='best', fontsize=10)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/global_volcano_tumor_vs_normal.pdf")
        plt.close()
        
        # Create violin plots for top DEGs
        top_genes = list(significant_degs.sort_values('logfoldchanges', ascending=False).head(10)['names'])
        sc.pl.violin(adata, top_genes, groupby='Tissue_Type', 
                    rotation=90, save="_global_tumor_vs_normal_up.pdf")
        
        bottom_genes = list(significant_degs.sort_values('logfoldchanges', ascending=True).head(10)['names'])
        sc.pl.violin(adata, bottom_genes, groupby='Tissue_Type', 
                    rotation=90, save="_global_tumor_vs_normal_down.pdf")
        
        # Move violin plots to the right directory
        if os.path.exists("figures/violin_global_tumor_vs_normal_up.pdf"):
            os.rename("figures/violin_global_tumor_vs_normal_up.pdf", 
                     f"{output_dir}/violin_global_tumor_vs_normal_up.pdf")
        
        if os.path.exists("figures/violin_global_tumor_vs_normal_down.pdf"):
            os.rename("figures/violin_global_tumor_vs_normal_down.pdf", 
                     f"{output_dir}/violin_global_tumor_vs_normal_down.pdf")
    
    # 2. Age-specific tumor vs normal comparison
    age_groups = adata.obs['Age_Group'].unique()
    
    # Track DEG counts across age groups
    age_deg_counts = pd.DataFrame(index=['All_Cells'], columns=age_groups)
    
    for age in age_groups:
        print(f"\nAnalyzing tumor vs normal in age group: {age}")
        # Subset to this age group
        age_adata = adata[adata.obs['Age_Group'] == age].copy()
        
        # Check if we have enough cells in both tumor and normal groups
        n_tumor = np.sum(age_adata.obs['Tissue_Type'] == 'Tumor')
        n_normal = np.sum(age_adata.obs['Tissue_Type'] == 'Normal')
        
        if n_tumor < min_cells or n_normal < min_cells:
            print(f"  Skipping {age} - insufficient cells (normal: {n_normal}, tumor: {n_tumor})")
            age_deg_counts.loc['All_Cells', age] = 0
            continue
        
        # Run differential expression
        sc.tl.rank_genes_groups(age_adata, 'Tissue_Type', reference='Normal', method='wilcoxon')
        
        # Get results
        de_results = sc.get.rank_genes_groups_df(age_adata, group='Tumor')
        significant_degs = de_results[
            (de_results['pvals_adj'] < 0.05) & 
            (abs(de_results['logfoldchanges']) > 0.5)
        ]
        
        # Save to file
        age_dir = os.path.join(output_dir, age.replace(" ", "_").replace("/", "_"))
        os.makedirs(age_dir, exist_ok=True)
        
        de_results.to_csv(f"{age_dir}/tumor_vs_normal_all.csv", index=False)
        significant_degs.to_csv(f"{age_dir}/tumor_vs_normal_significant.csv", index=False)
        
        # Track DEG counts
        n_degs = len(significant_degs)
        age_deg_counts.loc['All_Cells', age] = n_degs
        
        # Create volcano plot
        if n_degs > 0:
            plt.figure(figsize=(10, 8))
            plt.scatter(de_results['logfoldchanges'], -np.log10(de_results['pvals']), 
                     alpha=0.5, s=20, color='gray')
            
            # Highlight significant genes
            sig_up = significant_degs[significant_degs['logfoldchanges'] > 0]
            sig_down = significant_degs[significant_degs['logfoldchanges'] < 0]
            
            plt.scatter(sig_up['logfoldchanges'], -np.log10(sig_up['pvals']), 
                     alpha=0.8, s=25, color='red', label=f'Up in Tumor ({len(sig_up)})')
            plt.scatter(sig_down['logfoldchanges'], -np.log10(sig_down['pvals']), 
                     alpha=0.8, s=25, color='blue', label=f'Up in Normal ({len(sig_down)})')
            
            # Add labels for top genes
            top_up_genes = sig_up.sort_values('logfoldchanges', ascending=False).head(10)
            top_down_genes = sig_down.sort_values('logfoldchanges', ascending=True).head(10)
            
            for _, gene in top_up_genes.head(5).iterrows():
                plt.text(gene['logfoldchanges'], -np.log10(gene['pvals']), 
                        gene['names'], fontsize=9, ha='center', va='bottom')
                
            for _, gene in top_down_genes.head(5).iterrows():
                plt.text(gene['logfoldchanges'], -np.log10(gene['pvals']), 
                        gene['names'], fontsize=9, ha='center', va='bottom')
                
            plt.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
            plt.axvline(0.5, color='gray', linestyle='--', alpha=0.5)
            plt.axvline(-0.5, color='gray', linestyle='--', alpha=0.5)
            
            plt.xlabel('Log2 Fold Change')
            plt.ylabel('-log10(p-value)')
            plt.title(f'Tumor vs Normal: {age}')
            plt.legend(loc='best')
            plt.tight_layout()
            plt.savefig(f"{age_dir}/volcano_tumor_vs_normal.pdf")
            plt.close()
            
            # Add violin plots for top DEGs
            age_subset = adata[adata.obs['Age_Group'] == age].copy()
            
            # Top up-regulated genes (higher in tumor)
            top_genes = list(significant_degs.sort_values('logfoldchanges', ascending=False).head(10)['names'])
            if top_genes:
                sc.pl.violin(age_subset, top_genes, groupby='Tissue_Type', 
                            rotation=90, save=f"_{age.replace(' ', '_')}_tumor_up.pdf")
                violin_file = f"figures/violin_{age.replace(' ', '_')}_tumor_up.pdf"
                if os.path.exists(violin_file):
                    os.rename(violin_file, f"{age_dir}/violin_tumor_up.pdf")
                    
            # Top down-regulated genes (higher in normal)
            bottom_genes = list(significant_degs.sort_values('logfoldchanges', ascending=True).head(10)['names'])
            if bottom_genes:
                sc.pl.violin(age_subset, bottom_genes, groupby='Tissue_Type', 
                            rotation=90, save=f"_{age.replace(' ', '_')}_normal_up.pdf")
                violin_file = f"figures/violin_{age.replace(' ', '_')}_normal_up.pdf"
                if os.path.exists(violin_file):
                    os.rename(violin_file, f"{age_dir}/violin_normal_up.pdf")
    
    # 3. Cell type-specific tumor vs normal comparison
    cell_types = adata.obs['predicted_cell_type'].unique()
    
    # Initialize tracking for cell type-specific DEGs
    celltype_deg_counts = pd.DataFrame(index=cell_types, columns=age_groups)
    
    for ct in cell_types:
        print(f"\nAnalyzing cell type: {ct}")
        # Subset to this cell type
        ct_adata = adata[adata.obs['predicted_cell_type'] == ct].copy()
        
        if ct_adata.n_obs < min_cells:
            print(f"  Skipping {ct} - too few cells overall ({ct_adata.n_obs})")
            continue
        
        # Create directory for cell type
        ct_dir = os.path.join(output_dir, ct.replace(" ", "_").replace("/", "_"))
        os.makedirs(ct_dir, exist_ok=True)
        
        # Check tumor vs normal counts for this cell type
        ct_tissue_counts = ct_adata.obs['Tissue_Type'].value_counts()
        pd.DataFrame(ct_tissue_counts).to_csv(f"{ct_dir}/tissue_counts.csv")
        
        if 'Normal' not in ct_tissue_counts or 'Tumor' not in ct_tissue_counts:
            print(f"  Skipping {ct} - missing either tumor or normal cells")
            continue
            
        if ct_tissue_counts['Normal'] < min_cells or ct_tissue_counts['Tumor'] < min_cells:
            print(f"  Skipping {ct} - insufficient cells in tumor or normal")
            continue
        
        # Global analysis for this cell type (all ages)
        sc.tl.rank_genes_groups(ct_adata, 'Tissue_Type', reference='Normal', method='wilcoxon')
        
        # Get results
        de_results = sc.get.rank_genes_groups_df(ct_adata, group='Tumor')
        significant_degs = de_results[
            (de_results['pvals_adj'] < 0.05) & 
            (abs(de_results['logfoldchanges']) > 0.5)
        ]
        
        # Save to file
        de_results.to_csv(f"{ct_dir}/tumor_vs_normal_all.csv", index=False)
        significant_degs.to_csv(f"{ct_dir}/tumor_vs_normal_significant.csv", index=False)
        
        # Create volcano plot
        if len(significant_degs) > 0:
            plt.figure(figsize=(10, 8))
            plt.scatter(de_results['logfoldchanges'], -np.log10(de_results['pvals']), 
                     alpha=0.5, s=20, color='gray')
            
            # Highlight significant genes
            sig_up = significant_degs[significant_degs['logfoldchanges'] > 0]
            sig_down = significant_degs[significant_degs['logfoldchanges'] < 0]
            
            plt.scatter(sig_up['logfoldchanges'], -np.log10(sig_up['pvals']), 
                     alpha=0.8, s=25, color='red', label=f'Up in Tumor ({len(sig_up)})')
            plt.scatter(sig_down['logfoldchanges'], -np.log10(sig_down['pvals']), 
                     alpha=0.8, s=25, color='blue', label=f'Up in Normal ({len(sig_down)})')
            
            # Add labels for top genes
            top_up_genes = sig_up.sort_values('logfoldchanges', ascending=False).head(10)
            top_down_genes = sig_down.sort_values('logfoldchanges', ascending=True).head(10)
            
            for _, gene in top_up_genes.head(5).iterrows():
                plt.text(gene['logfoldchanges'], -np.log10(gene['pvals']), 
                        gene['names'], fontsize=9, ha='center', va='bottom')
                
            for _, gene in top_down_genes.head(5).iterrows():
                plt.text(gene['logfoldchanges'], -np.log10(gene['pvals']), 
                        gene['names'], fontsize=9, ha='center', va='bottom')
            
            plt.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
            plt.axvline(0.5, color='gray', linestyle='--', alpha=0.5)
            plt.axvline(-0.5, color='gray', linestyle='--', alpha=0.5)
            
            plt.xlabel('Log2 Fold Change')
            plt.ylabel('-log10(p-value)')
            plt.title(f'{ct}: Tumor vs Normal (All Ages)')
            plt.legend(loc='best')
            plt.tight_layout()
            plt.savefig(f"{ct_dir}/volcano_tumor_vs_normal.pdf")
            plt.close()
            
            # Add violin plots for top DEGs
            # Top up-regulated genes (higher in tumor)
            top_genes = list(significant_degs.sort_values('logfoldchanges', ascending=False).head(10)['names'])
            if top_genes:
                sc.pl.violin(ct_adata, top_genes, groupby='Tissue_Type', 
                            rotation=90, save=f"_{ct.replace(' ', '_')}_tumor_up.pdf")
                violin_file = f"figures/violin_{ct.replace(' ', '_')}_tumor_up.pdf"
                if os.path.exists(violin_file):
                    os.rename(violin_file, f"{ct_dir}/violin_tumor_up.pdf")
                    
            # Top down-regulated genes (higher in normal)
            bottom_genes = list(significant_degs.sort_values('logfoldchanges', ascending=True).head(10)['names'])
            if bottom_genes:
                sc.pl.violin(ct_adata, bottom_genes, groupby='Tissue_Type', 
                            rotation=90, save=f"_{ct.replace(' ', '_')}_normal_up.pdf")
                violin_file = f"figures/violin_{ct.replace(' ', '_')}_normal_up.pdf"
                if os.path.exists(violin_file):
                    os.rename(violin_file, f"{ct_dir}/violin_normal_up.pdf")
        
        # Now do age-specific analysis
        for age in age_groups:
            # Subset to this age group within this cell type
            age_ct_adata = ct_adata[ct_adata.obs['Age_Group'] == age].copy()
            
            # Check if we have enough cells
            n_tumor = np.sum(age_ct_adata.obs['Tissue_Type'] == 'Tumor')
            n_normal = np.sum(age_ct_adata.obs['Tissue_Type'] == 'Normal')
            
            if n_tumor < min_cells or n_normal < min_cells:
                print(f"  Skipping {ct} in {age} - insufficient cells (normal: {n_normal}, tumor: {n_tumor})")
                celltype_deg_counts.loc[ct, age] = 0
                continue
            
            # Run differential expression
            sc.tl.rank_genes_groups(age_ct_adata, 'Tissue_Type', reference='Normal', method='wilcoxon')
            
            # Get results
            de_results = sc.get.rank_genes_groups_df(age_ct_adata, group='Tumor')
            significant_degs = de_results[
                (de_results['pvals_adj'] < 0.05) & 
                (abs(de_results['logfoldchanges']) > 0.5)
            ]
            
            # Save to file
            age_ct_dir = os.path.join(ct_dir, age.replace(" ", "_").replace("/", "_"))
            os.makedirs(age_ct_dir, exist_ok=True)
            
            de_results.to_csv(f"{age_ct_dir}/tumor_vs_normal_all.csv", index=False)
            significant_degs.to_csv(f"{age_ct_dir}/tumor_vs_normal_significant.csv", index=False)
            
            # Track DEG counts
            n_degs = len(significant_degs)
            celltype_deg_counts.loc[ct, age] = n_degs
            
            # Create volcano plot if we have significant DEGs
            if n_degs > 0:
                plt.figure(figsize=(10, 8))
                plt.scatter(de_results['logfoldchanges'], -np.log10(de_results['pvals']), 
                         alpha=0.5, s=20, color='gray')
                
                # Highlight significant genes
                sig_up = significant_degs[significant_degs['logfoldchanges'] > 0]
                sig_down = significant_degs[significant_degs['logfoldchanges'] < 0]
                
                plt.scatter(sig_up['logfoldchanges'], -np.log10(sig_up['pvals']), 
                         alpha=0.8, s=25, color='red', label=f'Up in Tumor ({len(sig_up)})')
                plt.scatter(sig_down['logfoldchanges'], -np.log10(sig_down['pvals']), 
                         alpha=0.8, s=25, color='blue', label=f'Up in Normal ({len(sig_down)})')
                
                # Add labels for top genes
                top_up_genes = sig_up.sort_values('logfoldchanges', ascending=False).head(10)
                top_down_genes = sig_down.sort_values('logfoldchanges', ascending=True).head(10)
                
                for _, gene in top_up_genes.head(5).iterrows():
                    plt.text(gene['logfoldchanges'], -np.log10(gene['pvals']), 
                            gene['names'], fontsize=9, ha='center', va='bottom')
                    
                for _, gene in top_down_genes.head(5).iterrows():
                    plt.text(gene['logfoldchanges'], -np.log10(gene['pvals']), 
                            gene['names'], fontsize=9, ha='center', va='bottom')
                
                plt.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
                plt.axvline(0.5, color='gray', linestyle='--', alpha=0.5)
                plt.axvline(-0.5, color='gray', linestyle='--', alpha=0.5)
                
                plt.xlabel('Log2 Fold Change')
                plt.ylabel('-log10(p-value)')
                plt.title(f'{ct}: Tumor vs Normal in {age}')
                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig(f"{age_ct_dir}/volcano_tumor_vs_normal.pdf")
                plt.close()
                
                # Add violin plots for top DEGs
                # We need to create a subset for just this cell type and age group
                age_ct_subset = ct_adata[ct_adata.obs['Age_Group'] == age].copy()
                
                # Top up-regulated genes (higher in tumor)
                top_genes = list(significant_degs.sort_values('logfoldchanges', ascending=False).head(10)['names'])
                if top_genes:
                    sc.pl.violin(age_ct_subset, top_genes, groupby='Tissue_Type', 
                                rotation=90, save=f"_{ct.replace(' ', '_')}_{age.replace(' ', '_')}_tumor_up.pdf")
                    violin_file = f"figures/violin_{ct.replace(' ', '_')}_{age.replace(' ', '_')}_tumor_up.pdf"
                    if os.path.exists(violin_file):
                        os.rename(violin_file, f"{age_ct_dir}/violin_tumor_up.pdf")
                        
                # Top down-regulated genes (higher in normal)
                bottom_genes = list(significant_degs.sort_values('logfoldchanges', ascending=True).head(10)['names'])
                if bottom_genes:
                    sc.pl.violin(age_ct_subset, bottom_genes, groupby='Tissue_Type', 
                                rotation=90, save=f"_{ct.replace(' ', '_')}_{age.replace(' ', '_')}_normal_up.pdf")
                    violin_file = f"figures/violin_{ct.replace(' ', '_')}_{age.replace(' ', '_')}_normal_up.pdf"
                    if os.path.exists(violin_file):
                        os.rename(violin_file, f"{age_ct_dir}/violin_normal_up.pdf")
    
    # Save summary tables
    age_deg_counts.to_csv(f"{output_dir}/age_specific_deg_counts.csv")
    celltype_deg_counts.to_csv(f"{output_dir}/celltype_age_deg_counts.csv")
    
    # Create summary heatmap for cell types
    celltype_deg_counts = celltype_deg_counts.fillna(0).astype(int)
    celltype_deg_counts = celltype_deg_counts.loc[celltype_deg_counts.sum(axis=1) > 0]
    
    if len(celltype_deg_counts) > 0:
        plt.figure(figsize=(12, 10))
        import seaborn as sns
        sns.heatmap(celltype_deg_counts, annot=True, fmt="d", cmap="YlOrRd")
        plt.title("Number of Tumor vs Normal DEGs by Cell Type and Age Group")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/celltype_age_degs_heatmap.pdf")
    
    return celltype_deg_counts

def main():
    parser = argparse.ArgumentParser(description='Tumor vs normal differential expression analysis')
    parser.add_argument('-i', '--input', default="./age_analysis/age_celltype_annotation.h5ad",
                       help="Path to annotated AnnData file (.h5ad)")
    parser.add_argument('-o', '--output', default="./age_analysis/tumor_vs_normal",
                       help="Output directory for results")
    parser.add_argument('-m', '--min_cells', type=int, default=20,
                       help="Minimum cells required for each group")
    args = parser.parse_args()

    # Load the annotated data
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells with {adata.n_vars} genes")
    
    # Run tumor vs normal analysis
    summary = analyze_tumor_vs_normal(adata, args.output, args.min_cells)
    
    print("\nAnalysis complete!")
    if len(summary) > 0:
        print("\nSummary of differential expression analysis:")
        print(summary)
    else:
        print("No significant differentially expressed genes found between tumor and normal tissues.")

if __name__ == "__main__":
    main()
