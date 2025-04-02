#!/usr/bin/env python3
"""
Sex-specific differential expression analysis across age groups and tissue types.
Identifies genes that are differentially expressed between males and females,
stratified by age group, cell type, and tissue origin (tumor/normal).
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from tumor_normal_degs import classify_tissue_type  # Reuse tissue type classifier

def analyze_sex_differences(adata, output_dir, tissue_type=None, min_cells=20):
    """
    Perform differential expression analysis between sexes.
    Can be filtered by tissue type (tumor or normal).
    
    Parameters:
    -----------
    adata : AnnData
        Annotated single-cell data
    output_dir : str
        Directory to save results
    tissue_type : str or None
        If specified, limits analysis to 'Tumor' or 'Normal' tissues
    min_cells : int
        Minimum number of cells required in each group for comparison
    
    Returns:
    --------
    DataFrame
        Summary of DEG counts across conditions
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Add/validate tissue type classification
    if 'Tissue_Type' not in adata.obs.columns:
        adata.obs['Tissue_Type'] = adata.obs['Tissue origins'].apply(classify_tissue_type)
    
    # Print distribution info
    print(f"Sex distribution:\n{adata.obs['Sex'].value_counts()}")
    print(f"Tissue type distribution:\n{adata.obs['Tissue_Type'].value_counts()}")
    
    sex_counts = pd.crosstab(adata.obs['Sex'], adata.obs['Tissue_Type'])
    sex_counts.to_csv(f"{output_dir}/sex_tissue_distribution.csv")
    
    # Filter by tissue type if specified
    if tissue_type:
        adata = adata[adata.obs['Tissue_Type'] == tissue_type].copy()
        print(f"Filtered to {tissue_type} tissue: {adata.n_obs} cells remaining")
        if adata.n_obs < min_cells:
            print(f"Too few cells after filtering to {tissue_type} tissue")
            return pd.DataFrame()
    
    # Initialize results tracking
    results = {}
    
    # 1. Global sex comparison (all age groups, all cell types)
    print("\nAnalyzing overall sex differences...")
    # Verify we have enough cells in each sex group
    sex_cell_counts = adata.obs['Sex'].value_counts()
    
    if len(sex_cell_counts) < 2 or 'Male' not in sex_cell_counts or 'Female' not in sex_cell_counts:
        print("Error: Both male and female groups are required")
        return pd.DataFrame()
        
    if sex_cell_counts['Male'] < min_cells or sex_cell_counts['Female'] < min_cells:
        print(f"Error: Insufficient cells in one or both sex groups (Male: {sex_cell_counts.get('Male', 0)}, Female: {sex_cell_counts.get('Female', 0)})")
        return pd.DataFrame()
    
    # Run differential expression (Female vs Male)
    sc.tl.rank_genes_groups(adata, 'Sex', reference='Male', method='wilcoxon')
    
    # Save results
    de_results = sc.get.rank_genes_groups_df(adata, group='Female')
    significant_degs = de_results[
        (de_results['pvals_adj'] < 0.05) & 
        (abs(de_results['logfoldchanges']) > 0.5)
    ]
    
    de_results.to_csv(f"{output_dir}/global_female_vs_male_all.csv", index=False)
    significant_degs.to_csv(f"{output_dir}/global_female_vs_male_significant.csv", index=False)
    
    # Create global volcano plot
    if len(significant_degs) > 0:
        plt.figure(figsize=(12, 10))
        plt.scatter(de_results['logfoldchanges'], -np.log10(de_results['pvals']), 
                 alpha=0.5, s=20, color='gray')
        
        # Highlight significant genes
        sig_up = significant_degs[significant_degs['logfoldchanges'] > 0]
        sig_down = significant_degs[significant_degs['logfoldchanges'] < 0]
        
        plt.scatter(sig_up['logfoldchanges'], -np.log10(sig_up['pvals']), 
                 alpha=0.8, s=25, color='red', label=f'Up in Female ({len(sig_up)})')
        plt.scatter(sig_down['logfoldchanges'], -np.log10(sig_down['pvals']), 
                 alpha=0.8, s=25, color='blue', label=f'Up in Male ({len(sig_down)})')
        
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
        title_suffix = f" ({tissue_type})" if tissue_type else ""
        plt.title(f'Global Female vs Male Comparison{title_suffix}', fontsize=14)
        plt.legend(loc='best', fontsize=10)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/global_volcano_female_vs_male.pdf")
        plt.close()
        
        # Create violin plots for top DEGs
        top_genes = list(significant_degs.sort_values('logfoldchanges', ascending=False).head(10)['names'])
        sc.pl.violin(adata, top_genes, groupby='Sex', 
                    rotation=90, save="_global_female_vs_male_up.pdf")
        
        bottom_genes = list(significant_degs.sort_values('logfoldchanges', ascending=True).head(10)['names'])
        sc.pl.violin(adata, bottom_genes, groupby='Sex', 
                    rotation=90, save="_global_female_vs_male_down.pdf")
        
        # Move violin plots to the right directory
        if os.path.exists("figures/violin_global_female_vs_male_up.pdf"):
            os.rename("figures/violin_global_female_vs_male_up.pdf", 
                     f"{output_dir}/violin_global_female_vs_male_up.pdf")
        
        if os.path.exists("figures/violin_global_female_vs_male_down.pdf"):
            os.rename("figures/violin_global_female_vs_male_down.pdf", 
                     f"{output_dir}/violin_global_female_vs_male_down.pdf")
    
    # 2. Age-specific sex comparison
    age_groups = adata.obs['Age_Group'].unique()
    
    # Track DEG counts across age groups
    age_deg_counts = pd.DataFrame(index=['All_Cells'], columns=age_groups)
    
    for age in age_groups:
        print(f"\nAnalyzing sex differences in age group: {age}")
        # Subset to this age group
        age_adata = adata[adata.obs['Age_Group'] == age].copy()
        
        # Check if we have enough cells in both sex groups
        n_male = np.sum(age_adata.obs['Sex'] == 'Male')
        n_female = np.sum(age_adata.obs['Sex'] == 'Female')
        
        if n_male < min_cells or n_female < min_cells:
            print(f"  Skipping {age} - insufficient cells (male: {n_male}, female: {n_female})")
            age_deg_counts.loc['All_Cells', age] = 0
            continue
        
        # Run differential expression
        sc.tl.rank_genes_groups(age_adata, 'Sex', reference='Male', method='wilcoxon')
        
        # Get results
        de_results = sc.get.rank_genes_groups_df(age_adata, group='Female')
        significant_degs = de_results[
            (de_results['pvals_adj'] < 0.05) & 
            (abs(de_results['logfoldchanges']) > 0.5)
        ]
        
        # Save to file
        age_dir = os.path.join(output_dir, age.replace(" ", "_").replace("/", "_"))
        os.makedirs(age_dir, exist_ok=True)
        
        de_results.to_csv(f"{age_dir}/female_vs_male_all.csv", index=False)
        significant_degs.to_csv(f"{age_dir}/female_vs_male_significant.csv", index=False)
        
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
                     alpha=0.8, s=25, color='red', label=f'Up in Female ({len(sig_up)})')
            plt.scatter(sig_down['logfoldchanges'], -np.log10(sig_down['pvals']), 
                     alpha=0.8, s=25, color='blue', label=f'Up in Male ({len(sig_down)})')
            
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
            title_suffix = f" ({tissue_type})" if tissue_type else ""
            plt.title(f'Female vs Male: {age}{title_suffix}')
            plt.legend(loc='best')
            plt.tight_layout()
            plt.savefig(f"{age_dir}/volcano_female_vs_male.pdf")
            plt.close()
            
            # Add violin plots for top DEGs
            age_subset = adata[adata.obs['Age_Group'] == age].copy()
            
            # Top up-regulated genes (higher in female)
            top_genes = list(significant_degs.sort_values('logfoldchanges', ascending=False).head(10)['names'])
            if top_genes:
                sc.pl.violin(age_subset, top_genes, groupby='Sex', 
                            rotation=90, save=f"_{age.replace(' ', '_')}_female_up.pdf")
                violin_file = f"figures/violin_{age.replace(' ', '_')}_female_up.pdf"
                if os.path.exists(violin_file):
                    os.rename(violin_file, f"{age_dir}/violin_female_up.pdf")
                    
            # Top down-regulated genes (higher in male)
            bottom_genes = list(significant_degs.sort_values('logfoldchanges', ascending=True).head(10)['names'])
            if bottom_genes:
                sc.pl.violin(age_subset, bottom_genes, groupby='Sex', 
                            rotation=90, save=f"_{age.replace(' ', '_')}_male_up.pdf")
                violin_file = f"figures/violin_{age.replace(' ', '_')}_male_up.pdf"
                if os.path.exists(violin_file):
                    os.rename(violin_file, f"{age_dir}/violin_male_up.pdf")
    
    # 3. Cell type-specific sex comparison
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
        
        # Check sex counts for this cell type
        ct_sex_counts = ct_adata.obs['Sex'].value_counts()
        pd.DataFrame(ct_sex_counts).to_csv(f"{ct_dir}/sex_counts.csv")
        
        if 'Male' not in ct_sex_counts or 'Female' not in ct_sex_counts:
            print(f"  Skipping {ct} - missing either male or female cells")
            continue
            
        if ct_sex_counts['Male'] < min_cells or ct_sex_counts['Female'] < min_cells:
            print(f"  Skipping {ct} - insufficient cells in one sex group")
            continue
        
        # Global analysis for this cell type (all ages)
        sc.tl.rank_genes_groups(ct_adata, 'Sex', reference='Male', method='wilcoxon')
        
        # Get results
        de_results = sc.get.rank_genes_groups_df(ct_adata, group='Female')
        significant_degs = de_results[
            (de_results['pvals_adj'] < 0.05) & 
            (abs(de_results['logfoldchanges']) > 0.5)
        ]
        
        # Save to file
        de_results.to_csv(f"{ct_dir}/female_vs_male_all.csv", index=False)
        significant_degs.to_csv(f"{ct_dir}/female_vs_male_significant.csv", index=False)
        
        # Create volcano plot
        if len(significant_degs) > 0:
            plt.figure(figsize=(10, 8))
            plt.scatter(de_results['logfoldchanges'], -np.log10(de_results['pvals']), 
                     alpha=0.5, s=20, color='gray')
            
            # Highlight significant genes
            sig_up = significant_degs[significant_degs['logfoldchanges'] > 0]
            sig_down = significant_degs[significant_degs['logfoldchanges'] < 0]
            
            plt.scatter(sig_up['logfoldchanges'], -np.log10(sig_up['pvals']), 
                     alpha=0.8, s=25, color='red', label=f'Up in Female ({len(sig_up)})')
            plt.scatter(sig_down['logfoldchanges'], -np.log10(sig_down['pvals']), 
                     alpha=0.8, s=25, color='blue', label=f'Up in Male ({len(sig_down)})')
            
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
            title_suffix = f" ({tissue_type})" if tissue_type else ""
            plt.title(f'{ct}: Female vs Male (All Ages){title_suffix}')
            plt.legend(loc='best')
            plt.tight_layout()
            plt.savefig(f"{ct_dir}/volcano_female_vs_male.pdf")
            plt.close()
            
            # Add violin plots for top DEGs
            # Top up-regulated genes (higher in female)
            top_genes = list(significant_degs.sort_values('logfoldchanges', ascending=False).head(10)['names'])
            if top_genes:
                sc.pl.violin(ct_adata, top_genes, groupby='Sex', 
                            rotation=90, save=f"_{ct.replace(' ', '_')}_female_up.pdf")
                violin_file = f"figures/violin_{ct.replace(' ', '_')}_female_up.pdf"
                if os.path.exists(violin_file):
                    os.rename(violin_file, f"{ct_dir}/violin_female_up.pdf")
                    
            # Top down-regulated genes (higher in male)
            bottom_genes = list(significant_degs.sort_values('logfoldchanges', ascending=True).head(10)['names'])
            if bottom_genes:
                sc.pl.violin(ct_adata, bottom_genes, groupby='Sex', 
                            rotation=90, save=f"_{ct.replace(' ', '_')}_male_up.pdf")
                violin_file = f"figures/violin_{ct.replace(' ', '_')}_male_up.pdf"
                if os.path.exists(violin_file):
                    os.rename(violin_file, f"{ct_dir}/violin_male_up.pdf")
        
        # Now do age-specific analysis
        for age in age_groups:
            # Subset to this age group within this cell type
            age_ct_adata = ct_adata[ct_adata.obs['Age_Group'] == age].copy()
            
            # Check if we have enough cells
            n_male = np.sum(age_ct_adata.obs['Sex'] == 'Male')
            n_female = np.sum(age_ct_adata.obs['Sex'] == 'Female')
            
            if n_male < min_cells or n_female < min_cells:
                print(f"  Skipping {ct} in {age} - insufficient cells (male: {n_male}, female: {n_female})")
                celltype_deg_counts.loc[ct, age] = 0
                continue
            
            # Run differential expression
            sc.tl.rank_genes_groups(age_ct_adata, 'Sex', reference='Male', method='wilcoxon')
            
            # Get results
            de_results = sc.get.rank_genes_groups_df(age_ct_adata, group='Female')
            significant_degs = de_results[
                (de_results['pvals_adj'] < 0.05) & 
                (abs(de_results['logfoldchanges']) > 0.5)
            ]
            
            # Save to file
            age_ct_dir = os.path.join(ct_dir, age.replace(" ", "_").replace("/", "_"))
            os.makedirs(age_ct_dir, exist_ok=True)
            
            de_results.to_csv(f"{age_ct_dir}/female_vs_male_all.csv", index=False)
            significant_degs.to_csv(f"{age_ct_dir}/female_vs_male_significant.csv", index=False)
            
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
                         alpha=0.8, s=25, color='red', label=f'Up in Female ({len(sig_up)})')
                plt.scatter(sig_down['logfoldchanges'], -np.log10(sig_down['pvals']), 
                         alpha=0.8, s=25, color='blue', label=f'Up in Male ({len(sig_down)})')
                
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
                title_suffix = f" ({tissue_type})" if tissue_type else ""
                plt.title(f'{ct}: Female vs Male in {age}{title_suffix}')
                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig(f"{age_ct_dir}/volcano_female_vs_male.pdf")
                plt.close()
                
                # Add violin plots for top DEGs
                # We need to create a subset for just this cell type and age group
                age_ct_subset = ct_adata[ct_adata.obs['Age_Group'] == age].copy()
                
                # Top up-regulated genes (higher in female)
                top_genes = list(significant_degs.sort_values('logfoldchanges', ascending=False).head(10)['names'])
                if top_genes:
                    sc.pl.violin(age_ct_subset, top_genes, groupby='Sex', 
                                rotation=90, save=f"_{ct.replace(' ', '_')}_{age.replace(' ', '_')}_female_up.pdf")
                    violin_file = f"figures/violin_{ct.replace(' ', '_')}_{age.replace(' ', '_')}_female_up.pdf"
                    if os.path.exists(violin_file):
                        os.rename(violin_file, f"{age_ct_dir}/violin_female_up.pdf")
                        
                # Top down-regulated genes (higher in male)
                bottom_genes = list(significant_degs.sort_values('logfoldchanges', ascending=True).head(10)['names'])
                if bottom_genes:
                    sc.pl.violin(age_ct_subset, bottom_genes, groupby='Sex', 
                                rotation=90, save=f"_{ct.replace(' ', '_')}_{age.replace(' ', '_')}_male_up.pdf")
                    violin_file = f"figures/violin_{ct.replace(' ', '_')}_{age.replace(' ', '_')}_male_up.pdf"
                    if os.path.exists(violin_file):
                        os.rename(violin_file, f"{age_ct_dir}/violin_male_up.pdf")
    
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
        title_suffix = f" ({tissue_type})" if tissue_type else ""
        plt.title(f"Number of Female vs Male DEGs by Cell Type and Age Group{title_suffix}")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/celltype_age_degs_heatmap.pdf")
    
    return celltype_deg_counts

def main():
    parser = argparse.ArgumentParser(description='Sex-specific differential expression analysis')
    parser.add_argument('-i', '--input', default="./age_analysis/age_celltype_annotation.h5ad",
                       help="Path to annotated AnnData file (.h5ad)")
    parser.add_argument('-o', '--output', default="./age_analysis/sex_specific",
                       help="Output directory for results")
    parser.add_argument('-t', '--tissue_type', default=None, choices=['Tumor', 'Normal'],
                       help="Limit analysis to specific tissue type (Tumor or Normal)")
    parser.add_argument('-m', '--min_cells', type=int, default=20,
                       help="Minimum cells required for each group")
    args = parser.parse_args()

    # Load the annotated data
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells with {adata.n_vars} genes")
    
    # Adjust output directory if tissue-specific
    if args.tissue_type:
        args.output = os.path.join(args.output, args.tissue_type.lower())
    
    # Run sex-specific analysis
    summary = analyze_sex_differences(adata, args.output, args.tissue_type, args.min_cells)
    
    print("\nAnalysis complete!")
    if len(summary) > 0:
        print("\nSummary of differential expression analysis:")
        print(summary)
    else:
        print("No significant differentially expressed genes found between sexes.")

if __name__ == "__main__":
    main()
