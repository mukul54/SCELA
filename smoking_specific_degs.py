#!/usr/bin/env python3
"""
Smoking-specific differential expression analysis with gender stratification.
Identifies genes differentially expressed between smokers and non-smokers,
stratified by gender, age group, cell type, and tissue origin.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse
from tumor_normal_degs import classify_tissue_type

def analyze_smoking_differences(adata, output_dir, gender=None, tissue_type=None, min_cells=40):
    """
    Perform differential expression analysis between smokers and non-smokers.
    Can be filtered by tissue type and gender.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated single-cell data
    output_dir : str
        Directory to save results
    gender : str or None
        If specified, limits analysis to 'Male' or 'Female'
    tissue_type : str or None
        If specified, limits analysis to 'Tumor' or 'Normal' tissues
    min_cells : int
        Minimum number of cells required in each group for comparison
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Ensure tissue type classification
    if 'Tissue_Type' not in adata.obs.columns:
        adata.obs['Tissue_Type'] = adata.obs['Tissue origins'].apply(classify_tissue_type)
    
    # Print distribution info
    print(f"Smoking distribution:\n{adata.obs['Smoking'].value_counts()}")
    print(f"Gender distribution:\n{adata.obs['Sex'].value_counts()}")
    print(f"Tissue type distribution:\n{adata.obs['Tissue_Type'].value_counts()}")
    
    # Cross-tab of smoking by gender and tissue type
    smoking_crosstab = pd.crosstab([adata.obs['Sex'], adata.obs['Tissue_Type']], 
                               adata.obs['Smoking'])
    smoking_crosstab.to_csv(f"{output_dir}/smoking_distribution.csv")
    
    # Apply filters if specified
    if gender:
        adata = adata[adata.obs['Sex'] == gender].copy()
        print(f"Filtered to {gender}: {adata.n_obs} cells remaining")
    
    if tissue_type:
        adata = adata[adata.obs['Tissue_Type'] == tissue_type].copy()
        print(f"Filtered to {tissue_type} tissue: {adata.n_obs} cells remaining")
    
    # Initialize tracking dataframes
    summary_all = {}
    age_deg_counts = pd.DataFrame()
    celltype_deg_counts = pd.DataFrame()
    
    # 1. Global smoking comparison
    print("\nAnalyzing overall smoking differences...")
    smoking_counts = adata.obs['Smoking'].value_counts()
    
    # Check smoking status values
    print(f"Available smoking statuses: {smoking_counts.index.tolist()}")
    
    # Assuming 'yes'/'no' format, adjust if needed based on actual values
    if len(smoking_counts) < 2:
        print("Error: Need at least two smoking categories")
        return pd.DataFrame()
    
    # Get the most common smoking statuses as reference/comparison
    if len(smoking_counts) >= 2:
        comparison_group = smoking_counts.index[0]  # Most common group (e.g., 'yes')
        reference_group = smoking_counts.index[1]   # Second most common (e.g., 'no')
        
        # Check if we have enough cells and at least 2 samples in each group
        group_sample_counts = {}
        for group in [comparison_group, reference_group]:
            # Count number of distinct samples in this group
            if 'Sample' in adata.obs.columns:
                group_sample_counts[group] = adata[adata.obs['Smoking'] == group].obs['Sample'].nunique()
            else:
                # If no Sample column, assume each cell is a unique sample (conservative approach)
                group_sample_counts[group] = 2 if smoking_counts[group] >= min_cells else 1
        
        # Check if we have enough cells and samples
        if (smoking_counts[comparison_group] < min_cells or 
            smoking_counts[reference_group] < min_cells or
            group_sample_counts[comparison_group] < 2 or 
            group_sample_counts[reference_group] < 2):
            print(f"  Skipping - Insufficient data in smoking groups")
            print(f"  {comparison_group}: {smoking_counts[comparison_group]} cells, {group_sample_counts[comparison_group]} samples")
            print(f"  {reference_group}: {smoking_counts[reference_group]} cells, {group_sample_counts[reference_group]} samples")
            return pd.DataFrame()
            
        print(f"Comparing {comparison_group} (n={smoking_counts[comparison_group]}) vs {reference_group} (n={smoking_counts[reference_group]})")
    
        # Run DE analysis (smokers vs non-smokers)
        sc.tl.rank_genes_groups(adata, 'Smoking', reference=reference_group, method='wilcoxon')
        de_results = sc.get.rank_genes_groups_df(adata, group=comparison_group)
    
    # Filter significant DEGs
    significant_degs = de_results[
        (de_results['pvals_adj'] < 0.05) & 
        (abs(de_results['logfoldchanges']) > 0.5)
    ]
    
    # Save results
    de_results.to_csv(f"{output_dir}/global_smoker_vs_nonsmoker_all.csv", index=False)
    significant_degs.to_csv(f"{output_dir}/global_smoker_vs_nonsmoker_significant.csv", index=False)
    
    # Create volcano plot
    if len(significant_degs) > 0:
        plt.figure(figsize=(10, 8))
        plt.scatter(de_results['logfoldchanges'], -np.log10(de_results['pvals']), 
                 alpha=0.5, s=20, color='gray')
        
        # Highlight significant genes
        sig_up = significant_degs[significant_degs['logfoldchanges'] > 0]
        sig_down = significant_degs[significant_degs['logfoldchanges'] < 0]
        
        plt.scatter(sig_up['logfoldchanges'], -np.log10(sig_up['pvals']), 
                 alpha=0.8, s=25, color='red', label=f'Up in Smokers ({len(sig_up)})')
        plt.scatter(sig_down['logfoldchanges'], -np.log10(sig_down['pvals']), 
                 alpha=0.8, s=25, color='blue', label=f'Up in Non-smokers ({len(sig_down)})')
        
        # Add labels for top genes
        for _, gene in sig_up.sort_values('logfoldchanges', ascending=False).head(5).iterrows():
            plt.text(gene['logfoldchanges'], -np.log10(gene['pvals']), 
                    gene['names'], fontsize=9, ha='center', va='bottom')
            
        for _, gene in sig_down.sort_values('logfoldchanges', ascending=True).head(5).iterrows():
            plt.text(gene['logfoldchanges'], -np.log10(gene['pvals']), 
                    gene['names'], fontsize=9, ha='center', va='bottom')
        
        plt.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
        plt.axvline(0.5, color='gray', linestyle='--', alpha=0.5)
        plt.axvline(-0.5, color='gray', linestyle='--', alpha=0.5)
        
        plt.xlabel('Log2 Fold Change')
        plt.ylabel('-log10(p-value)')
        title = f'Smokers vs Non-smokers'
        if gender: title += f' ({gender})'
        if tissue_type: title += f' in {tissue_type} tissue'
        plt.title(title)
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/global_volcano_smoker_vs_nonsmoker.pdf")
        plt.close()
        
        # Create violin plots for top DEGs
        top_genes = list(significant_degs.sort_values('logfoldchanges', ascending=False).head(10)['names'])
        if top_genes:
            sc.pl.violin(adata, top_genes, groupby='Smoking', 
                        rotation=90, save="_global_smoker_up.pdf")
            if os.path.exists("figures/violin_global_smoker_up.pdf"):
                os.rename("figures/violin_global_smoker_up.pdf", 
                         f"{output_dir}/violin_global_{comparison_group}_up.pdf")
        
        bottom_genes = list(significant_degs.sort_values('logfoldchanges', ascending=True).head(10)['names'])
        if bottom_genes:
            sc.pl.violin(adata, bottom_genes, groupby='Smoking', 
                        rotation=90, save="_global_nonsmoker_up.pdf")
            if os.path.exists("figures/violin_global_nonsmoker_up.pdf"):
                os.rename("figures/violin_global_nonsmoker_up.pdf", 
                         f"{output_dir}/violin_global_nonsmoker_up.pdf")
    
    # 2. Age-specific smoking comparison
    age_groups = adata.obs['Age_Group'].unique()
    age_deg_counts = pd.DataFrame(index=['All_Cells'], columns=age_groups)
    
    for age in age_groups:
        print(f"\nAnalyzing smoking differences in age group: {age}")
        age_adata = adata[adata.obs['Age_Group'] == age].copy()
        
        # Check if we have enough cells in both smoking groups
        comparison_values = age_adata.obs['Smoking'].value_counts()
        if len(comparison_values) < 2:
            print(f"  Skipping {age} - insufficient smoking categories")
            age_deg_counts.loc['All_Cells', age] = 0
            continue
            
        comparison_group = comparison_values.index[0]  # e.g., 'yes'
        reference_group = comparison_values.index[1]   # e.g., 'no'
        n_comparison = comparison_values[comparison_group]
        n_reference = comparison_values[reference_group]
        
        # Check for sufficient cells and samples
        group_sample_counts = {}
        for group in [comparison_group, reference_group]:
            # Count number of distinct samples in this group
            if 'Sample' in age_adata.obs.columns:
                group_sample_counts[group] = age_adata[age_adata.obs['Smoking'] == group].obs['Sample'].nunique()
            else:
                # If no Sample column, estimate from cell count
                group_sample_counts[group] = 2 if comparison_values[group] >= min_cells else 1
        
        if (n_comparison < min_cells or n_reference < min_cells or
            group_sample_counts[comparison_group] < 2 or group_sample_counts[reference_group] < 2):
            print(f"  Skipping {age} - insufficient cells or samples")
            print(f"  {comparison_group}: {n_comparison} cells, {group_sample_counts.get(comparison_group, 0)} samples")
            print(f"  {reference_group}: {n_reference} cells, {group_sample_counts.get(reference_group, 0)} samples")
            age_deg_counts.loc['All_Cells', age] = 0
            continue
        
        # Run differential expression
        sc.tl.rank_genes_groups(age_adata, 'Smoking', reference=reference_group, method='wilcoxon')
        de_results = sc.get.rank_genes_groups_df(age_adata, group=comparison_group)
        significant_degs = de_results[
            (de_results['pvals_adj'] < 0.05) & 
            (abs(de_results['logfoldchanges']) > 0.5)
        ]
        
        # Save results
        age_dir = os.path.join(output_dir, age.replace(" ", "_").replace("/", "_"))
        os.makedirs(age_dir, exist_ok=True)
        de_results.to_csv(f"{age_dir}/smoker_vs_nonsmoker_all.csv", index=False)
        significant_degs.to_csv(f"{age_dir}/smoker_vs_nonsmoker_significant.csv", index=False)
        
        # Track DEG counts
        age_deg_counts.loc['All_Cells', age] = len(significant_degs)
        
        # Create volcano plot
        if len(significant_degs) > 0:
            plt.figure(figsize=(10, 8))
            plt.scatter(de_results['logfoldchanges'], -np.log10(de_results['pvals']), 
                     alpha=0.5, s=20, color='gray')
            
            # Highlight significant genes
            sig_up = significant_degs[significant_degs['logfoldchanges'] > 0]
            sig_down = significant_degs[significant_degs['logfoldchanges'] < 0]
            
            plt.scatter(sig_up['logfoldchanges'], -np.log10(sig_up['pvals']), 
                     alpha=0.8, s=25, color='red', label=f'Up in Smokers ({len(sig_up)})')
            plt.scatter(sig_down['logfoldchanges'], -np.log10(sig_down['pvals']), 
                     alpha=0.8, s=25, color='blue', label=f'Up in Non-smokers ({len(sig_down)})')
            
            # Label top genes
            for _, gene in sig_up.sort_values('logfoldchanges', ascending=False).head(5).iterrows():
                plt.text(gene['logfoldchanges'], -np.log10(gene['pvals']), 
                        gene['names'], fontsize=9, ha='center', va='bottom')
                
            for _, gene in sig_down.sort_values('logfoldchanges', ascending=True).head(5).iterrows():
                plt.text(gene['logfoldchanges'], -np.log10(gene['pvals']), 
                        gene['names'], fontsize=9, ha='center', va='bottom')
            
            plt.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
            plt.axvline(0.5, color='gray', linestyle='--', alpha=0.5)
            plt.axvline(-0.5, color='gray', linestyle='--', alpha=0.5)
            
            plt.xlabel('Log2 Fold Change')
            plt.ylabel('-log10(p-value)')
            title = f'Smokers vs Non-smokers: {age}'
            if gender: title += f' ({gender})'
            if tissue_type: title += f' in {tissue_type} tissue'
            plt.title(title)
            plt.legend(loc='best')
            plt.tight_layout()
            plt.savefig(f"{age_dir}/volcano_smoker_vs_nonsmoker.pdf")
            plt.close()
            
            # Create violin plots
            sc.pl.violin(age_adata, 
                        list(significant_degs.sort_values('logfoldchanges', ascending=False).head(10)['names']), 
                        groupby='Smoking', rotation=90, 
                        save=f"_{age.replace(' ', '_')}_smoker_up.pdf")
            
            violin_file = f"figures/violin_{age.replace(' ', '_')}_smoker_up.pdf"
            if os.path.exists(violin_file):
                os.rename(violin_file, f"{age_dir}/violin_smoker_up.pdf")
    
    # 3. Cell type-specific smoking comparison
    cell_types = adata.obs['predicted_cell_type'].unique()
    celltype_deg_counts = pd.DataFrame(index=cell_types, columns=age_groups)
    
    for ct in cell_types:
        print(f"\nAnalyzing cell type: {ct}")
        ct_adata = adata[adata.obs['predicted_cell_type'] == ct].copy()
        
        if ct_adata.n_obs < min_cells:
            print(f"  Skipping {ct} - too few cells overall ({ct_adata.n_obs})")
            continue
        
        # Create directory for cell type
        ct_dir = os.path.join(output_dir, ct.replace(" ", "_").replace("/", "_"))
        os.makedirs(ct_dir, exist_ok=True)
        
        # Check smoking counts
        ct_smoking_counts = ct_adata.obs['Smoking'].value_counts()
        pd.DataFrame(ct_smoking_counts).to_csv(f"{ct_dir}/smoking_counts.csv")
        
        if len(ct_smoking_counts) < 2:
            print(f"  Skipping {ct} - insufficient smoking categories")
            continue
            
        comparison_group = ct_smoking_counts.index[0]  # e.g., 'yes'
        reference_group = ct_smoking_counts.index[1]   # e.g., 'no'
            
        if ct_smoking_counts[comparison_group] < min_cells or ct_smoking_counts[reference_group] < min_cells:
            print(f"  Skipping {ct} - insufficient cells in smoking groups ({comparison_group}: {ct_smoking_counts[comparison_group]}, {reference_group}: {ct_smoking_counts[reference_group]})")
            continue
        
        # Global analysis for this cell type
        sc.tl.rank_genes_groups(ct_adata, 'Smoking', reference=reference_group, method='wilcoxon')
        de_results = sc.get.rank_genes_groups_df(ct_adata, group=comparison_group)
        significant_degs = de_results[
            (de_results['pvals_adj'] < 0.05) & 
            (abs(de_results['logfoldchanges']) > 0.5)
        ]
        
        # Save results
        de_results.to_csv(f"{ct_dir}/smoker_vs_nonsmoker_all.csv", index=False)
        significant_degs.to_csv(f"{ct_dir}/smoker_vs_nonsmoker_significant.csv", index=False)
        
        # Create volcano plot
        if len(significant_degs) > 0:
            plt.figure(figsize=(10, 8))
            plt.scatter(de_results['logfoldchanges'], -np.log10(de_results['pvals']), 
                     alpha=0.5, s=20, color='gray')
            
            # Highlight significant genes
            sig_up = significant_degs[significant_degs['logfoldchanges'] > 0]
            sig_down = significant_degs[significant_degs['logfoldchanges'] < 0]
            
            plt.scatter(sig_up['logfoldchanges'], -np.log10(sig_up['pvals']), 
                     alpha=0.8, s=25, color='red', label=f'Up in Smokers ({len(sig_up)})')
            plt.scatter(sig_down['logfoldchanges'], -np.log10(sig_down['pvals']), 
                     alpha=0.8, s=25, color='blue', label=f'Up in Non-smokers ({len(sig_down)})')
            
            # Label top genes
            for _, gene in sig_up.sort_values('logfoldchanges', ascending=False).head(5).iterrows():
                plt.text(gene['logfoldchanges'], -np.log10(gene['pvals']), 
                        gene['names'], fontsize=9, ha='center', va='bottom')
                
            for _, gene in sig_down.sort_values('logfoldchanges', ascending=True).head(5).iterrows():
                plt.text(gene['logfoldchanges'], -np.log10(gene['pvals']), 
                        gene['names'], fontsize=9, ha='center', va='bottom')
            
            plt.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
            plt.axvline(0.5, color='gray', linestyle='--', alpha=0.5)
            plt.axvline(-0.5, color='gray', linestyle='--', alpha=0.5)
            
            plt.xlabel('Log2 Fold Change')
            plt.ylabel('-log10(p-value)')
            title = f'{ct}: Smokers vs Non-smokers'
            if gender: title += f' ({gender})'
            if tissue_type: title += f' in {tissue_type} tissue'
            plt.title(title)
            plt.legend(loc='best')
            plt.tight_layout()
            plt.savefig(f"{ct_dir}/volcano_smoker_vs_nonsmoker.pdf")
            plt.close()
            
            # Create violin plots
            sc.pl.violin(ct_adata, 
                        list(significant_degs.sort_values('logfoldchanges', ascending=False).head(10)['names']), 
                        groupby='Smoking', rotation=90, 
                        save=f"_{ct.replace(' ', '_')}_smoker_up.pdf")
            
            violin_file = f"figures/violin_{ct.replace(' ', '_')}_smoker_up.pdf"
            if os.path.exists(violin_file):
                os.rename(violin_file, f"{ct_dir}/violin_smoker_up.pdf")
        
        # Age-specific analysis for this cell type
        for age in age_groups:
            age_ct_adata = ct_adata[ct_adata.obs['Age_Group'] == age].copy()
            
            # Check cell counts
            age_ct_counts = age_ct_adata.obs['Smoking'].value_counts()
            if len(age_ct_counts) < 2:
                print(f"  Skipping {ct} in {age} - insufficient smoking categories")
                celltype_deg_counts.loc[ct, age] = 0
                continue
                
            comparison_group = age_ct_counts.index[0]  # Most common value
            reference_group = age_ct_counts.index[1]   # Second most common
            
            # Check for sufficient cells and samples
            group_sample_counts = {}
            for group in [comparison_group, reference_group]:
                # Count number of distinct samples in this group
                if 'Sample' in age_ct_adata.obs.columns:
                    group_sample_counts[group] = age_ct_adata[age_ct_adata.obs['Smoking'] == group].obs['Sample'].nunique()
                else:
                    # If no Sample column, estimate from cell count
                    group_sample_counts[group] = 2 if age_ct_counts[group] >= min_cells else 1
            
            if (age_ct_counts[comparison_group] < min_cells or 
                age_ct_counts[reference_group] < min_cells or
                group_sample_counts[comparison_group] < 2 or 
                group_sample_counts[reference_group] < 2):
                print(f"  Skipping {ct} in {age} - insufficient cells ({comparison_group}: {age_ct_counts[comparison_group]}, {reference_group}: {age_ct_counts[reference_group]})")
                celltype_deg_counts.loc[ct, age] = 0
                continue
            
            # Run differential expression with additional safety check
            try:
                sc.tl.rank_genes_groups(age_ct_adata, 'Smoking', reference=reference_group, method='wilcoxon')
                de_results = sc.get.rank_genes_groups_df(age_ct_adata, group=comparison_group)
            except ValueError as e:
                if "only contain one sample" in str(e):
                    print(f"  Skipping {ct} in {age} - one or more groups only contain one sample")
                    celltype_deg_counts.loc[ct, age] = 0
                    continue
                else:
                    raise
            significant_degs = de_results[
                (de_results['pvals_adj'] < 0.05) & 
                (abs(de_results['logfoldchanges']) > 0.5)
            ]
            
            # Save results
            age_ct_dir = os.path.join(ct_dir, age.replace(" ", "_").replace("/", "_"))
            os.makedirs(age_ct_dir, exist_ok=True)
            de_results.to_csv(f"{age_ct_dir}/smoker_vs_nonsmoker_all.csv", index=False)
            significant_degs.to_csv(f"{age_ct_dir}/smoker_vs_nonsmoker_significant.csv", index=False)
            
            # Track DEG counts
            celltype_deg_counts.loc[ct, age] = len(significant_degs)
            
            # Create volcano plot if we have significant DEGs
            if len(significant_degs) > 0:
                plt.figure(figsize=(10, 8))
                plt.scatter(de_results['logfoldchanges'], -np.log10(de_results['pvals']), 
                         alpha=0.5, s=20, color='gray')
                
                # Highlight significant genes
                sig_up = significant_degs[significant_degs['logfoldchanges'] > 0]
                sig_down = significant_degs[significant_degs['logfoldchanges'] < 0]
                
                plt.scatter(sig_up['logfoldchanges'], -np.log10(sig_up['pvals']), 
                         alpha=0.8, s=25, color='red', label=f'Up in Smokers ({len(sig_up)})')
                plt.scatter(sig_down['logfoldchanges'], -np.log10(sig_down['pvals']), 
                         alpha=0.8, s=25, color='blue', label=f'Up in Non-smokers ({len(sig_down)})')
                
                # Add labels for top genes
                for _, gene in sig_up.sort_values('logfoldchanges', ascending=False).head(5).iterrows():
                    plt.text(gene['logfoldchanges'], -np.log10(gene['pvals']), 
                            gene['names'], fontsize=9, ha='center', va='bottom')
                    
                for _, gene in sig_down.sort_values('logfoldchanges', ascending=True).head(5).iterrows():
                    plt.text(gene['logfoldchanges'], -np.log10(gene['pvals']), 
                            gene['names'], fontsize=9, ha='center', va='bottom')
                
                plt.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
                plt.axvline(0.5, color='gray', linestyle='--', alpha=0.5)
                plt.axvline(-0.5, color='gray', linestyle='--', alpha=0.5)
                
                plt.xlabel('Log2 Fold Change')
                plt.ylabel('-log10(p-value)')
                title = f'{ct}: Smokers vs Non-smokers in {age}'
                if gender: title += f' ({gender})'
                if tissue_type: title += f' in {tissue_type} tissue'
                plt.title(title)
                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig(f"{age_ct_dir}/volcano_smoker_vs_nonsmoker.pdf")
                plt.close()
                
                # Create violin plots
                sc.pl.violin(age_ct_adata, 
                            list(significant_degs.sort_values('logfoldchanges', ascending=False).head(10)['names']), 
                            groupby='Smoking', rotation=90, 
                            save=f"_{ct.replace(' ', '_')}_{age.replace(' ', '_')}_smoker_up.pdf")
                
                violin_file = f"figures/violin_{ct.replace(' ', '_')}_{age.replace(' ', '_')}_smoker_up.pdf"
                if os.path.exists(violin_file):
                    os.rename(violin_file, f"{age_ct_dir}/violin_smoker_up.pdf")
    
    # Save summary tables
    age_deg_counts.to_csv(f"{output_dir}/age_specific_deg_counts.csv")
    celltype_deg_counts.to_csv(f"{output_dir}/celltype_age_deg_counts.csv")
    
    # Create summary heatmap
    celltype_deg_counts = celltype_deg_counts.fillna(0).astype(int)
    celltype_deg_counts = celltype_deg_counts.loc[celltype_deg_counts.sum(axis=1) > 0]
    
    if len(celltype_deg_counts) > 0:
        plt.figure(figsize=(12, 10))
        sns.heatmap(celltype_deg_counts, annot=True, fmt="d", cmap="YlOrRd")
        title = "Number of Smoker vs Non-smoker DEGs by Cell Type and Age Group"
        if gender: title += f" ({gender})"
        if tissue_type: title += f" in {tissue_type} tissue"
        plt.title(title)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/celltype_age_degs_heatmap.pdf")
    
    return celltype_deg_counts

def main():
    parser = argparse.ArgumentParser(description='Smoking-specific differential expression analysis')
    parser.add_argument('-i', '--input', default="./age_analysis/age_celltype_annotation.h5ad",
                       help="Path to annotated AnnData file")
    parser.add_argument('-o', '--output', default="./age_analysis/smoking_specific",
                       help="Output directory for results")
    parser.add_argument('-g', '--gender', default=None, choices=['Male', 'Female'],
                       help="Limit analysis to specific gender")
    parser.add_argument('-t', '--tissue_type', default=None, choices=['Tumor', 'Normal'],
                       help="Limit analysis to specific tissue type")
    parser.add_argument('-m', '--min_cells', type=int, default=40,
                       help="Minimum cells required for each group")
    args = parser.parse_args()

    # Load the annotated data
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells with {adata.n_vars} genes")
    
    # Set up output directory
    base_output = args.output
    
    # Run analysis for specified gender or both genders
    if args.gender:
        # Single gender mode
        gender_output = os.path.join(base_output, args.gender.lower())
        if args.tissue_type:
            gender_output = os.path.join(gender_output, args.tissue_type.lower())
            
        print(f"\nRunning smoking analysis for {args.gender}...")
        analyze_smoking_differences(adata, gender_output, args.gender, args.tissue_type, args.min_cells)
    else:
        # Run for both genders separately and then overall
        # 1. Overall (no gender filter)
        overall_output = base_output
        if args.tissue_type:
            overall_output = os.path.join(overall_output, args.tissue_type.lower())
            
        print("\nRunning overall smoking analysis...")
        analyze_smoking_differences(adata, overall_output, None, args.tissue_type, args.min_cells)
        
        # 2. Male-specific
        male_output = os.path.join(base_output, "male")
        if args.tissue_type:
            male_output = os.path.join(male_output, args.tissue_type.lower())
            
        print("\nRunning smoking analysis for males...")
        analyze_smoking_differences(adata, male_output, "Male", args.tissue_type, args.min_cells)
        
        # 3. Female-specific
        female_output = os.path.join(base_output, "female")
        if args.tissue_type:
            female_output = os.path.join(female_output, args.tissue_type.lower())
            
        print("\nRunning smoking analysis for females...")
        analyze_smoking_differences(adata, female_output, "Female", args.tissue_type, args.min_cells)
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()
