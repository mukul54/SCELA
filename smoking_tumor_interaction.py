#!/usr/bin/env python3
"""
Smoking-Tumor Interaction Analysis
Identifies genes that show differential expression between tumor and normal tissue
that are influenced by smoking status.
"""

import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

def analyze_smoking_tumor_interaction(adata, output_dir, min_cells=40, analyze_cell_types=True):
    """
    Analyze how smoking status influences tumor vs normal differential expression.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated single-cell data
    output_dir : str
        Directory to save results
    min_cells : int
        Minimum number of cells required for analysis
    analyze_cell_types : bool
        Whether to perform cell type-specific analysis
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/figures", exist_ok=True)
    
    # Ensure tissue type classification
    if 'Tissue_Type' not in adata.obs.columns:
        print("Adding Tissue_Type classification...")
        def classify_tissue_type(tissue_origin):
            """Classify tissue origins into tumor or normal categories"""
            if isinstance(tissue_origin, str):
                if tissue_origin.startswith('n'):
                    return 'Normal'
                elif tissue_origin.startswith('t') or tissue_origin.startswith('m') or tissue_origin == 'PE':
                    return 'Tumor'
            return 'Unknown'
        adata.obs['Tissue_Type'] = adata.obs['Tissue origins'].apply(classify_tissue_type)
        
        # Save tissue classification overview
        tissue_counts = pd.crosstab(adata.obs['Tissue origins'], adata.obs['Tissue_Type'])
        tissue_counts.to_csv(f"{output_dir}/tissue_classification.csv")
    
    # Print distribution info
    print(f"Smoking distribution:\n{adata.obs['Smoking'].value_counts()}")
    print(f"Tissue type distribution:\n{adata.obs['Tissue_Type'].value_counts()}")
    
    # Create crosstab
    smoking_tissue_crosstab = pd.crosstab(adata.obs['Smoking'], adata.obs['Tissue_Type'])
    smoking_tissue_crosstab.to_csv(f"{output_dir}/smoking_tissue_crosstab.csv")
    print(f"\nSmoking by Tissue Type distribution:")
    print(smoking_tissue_crosstab)
    
    # Define smoking groups to analyze
    smoking_groups = adata.obs['Smoking'].unique()
    
    # List to store results for comparison
    group_degs = {}
    all_significant_genes = set()
    
    # For each smoking group, perform tumor vs normal comparison
    for smoking in smoking_groups:
        print(f"\nAnalyzing Tumor vs Normal in {smoking} group")
        # Filter to current smoking group
        smoking_adata = adata[adata.obs['Smoking'] == smoking].copy()
        
        # Check if we have enough cells in both tumor and normal groups
        tumor_count = np.sum(smoking_adata.obs['Tissue_Type'] == 'Tumor')
        normal_count = np.sum(smoking_adata.obs['Tissue_Type'] == 'Normal')
        
        if tumor_count < min_cells or normal_count < min_cells:
            print(f"  Skipping {smoking} - insufficient cells (Tumor: {tumor_count}, Normal: {normal_count})")
            continue
        
        # Run differential expression
        sc.tl.rank_genes_groups(smoking_adata, 'Tissue_Type', reference='Normal', method='wilcoxon')
        de_results = sc.get.rank_genes_groups_df(smoking_adata, group='Tumor')
        
        # Filter significant DEGs
        significant_degs = de_results[
            (de_results['pvals_adj'] < 0.05) & 
            (abs(de_results['logfoldchanges']) > 0.5)
        ]
        
        # Save results
        smoke_dir = os.path.join(output_dir, smoking.replace(" ", "_"))
        os.makedirs(smoke_dir, exist_ok=True)
        de_results.to_csv(f"{smoke_dir}/tumor_vs_normal_all.csv", index=False)
        significant_degs.to_csv(f"{smoke_dir}/tumor_vs_normal_significant.csv", index=False)
        
        # Store for comparison
        group_degs[smoking] = significant_degs
        all_significant_genes.update(significant_degs['names'])
        
        # Create volcano plot and other visualizations only if we have significant DEGs
        if len(significant_degs) > 0:
                
            # Limit to only top 5 genes to avoid overlapping labels
            # Just top 5 by p-value, regardless of fold change direction
            genes_to_label = significant_degs.sort_values('pvals_adj').head(5)
            
            plt.figure(figsize=(12, 10))
            # Add epsilon to avoid log10(0)
            epsilon = 1e-300  # Small value to prevent log10(0)
            
            # Create safe log10 p-values for all data points
            plot_df = de_results.copy()
            plot_df['safe_log10_pval'] = -np.log10(np.maximum(plot_df['pvals'].values, epsilon))
            
            plt.scatter(plot_df['logfoldchanges'], plot_df['safe_log10_pval'], 
                       alpha=0.4, s=15, color='gray')
            
            # Highlight significant genes with safe log10 values
            sig_up = significant_degs[significant_degs['logfoldchanges'] > 0].copy()
            sig_down = significant_degs[significant_degs['logfoldchanges'] < 0].copy()
            
            # Add safe log10 values to avoid division by zero
            sig_up['safe_log10_pval'] = -np.log10(np.maximum(sig_up['pvals'].values, epsilon))
            sig_down['safe_log10_pval'] = -np.log10(np.maximum(sig_down['pvals'].values, epsilon))
            
            plt.scatter(sig_up['logfoldchanges'], sig_up['safe_log10_pval'], 
                       alpha=0.7, s=25, color='red', label=f'Up in Tumor ({len(sig_up)})')
            plt.scatter(sig_down['logfoldchanges'], sig_down['safe_log10_pval'], 
                       alpha=0.7, s=25, color='blue', label=f'Down in Tumor ({len(sig_down)})')
            
            # Add safe log10 values for genes to label
            genes_to_label['safe_log10_pval'] = -np.log10(np.maximum(genes_to_label['pvals'].values, epsilon))
            
            # Highlight genes we're labeling with larger points
            plt.scatter(genes_to_label['logfoldchanges'], genes_to_label['safe_log10_pval'],
                      s=60, color='black', alpha=0.7)
            
            # Apply dynamic label positioning to avoid overlap
            # Use adjustText library if available, otherwise use more separation
            from adjustText import adjust_text
            
            # Add epsilon to avoid log10(0)
            epsilon = 1e-300  # Small value to prevent log10(0)
            
            # Add gene labels with better placement to avoid overlap
            texts = []
            for idx, gene in genes_to_label.iterrows():
                # Create text elements for all genes with safe p-values - use the pre-calculated column
                text = plt.text(gene['logfoldchanges'], gene['safe_log10_pval'],
                             gene['names'], fontsize=12, fontweight='bold',
                             ha='center', va='center')
                texts.append(text)
            
            # Adjust text positions to minimize overlap
            try:
                adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', lw=0.5))
            except (ImportError, NameError):
                # If adjustText is not available, remove all but use simple offsets
                for text in texts:
                    text.remove()
                    
                # Use simple offsets with more space between labels
                for i, (idx, gene) in enumerate(genes_to_label.iterrows()):
                    # Generate offsets that spread out more
                    if i == 0:  # First gene - upper right
                        x_offset, y_offset = 0.3, 0.3
                    elif i == 1:  # Second gene - upper left
                        x_offset, y_offset = -0.3, 0.3
                    elif i == 2:  # Third gene - lower right
                        x_offset, y_offset = 0.3, -0.3
                    elif i == 3:  # Fourth gene - lower left
                        x_offset, y_offset = -0.3, -0.3
                    else:  # Fifth gene - center top with more offset
                        x_offset, y_offset = 0, 0.5
                        
                    plt.text(gene['logfoldchanges'] + x_offset, -np.log10(gene['pvals']) + y_offset,
                           gene['names'], fontsize=10, fontweight='bold',
                           ha='center', va='center', bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray', pad=2))
            
            plt.axhline(-np.log10(0.05), linestyle='--', color='gray')
            plt.axvline(-0.5, linestyle='--', color='gray')
            plt.axvline(0.5, linestyle='--', color='gray')
            
            plt.xlabel('Log2 Fold Change')
            plt.ylabel('-log10(p-value)')
            plt.title(f'Tumor vs Normal in {smoking} Group')
            plt.legend(loc='best')
            plt.tight_layout()
            plt.savefig(f"{smoke_dir}/volcano_tumor_vs_normal.pdf")
            plt.close()
            
            # Create violin plots of top SIGNIFICANT genes only
            if len(significant_degs) >= 5:
                # Get top up-regulated significant genes
                top_genes_up = significant_degs[
                    (significant_degs['pvals_adj'] < 0.05) & 
                    (significant_degs['logfoldchanges'] > 0.5)
                ].sort_values('logfoldchanges', ascending=False).head(5)['names'].tolist()
                
                # Get top down-regulated significant genes
                top_genes_down = significant_degs[
                    (significant_degs['pvals_adj'] < 0.05) & 
                    (significant_degs['logfoldchanges'] < -0.5)
                ].sort_values('logfoldchanges', ascending=True).head(5)['names'].tolist()
                
                # Only create violin plots if we have significant genes
                if top_genes_up:
                    sc.pl.violin(smoking_adata, top_genes_up, groupby='Tissue_Type', 
                               save=f"_{smoking}_tumor_up_genes.pdf")
                    os.rename(f"figures/violin_{smoking}_tumor_up_genes.pdf", 
                            f"{smoke_dir}/violin_tumor_up_genes.pdf")
                    print(f"  Created violin plots for {len(top_genes_up)} up-regulated genes in {smoking} group")
                
                if top_genes_down:
                    sc.pl.violin(smoking_adata, top_genes_down, groupby='Tissue_Type', 
                               save=f"_{smoking}_tumor_down_genes.pdf")
                    os.rename(f"figures/violin_{smoking}_tumor_down_genes.pdf", 
                            f"{smoke_dir}/violin_tumor_down_genes.pdf")
                    print(f"  Created violin plots for {len(top_genes_down)} down-regulated genes in {smoking} group")
                    
                if not top_genes_up and not top_genes_down:
                    print(f"  No genes meet the significance threshold for violin plots in {smoking} group")
    
    # Compare DEGs across smoking groups
    if len(group_degs) >= 2:
        print(f"\nComparing tumor vs normal DEGs across smoking groups")
        
        # Create comparison dataframe
        comparison_df = pd.DataFrame(index=list(all_significant_genes))
        
        for smoking, degs in group_degs.items():
            # Add log fold changes
            fold_changes = degs.set_index('names')['logfoldchanges']
            comparison_df[f"{smoking}_log2FC"] = fold_changes
            
            # Add significance
            pvals = degs.set_index('names')['pvals_adj']
            comparison_df[f"{smoking}_padj"] = pvals
        
        # Fill NaN with zeros for fold changes, high p-values for significance
        comparison_df = comparison_df.fillna({col: 0 if 'log2FC' in col else 1 for col in comparison_df.columns})
        
        # Save comparison
        comparison_df.to_csv(f"{output_dir}/smoking_group_deg_comparison.csv")
        
        # Find smoking-specific tumor markers
        smoking_keys = list(group_degs.keys())
        
        # Find genes with different tumor vs normal effects based on smoking
        differential_effect_genes = []
        
        for gene in all_significant_genes:
            if gene not in comparison_df.index:
                continue
                
            # Check if gene has significant tumor vs normal difference in at least one group
            # and the fold changes differ substantially between groups
            fold_changes = [comparison_df.loc[gene, f"{smoking}_log2FC"] for smoking in smoking_keys]
            padjs = [comparison_df.loc[gene, f"{smoking}_padj"] for smoking in smoking_keys]
            
            # Only include if significant in at least one group
            if not any(p < 0.05 for p in padjs):
                continue
                
            # Check for substantial difference in effect size between groups
            # Looking for genes with at least 1.0 log2FC difference
            max_diff = max(fold_changes) - min(fold_changes)
            if max_diff >= 1.0:
                differential_effect_genes.append({
                    'gene': gene,
                    'max_effect_diff': max_diff,
                    **{f"{smoking}_log2FC": comparison_df.loc[gene, f"{smoking}_log2FC"] for smoking in smoking_keys},
                    **{f"{smoking}_padj": comparison_df.loc[gene, f"{smoking}_padj"] for smoking in smoking_keys}
                })
        
        # Create dataframe of genes with differential effects
        if differential_effect_genes:
            diff_effect_df = pd.DataFrame(differential_effect_genes)
            diff_effect_df = diff_effect_df.sort_values('max_effect_diff', ascending=False)
            diff_effect_df.to_csv(f"{output_dir}/smoking_dependent_tumor_genes.csv", index=False)
            
            # Visualize top differential effect genes
            top_diff_genes = diff_effect_df['gene'].head(15).tolist()
            
            # Store this for later reference to avoid the NameError
            top5_genes = diff_effect_df['gene'].head(5).tolist()
            
            # Create heatmap of log fold changes
            fold_change_matrix = comparison_df.loc[top_diff_genes, [col for col in comparison_df.columns if 'log2FC' in col]]
            fold_change_matrix.columns = [col.replace('_log2FC', '') for col in fold_change_matrix.columns]
            
            # 1. Enhanced heatmap with clearer presentation
            plt.figure(figsize=(12, 10))
            cmap = sns.diverging_palette(240, 10, as_cmap=True)
            ax = sns.heatmap(fold_change_matrix, cmap=cmap, center=0, 
                      annot=True, fmt='.2f', linewidths=0.8, 
                      cbar_kws={'label': 'Log2 Fold Change (Tumor vs Normal)'})
            
            # Improve readability of gene names
            plt.yticks(fontsize=11, fontweight='bold')
            plt.xticks(fontsize=12, fontweight='bold')
            
            plt.title('Smoking-Dependent Tumor Gene Expression Patterns', fontsize=14, fontweight='bold')
            plt.ylabel('Genes', fontsize=12, fontweight='bold')
            plt.xlabel('Smoking Status', fontsize=12, fontweight='bold')
            
            plt.tight_layout()
            plt.savefig(f"{output_dir}/figures/smoking_dependent_tumor_genes_heatmap.pdf", bbox_inches='tight', dpi=300)
            
            # 2. Also create a simplified barplot of the top 5 genes showing their effect size differences
            plt.figure(figsize=(14, 8))
            
            # For the top 5 genes, show their log2FC across smoking statuses
            top5_matrix = fold_change_matrix.loc[top5_genes]
            top5_melted = top5_matrix.reset_index().melt(id_vars='index', var_name='Smoking', value_name='Log2FC')
            top5_melted.rename(columns={'index': 'Gene'}, inplace=True)
            
            # Create grouped bar plot
            sns.barplot(x='Gene', y='Log2FC', hue='Smoking', data=top5_melted, palette='Set2')
            
            plt.axhline(y=0, color='black', linestyle='-', alpha=0.3)
            plt.title('Top 5 Genes with Smoking-Dependent Tumor Effects', fontsize=14, fontweight='bold')
            plt.ylabel('Log2 Fold Change (Tumor vs Normal)', fontsize=12)
            plt.xlabel('Gene', fontsize=12)
            plt.legend(title='Smoking Status', title_fontsize=12)
            plt.xticks(fontsize=11, fontweight='bold')
            
            plt.tight_layout()
            plt.savefig(f"{output_dir}/figures/top5_smoking_dependent_genes_barplot.pdf", bbox_inches='tight', dpi=300)
            
            plt.close('all')
            
            print(f"\nFound {len(differential_effect_genes)} genes with smoking-dependent tumor effects")
            print(f"Top genes: {', '.join(top_diff_genes[:5])}")

def analyze_smoking_tumor_celltype_interactions(adata, output_dir, min_cells=40):
    """
    Analyze how smoking status influences tumor vs normal differential expression for each cell type.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated single-cell data
    output_dir : str
        Directory to save results
    min_cells : int
        Minimum number of cells required for analysis
    """
    print("\nPerforming cell type-specific smoking-tumor interaction analysis...")
    
    # Create cell type output directory
    celltype_dir = os.path.join(output_dir, "cell_types")
    os.makedirs(celltype_dir, exist_ok=True)
    
    # Get unique cell types and smoking groups
    cell_types = adata.obs['predicted_cell_type'].unique()
    smoking_groups = adata.obs['Smoking'].unique()
    
    # Initialize summary DataFrame to track results
    # Create a multi-index DataFrame with cell types and smoking status
    summary_columns = pd.MultiIndex.from_product(
        [smoking_groups, ['Up', 'Down', 'Total']], 
        names=['Smoking', 'Direction']
    )
    summary_df = pd.DataFrame(index=cell_types, columns=summary_columns)
    summary_df = summary_df.fillna(0)
    
    # Process each cell type
    for ct in cell_types:
        print(f"\nAnalyzing cell type: {ct}")
        # Subset to this cell type
        ct_adata = adata[adata.obs['predicted_cell_type'] == ct].copy()
        
        # Skip if too few cells
        if ct_adata.n_obs < min_cells:
            print(f"  Skipping {ct} - too few cells ({ct_adata.n_obs})")
            continue
        
        # Create output directory for this cell type
        ct_output_dir = os.path.join(celltype_dir, ct.replace(" ", "_").replace("/", "_"))
        os.makedirs(ct_output_dir, exist_ok=True)
        
        # Print distribution
        smoke_tissue_counts = pd.crosstab(
            ct_adata.obs['Smoking'], 
            ct_adata.obs['Tissue_Type']
        )
        print(f"  Distribution for {ct}:")
        print(smoke_tissue_counts)
        
        # Save distribution
        smoke_tissue_counts.to_csv(f"{ct_output_dir}/smoking_tissue_counts.csv")
        
        # For each smoking group, perform tumor vs normal analysis
        for smoking in smoking_groups:
            # Filter to current smoking group
            smoking_ct_adata = ct_adata[ct_adata.obs['Smoking'] == smoking].copy()
            
            # Check if we have enough cells in both tumor and normal groups
            tumor_count = np.sum(smoking_ct_adata.obs['Tissue_Type'] == 'Tumor')
            normal_count = np.sum(smoking_ct_adata.obs['Tissue_Type'] == 'Normal')
            
            if tumor_count < min_cells or normal_count < min_cells:
                print(f"  Skipping {smoking} in {ct} - insufficient cells (Tumor: {tumor_count}, Normal: {normal_count})")
                continue
            
            # Create output directory for this smoking group
            smoke_dir = os.path.join(ct_output_dir, smoking.replace(" ", "_"))
            os.makedirs(smoke_dir, exist_ok=True)
            
            print(f"  Analyzing {smoking} smokers: Tumor vs Normal in {ct} cells")
            
            # Run differential expression with additional safety check
            try:
                sc.tl.rank_genes_groups(smoking_ct_adata, 'Tissue_Type', reference='Normal', method='wilcoxon')
                de_results = sc.get.rank_genes_groups_df(smoking_ct_adata, group='Tumor')
            except ValueError as e:
                if "only contain one sample" in str(e):
                    print(f"  Skipping {smoking} in {ct} - insufficient independent samples")
                    continue
                else:
                    raise
            
            # Filter significant DEGs
            significant_degs = de_results[
                (de_results['pvals_adj'] < 0.05) & 
                (abs(de_results['logfoldchanges']) > 0.5)
            ]
            
            # Save results
            de_results.to_csv(f"{smoke_dir}/tumor_vs_normal_all.csv", index=False)
            significant_degs.to_csv(f"{smoke_dir}/tumor_vs_normal_significant.csv", index=False)
            
            # Count up and down regulated genes
            n_up = len(significant_degs[significant_degs['logfoldchanges'] > 0])
            n_down = len(significant_degs[significant_degs['logfoldchanges'] < 0])
            
            # Update summary
            summary_df.loc[ct, (smoking, 'Up')] = n_up
            summary_df.loc[ct, (smoking, 'Down')] = n_down
            summary_df.loc[ct, (smoking, 'Total')] = n_up + n_down
            
            # Only generate visualizations if we have significant DEGs
            if len(significant_degs) > 0:
                # Create volcano plot with limited gene labels (top 5 only)
                plt.figure(figsize=(10, 8))
                # Add epsilon to avoid log10(0)
                epsilon = 1e-300  # Small value to prevent log10(0)
                
                # Create safe log10 p-values for plotting
                plot_df = de_results.copy()
                plot_df['safe_log10_pval'] = -np.log10(np.maximum(plot_df['pvals'].values, epsilon))
                
                plt.scatter(plot_df['logfoldchanges'], plot_df['safe_log10_pval'], 
                         alpha=0.4, s=15, color='gray')
                
                # Highlight significant genes
                sig_up = significant_degs[significant_degs['logfoldchanges'] > 0].copy()
                sig_down = significant_degs[significant_degs['logfoldchanges'] < 0].copy()
                
                # Add safe log10 values
                sig_up['safe_log10_pval'] = -np.log10(np.maximum(sig_up['pvals'].values, epsilon))
                sig_down['safe_log10_pval'] = -np.log10(np.maximum(sig_down['pvals'].values, epsilon))
                
                plt.scatter(sig_up['logfoldchanges'], sig_up['safe_log10_pval'], 
                         alpha=0.7, s=25, color='red', label=f'Up in Tumor ({len(sig_up)})')
                plt.scatter(sig_down['logfoldchanges'], sig_down['safe_log10_pval'], 
                         alpha=0.7, s=25, color='blue', label=f'Down in Tumor ({len(sig_down)})')
                
                # Label just top 5 genes by p-value
                genes_to_label = significant_degs.sort_values('pvals_adj').head(5)
                
                # Create a DataFrame with safe log10 p-values for the genes we're labeling
                genes_to_label['safe_log10_pval'] = -np.log10(np.maximum(genes_to_label['pvals'].values, epsilon))
                
                # Highlight genes we're labeling
                plt.scatter(genes_to_label['logfoldchanges'], genes_to_label['safe_log10_pval'],
                          s=50, color='black', alpha=0.7)
                
                # Add gene labels with manual offsets for clarity
                for i, (idx, gene) in enumerate(genes_to_label.iterrows()):
                    # Generate offsets that spread out labels
                    if i == 0:  # First gene - upper right
                        x_offset, y_offset = 0.3, 0.3
                    elif i == 1:  # Second gene - upper left
                        x_offset, y_offset = -0.3, 0.3
                    elif i == 2:  # Third gene - lower right
                        x_offset, y_offset = 0.3, -0.3
                    elif i == 3:  # Fourth gene - lower left
                        x_offset, y_offset = -0.3, -0.3
                    else:  # Fifth gene - center top with more offset
                        x_offset, y_offset = 0, 0.5
                    
                    # Use the safe log10 p-value for text positioning
                    plt.text(gene['logfoldchanges'] + x_offset, gene['safe_log10_pval'] + y_offset,
                           gene['names'], fontsize=10, fontweight='bold',
                           ha='center', va='center', bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray', pad=2))
                
                plt.axhline(-np.log10(0.05), linestyle='--', color='gray')
                plt.axvline(-0.5, linestyle='--', color='gray')
                plt.axvline(0.5, linestyle='--', color='gray')
                
                plt.xlabel('Log2 Fold Change')
                plt.ylabel('-log10(p-value)')
                plt.title(f'{ct}: Tumor vs Normal in {smoking} Smokers')
                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig(f"{smoke_dir}/volcano_plot.pdf")
                plt.close()
                
                # Generate violin plots for top genes
                # Get top up-regulated significant genes
                top_genes_up = significant_degs[
                    (significant_degs['pvals_adj'] < 0.05) & 
                    (significant_degs['logfoldchanges'] > 0.5)
                ].sort_values('logfoldchanges', ascending=False).head(5)['names'].tolist()
                
                # Get top down-regulated significant genes
                top_genes_down = significant_degs[
                    (significant_degs['pvals_adj'] < 0.05) & 
                    (significant_degs['logfoldchanges'] < -0.5)
                ].sort_values('logfoldchanges', ascending=True).head(5)['names'].tolist()
                
                # Only create violin plots if we have significant genes
                if top_genes_up:
                    sc.pl.violin(smoking_ct_adata, top_genes_up, groupby='Tissue_Type', 
                               save=f"_{ct}_{smoking}_tumor_up_genes.pdf")
                    if os.path.exists(f"figures/violin_{ct}_{smoking}_tumor_up_genes.pdf"):
                        os.rename(f"figures/violin_{ct}_{smoking}_tumor_up_genes.pdf", 
                                f"{smoke_dir}/violin_tumor_up_genes.pdf")
                    print(f"  Created violin plots for {len(top_genes_up)} up-regulated genes in {smoking} group")
                
                if top_genes_down:
                    sc.pl.violin(smoking_ct_adata, top_genes_down, groupby='Tissue_Type', 
                               save=f"_{ct}_{smoking}_tumor_down_genes.pdf")
                    if os.path.exists(f"figures/violin_{ct}_{smoking}_tumor_down_genes.pdf"):
                        os.rename(f"figures/violin_{ct}_{smoking}_tumor_down_genes.pdf", 
                                f"{smoke_dir}/violin_tumor_down_genes.pdf")
                    print(f"  Created violin plots for {len(top_genes_down)} down-regulated genes in {smoking} group")
    
    # Save summary table
    summary_df.to_csv(os.path.join(celltype_dir, "smoking_tumor_celltype_degs_summary.csv"))
    
    # Create summary heatmap
    # Extract just the total columns for the heatmap
    heatmap_data = summary_df.xs('Total', level='Direction', axis=1)
    
    if not heatmap_data.empty and not heatmap_data.isna().all().all():
        plt.figure(figsize=(12, 10))
        # Filter for rows with data
        heatmap_data = heatmap_data.loc[heatmap_data.sum(axis=1) > 0]
        if len(heatmap_data) > 0:
            sns.heatmap(heatmap_data, annot=True, fmt=".0f", cmap="YlOrRd")
            plt.title("Number of Tumor-Normal DEGs by Cell Type and Smoking Status")
            plt.ylabel('Cell Type')
            plt.xlabel('Smoking Status')
            plt.tight_layout()
            plt.savefig(os.path.join(celltype_dir, "smoking_tumor_celltype_degs_heatmap.pdf"))
            plt.close()
    
    return summary_df

def main():
    parser = argparse.ArgumentParser(description='Analyze smoking-tumor interactions')
    parser.add_argument('-i', '--input', required=True, help='Input h5ad file')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    parser.add_argument('-m', '--min_cells', type=int, default=40, help='Minimum cells required')
    parser.add_argument('-c', '--cell_types', action='store_true', help='Perform cell type-specific analysis')
    args = parser.parse_args()
    
    # Load data
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells with {adata.n_vars} genes")
    
    # Run global analysis
    analyze_smoking_tumor_interaction(adata, args.output_dir, args.min_cells, analyze_cell_types=args.cell_types)
    
    # Run cell type-specific analysis if requested
    if args.cell_types:
        analyze_smoking_tumor_celltype_interactions(adata, args.output_dir, args.min_cells)

if __name__ == "__main__":
    main()
