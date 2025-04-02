#!/usr/bin/env python3
"""
Sex-Tumor Interaction Analysis
Identifies genes that show differential expression between tumor and normal tissue
that are influenced by sex (male vs female).
"""

import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

def analyze_sex_tumor_interaction(adata, output_dir, min_cells=40, fold_change_threshold=1.0):
    """Analyze differential tumor vs normal expression based on sex"""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/figures", exist_ok=True)
    
    # Ensure tissue type classification
    if 'Tissue_Type' not in adata.obs.columns:
        def classify_tissue_type(tissue_origin):
            if isinstance(tissue_origin, str):
                if tissue_origin.startswith('n'):
                    return 'Normal'
                elif tissue_origin.startswith('t') or tissue_origin.startswith('m') or tissue_origin == 'PE':
                    return 'Tumor'
            return 'Unknown'
        adata.obs['Tissue_Type'] = adata.obs['Tissue origins'].apply(classify_tissue_type)
    
    # Print distribution 
    print(f"Sex distribution:\n{adata.obs['Sex'].value_counts()}")
    print(f"Tissue type distribution:\n{adata.obs['Tissue_Type'].value_counts()}")
    
    # Crosstab of sex and tissue type
    sex_tissue_crosstab = pd.crosstab(adata.obs['Sex'], adata.obs['Tissue_Type'])
    sex_tissue_crosstab.to_csv(f"{output_dir}/sex_tissue_crosstab.csv")
    print(f"\nSex by Tissue Type distribution:")
    print(sex_tissue_crosstab)
    
    # Define sex groups to analyze
    sex_groups = ['Male', 'Female']
    
    # Store results for comparison
    group_degs = {}
    all_significant_genes = set()
    
    # For each sex, perform tumor vs normal comparison
    for sex in sex_groups:
        print(f"\nAnalyzing Tumor vs Normal in {sex} group")
        # Filter to current sex group
        sex_adata = adata[adata.obs['Sex'] == sex].copy()
        
        # Check if we have enough cells in both tumor and normal groups
        tumor_count = np.sum(sex_adata.obs['Tissue_Type'] == 'Tumor')
        normal_count = np.sum(sex_adata.obs['Tissue_Type'] == 'Normal')
        
        if tumor_count < min_cells or normal_count < min_cells:
            print(f"  Skipping {sex} - insufficient cells (Tumor: {tumor_count}, Normal: {normal_count})")
            continue
        
        # Run differential expression
        sc.tl.rank_genes_groups(sex_adata, 'Tissue_Type', reference='Normal', method='wilcoxon')
        de_results = sc.get.rank_genes_groups_df(sex_adata, group='Tumor')
        
        # Filter significant DEGs
        significant_degs = de_results[
            (de_results['pvals_adj'] < 0.05) & 
            (abs(de_results['logfoldchanges']) > 0.5)
        ]
        
        # Save results
        sex_dir = os.path.join(output_dir, sex)
        os.makedirs(sex_dir, exist_ok=True)
        de_results.to_csv(f"{sex_dir}/tumor_vs_normal_all.csv", index=False)
        significant_degs.to_csv(f"{sex_dir}/tumor_vs_normal_significant.csv", index=False)
        
        # Store for comparison
        group_degs[sex] = significant_degs
        all_significant_genes.update(significant_degs['names'])
        
        # Create volcano plot
        if len(significant_degs) > 0:
            plt.figure(figsize=(10, 8))
            # Add epsilon to avoid log10(0)
            epsilon = 1e-300
            
            # Create safe log10 p-values for plotting
            plot_df = de_results.copy()
            plot_df['safe_log10_pval'] = -np.log10(np.maximum(plot_df['pvals'].values, epsilon))
            
            # Plot all genes
            plt.scatter(plot_df['logfoldchanges'], plot_df['safe_log10_pval'], 
                      alpha=0.4, s=15, color='gray')
            
            # Highlight significant genes
            sig_up = significant_degs[significant_degs['logfoldchanges'] > 0].copy()
            sig_down = significant_degs[significant_degs['logfoldchanges'] < 0].copy()
            
            sig_up['safe_log10_pval'] = -np.log10(np.maximum(sig_up['pvals'].values, epsilon))
            sig_down['safe_log10_pval'] = -np.log10(np.maximum(sig_down['pvals'].values, epsilon))
            
            plt.scatter(sig_up['logfoldchanges'], sig_up['safe_log10_pval'], 
                      alpha=0.7, s=25, color='red', label=f'Up in Tumor ({len(sig_up)})')
            plt.scatter(sig_down['logfoldchanges'], sig_down['safe_log10_pval'], 
                      alpha=0.7, s=25, color='blue', label=f'Down in Tumor ({len(sig_down)})')
            
            # Label top 5 genes by significance
            genes_to_label = significant_degs.sort_values('pvals_adj').head(5)
            genes_to_label['safe_log10_pval'] = -np.log10(np.maximum(genes_to_label['pvals'].values, epsilon))
            
            # Add labels with fixed offsets to avoid overlaps
            for i, (idx, gene) in enumerate(genes_to_label.iterrows()):
                offset_x = 0.2 if i % 2 == 0 else -0.2
                offset_y = 0.2 if i < 2 else -0.2
                plt.text(gene['logfoldchanges'] + offset_x, gene['safe_log10_pval'] + offset_y,
                       gene['names'], fontsize=10, fontweight='bold',
                       ha='center', va='center', bbox=dict(facecolor='white', alpha=0.8))
            
            plt.axhline(-np.log10(0.05), linestyle='--', color='gray')
            plt.axvline(-0.5, linestyle='--', color='gray')
            plt.axvline(0.5, linestyle='--', color='gray')
            
            plt.xlabel('Log2 Fold Change')
            plt.ylabel('-log10(p-value)')
            plt.title(f'Tumor vs Normal in {sex}')
            plt.legend(loc='best')
            plt.tight_layout()
            plt.savefig(f"{sex_dir}/volcano_tumor_vs_normal.pdf")
            plt.close()
            
            # Create violin plots for top DEGs
            if len(significant_degs) >= 5:
                # Get top up-regulated genes
                top_genes_up = significant_degs[
                    significant_degs['logfoldchanges'] > 0
                ].sort_values('logfoldchanges', ascending=False).head(5)['names'].tolist()
                
                # Get top down-regulated genes
                top_genes_down = significant_degs[
                    significant_degs['logfoldchanges'] < 0
                ].sort_values('logfoldchanges', ascending=True).head(5)['names'].tolist()
                
                # Create violin plots if we have significant genes
                if top_genes_up:
                    sc.pl.violin(sex_adata, top_genes_up, groupby='Tissue_Type', 
                               save=f"_{sex}_tumor_up_genes.pdf")
                    if os.path.exists(f"figures/violin_{sex}_tumor_up_genes.pdf"):
                        os.rename(f"figures/violin_{sex}_tumor_up_genes.pdf", 
                                f"{sex_dir}/violin_tumor_up_genes.pdf")
                
                if top_genes_down:
                    sc.pl.violin(sex_adata, top_genes_down, groupby='Tissue_Type', 
                               save=f"_{sex}_tumor_down_genes.pdf")
                    if os.path.exists(f"figures/violin_{sex}_tumor_down_genes.pdf"):
                        os.rename(f"figures/violin_{sex}_tumor_down_genes.pdf", 
                                f"{sex_dir}/violin_tumor_down_genes.pdf")
    
    # Create Venn diagram to show overlap of DEGs between sexes
    if len(group_degs) >= 2:
        try:
            from matplotlib_venn import venn2
            
            # Get sets of significant genes for each sex
            male_genes = set(group_degs.get('Male', pd.DataFrame())['names']) if 'Male' in group_degs else set()
            female_genes = set(group_degs.get('Female', pd.DataFrame())['names']) if 'Female' in group_degs else set()
            
            # Create Venn diagram
            plt.figure(figsize=(8, 8))
            v = venn2([male_genes, female_genes], ('Male', 'Female'))
            plt.title('Overlap of Tumor vs Normal DEGs Between Sexes')
            plt.savefig(f"{output_dir}/figures/sex_degs_venn.pdf")
            plt.close()
            
            # Save intersection and unique genes
            intersection = male_genes.intersection(female_genes)
            male_unique = male_genes - female_genes
            female_unique = female_genes - male_genes
            
            pd.DataFrame({'gene': list(intersection)}).to_csv(
                f"{output_dir}/shared_tumor_degs.csv", index=False)
            pd.DataFrame({'gene': list(male_unique)}).to_csv(
                f"{output_dir}/male_specific_tumor_degs.csv", index=False)
            pd.DataFrame({'gene': list(female_unique)}).to_csv(
                f"{output_dir}/female_specific_tumor_degs.csv", index=False)
            
            print(f"\nDEG overlap analysis:")
            print(f"  Male-specific: {len(male_unique)} genes")
            print(f"  Female-specific: {len(female_unique)} genes")
            print(f"  Shared: {len(intersection)} genes")
        except ImportError:
            print("\nNote: matplotlib_venn not available for Venn diagram visualization")
            print("Install with: pip install matplotlib-venn")
    
    # Compare DEGs across sex groups
    if len(group_degs) >= 2:
        print(f"\nComparing tumor vs normal DEGs across sex groups")
        
        # Create comparison dataframe
        comparison_df = pd.DataFrame(index=list(all_significant_genes))
        
        for sex, degs in group_degs.items():
            # Add log fold changes
            fold_changes = degs.set_index('names')['logfoldchanges']
            comparison_df[f"{sex}_log2FC"] = fold_changes
            
            # Add significance
            pvals = degs.set_index('names')['pvals_adj']
            comparison_df[f"{sex}_padj"] = pvals
        
        # Fill NaN
        comparison_df = comparison_df.fillna({col: 0 if 'log2FC' in col else 1 for col in comparison_df.columns})
        
        # Save comparison
        comparison_df.to_csv(f"{output_dir}/sex_group_deg_comparison.csv")
        
        # Find sex-specific tumor markers
        sex_keys = list(group_degs.keys())
        
        # Find genes with different tumor vs normal effects based on sex
        differential_effect_genes = []
        
        for gene in all_significant_genes:
            if gene not in comparison_df.index:
                continue
                
            # Get fold changes and p-values for each sex
            fold_changes = [comparison_df.loc[gene, f"{sex}_log2FC"] for sex in sex_keys]
            padjs = [comparison_df.loc[gene, f"{sex}_padj"] for sex in sex_keys]
            
            # Only include if significant in at least one group
            if not any(p < 0.05 for p in padjs):
                continue
                
            # Check for substantial difference between male and female
            # Using the specified fold change difference threshold
            max_diff = max(fold_changes) - min(fold_changes)
            if max_diff >= fold_change_threshold:
                differential_effect_genes.append({
                    'gene': gene,
                    'max_effect_diff': max_diff,
                    **{f"{sex}_log2FC": comparison_df.loc[gene, f"{sex}_log2FC"] for sex in sex_keys},
                    **{f"{sex}_padj": comparison_df.loc[gene, f"{sex}_padj"] for sex in sex_keys}
                })
        
        # Create dataframe of genes with differential effects
        if differential_effect_genes:
            diff_effect_df = pd.DataFrame(differential_effect_genes)
            diff_effect_df = diff_effect_df.sort_values('max_effect_diff', ascending=False)
            diff_effect_df.to_csv(f"{output_dir}/sex_dependent_tumor_genes.csv", index=False)
            
            # Create heatmap of top differential effect genes
            top_diff_genes = diff_effect_df['gene'].head(15).tolist()
            
            # Fold change heatmap
            fold_change_matrix = comparison_df.loc[top_diff_genes, [col for col in comparison_df.columns if 'log2FC' in col]]
            fold_change_matrix.columns = [col.replace('_log2FC', '') for col in fold_change_matrix.columns]
            
            # Plot heatmap
            plt.figure(figsize=(10, 8))
            sns.heatmap(fold_change_matrix, cmap=sns.diverging_palette(240, 10, as_cmap=True), 
                      center=0, annot=True, fmt='.2f')
            plt.title('Sex-Dependent Tumor Gene Expression Patterns')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/figures/sex_dependent_tumor_genes_heatmap.pdf")
            plt.close()
            
            # Also create bar plots for top 5 genes
            top5_genes = diff_effect_df['gene'].head(5).tolist()
            top5_matrix = fold_change_matrix.loc[top5_genes]
            top5_melted = top5_matrix.reset_index().melt(id_vars='index', var_name='Sex', value_name='Log2FC')
            top5_melted.rename(columns={'index': 'Gene'}, inplace=True)
            
            plt.figure(figsize=(12, 6))
            sns.barplot(x='Gene', y='Log2FC', hue='Sex', data=top5_melted, palette=['#1f77b4', '#d62728'])
            plt.axhline(y=0, color='black', linestyle='-', alpha=0.3)
            plt.title('Top 5 Genes with Sex-Dependent Tumor Effects')
            plt.ylabel('Log2 Fold Change (Tumor vs Normal)')
            plt.xticks(fontsize=10, fontweight='bold')
            plt.legend(title='Sex')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/figures/top5_sex_dependent_genes_barplot.pdf")
            plt.close()
            
            print(f"\nFound {len(differential_effect_genes)} genes with sex-dependent tumor effects")
            print(f"Top genes: {', '.join(top_diff_genes[:5])}")
    
    return group_degs

def analyze_sex_tumor_celltype_interactions(adata, output_dir, min_cells=40):
    """Analyze cell type-specific sex-tumor interactions"""
    print("\nPerforming cell type-specific sex-tumor interaction analysis...")
    
    # Create output directory
    celltype_dir = os.path.join(output_dir, "cell_types")
    os.makedirs(celltype_dir, exist_ok=True)
    
    # Get unique cell types and sex groups
    cell_types = adata.obs['predicted_cell_type'].unique()
    sex_groups = ['Male', 'Female']
    
    # Summary tracking
    summary_df = pd.DataFrame(index=cell_types, columns=pd.MultiIndex.from_product(
        [sex_groups, ['Up', 'Down', 'Total']], names=['Sex', 'Direction']))
    summary_df = summary_df.fillna(0)
    
    # Process each cell type
    for ct in cell_types:
        print(f"\nAnalyzing cell type: {ct}")
        ct_adata = adata[adata.obs['predicted_cell_type'] == ct].copy()
        
        # Skip if too few cells
        if ct_adata.n_obs < min_cells:
            print(f"  Skipping {ct} - too few cells ({ct_adata.n_obs})")
            continue
        
        # Create output directory
        ct_output_dir = os.path.join(celltype_dir, ct.replace(" ", "_").replace("/", "_"))
        os.makedirs(ct_output_dir, exist_ok=True)
        
        # Distribution
        sex_tissue_counts = pd.crosstab(ct_adata.obs['Sex'], ct_adata.obs['Tissue_Type'])
        sex_tissue_counts.to_csv(f"{ct_output_dir}/sex_tissue_counts.csv")
        
        # Compare tumor vs normal within each sex
        for sex in sex_groups:
            # Filter to current sex
            sex_ct_adata = ct_adata[ct_adata.obs['Sex'] == sex].copy()
            
            # Check if enough cells
            tumor_count = np.sum(sex_ct_adata.obs['Tissue_Type'] == 'Tumor')
            normal_count = np.sum(sex_ct_adata.obs['Tissue_Type'] == 'Normal')
            
            if tumor_count < min_cells or normal_count < min_cells:
                print(f"  Skipping {sex} in {ct} - insufficient cells (Tumor: {tumor_count}, Normal: {normal_count})")
                continue
            
            # Create output directory
            sex_dir = os.path.join(ct_output_dir, sex)
            os.makedirs(sex_dir, exist_ok=True)
            
            # Run differential expression
            try:
                sc.tl.rank_genes_groups(sex_ct_adata, 'Tissue_Type', reference='Normal', method='wilcoxon')
                de_results = sc.get.rank_genes_groups_df(sex_ct_adata, group='Tumor')
            except ValueError as e:
                print(f"  Error analyzing {sex} in {ct}: {str(e)}")
                continue
            
            # Filter significant DEGs
            significant_degs = de_results[
                (de_results['pvals_adj'] < 0.05) & 
                (abs(de_results['logfoldchanges']) > 0.5)
            ]
            
            # Save results
            de_results.to_csv(f"{sex_dir}/tumor_vs_normal_all.csv", index=False)
            significant_degs.to_csv(f"{sex_dir}/tumor_vs_normal_significant.csv", index=False)
            
            # Count up/down regulated
            n_up = len(significant_degs[significant_degs['logfoldchanges'] > 0])
            n_down = len(significant_degs[significant_degs['logfoldchanges'] < 0])
            
            # Update summary
            summary_df.loc[ct, (sex, 'Up')] = n_up
            summary_df.loc[ct, (sex, 'Down')] = n_down
            summary_df.loc[ct, (sex, 'Total')] = n_up + n_down
            
            # Create volcano plot for this cell type
            if len(significant_degs) > 0:
                plt.figure(figsize=(8, 6))
                epsilon = 1e-300  # To avoid log10(0)
                
                # Create safe log10 p-values for plotting
                plot_df = de_results.copy()
                plot_df['safe_log10_pval'] = -np.log10(np.maximum(plot_df['pvals'].values, epsilon))
                
                # Plot all genes in gray
                plt.scatter(plot_df['logfoldchanges'], plot_df['safe_log10_pval'], 
                            alpha=0.4, s=10, color='gray')
                
                # Highlight significant genes
                sig_up = significant_degs[significant_degs['logfoldchanges'] > 0].copy()
                sig_down = significant_degs[significant_degs['logfoldchanges'] < 0].copy()
                
                sig_up['safe_log10_pval'] = -np.log10(np.maximum(sig_up['pvals'].values, epsilon))
                sig_down['safe_log10_pval'] = -np.log10(np.maximum(sig_down['pvals'].values, epsilon))
                
                plt.scatter(sig_up['logfoldchanges'], sig_up['safe_log10_pval'], 
                            alpha=0.7, s=20, color='red', label=f'Up in Tumor ({len(sig_up)})')
                plt.scatter(sig_down['logfoldchanges'], sig_down['safe_log10_pval'], 
                            alpha=0.7, s=20, color='blue', label=f'Down in Tumor ({len(sig_down)})')
                
                # Add reference lines
                plt.axhline(-np.log10(0.05), linestyle='--', color='gray')
                plt.axvline(-0.5, linestyle='--', color='gray')
                plt.axvline(0.5, linestyle='--', color='gray')
                
                plt.xlabel('Log2 Fold Change')
                plt.ylabel('-log10(p-value)')
                plt.title(f'{ct} in {sex}: Tumor vs Normal')
                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig(f"{sex_dir}/volcano_plot.pdf")
                plt.close()
                
                # Create violin plots for top DEGs in this cell type
                if len(significant_degs) >= 3:
                    # Get top genes (up to 3 up and 3 down)
                    top_up = significant_degs[
                        significant_degs['logfoldchanges'] > 0
                    ].sort_values('logfoldchanges', ascending=False).head(3)['names'].tolist()
                    
                    top_down = significant_degs[
                        significant_degs['logfoldchanges'] < 0
                    ].sort_values('logfoldchanges', ascending=True).head(3)['names'].tolist()
                    
                    # Create violin plots
                    if top_up:
                        sc.pl.violin(sex_ct_adata, top_up, groupby='Tissue_Type', 
                                    save=f"_{ct}_{sex}_up_genes.pdf")
                        if os.path.exists(f"figures/violin_{ct}_{sex}_up_genes.pdf"):
                            os.rename(f"figures/violin_{ct}_{sex}_up_genes.pdf", 
                                    f"{sex_dir}/violin_up_genes.pdf")
                    
                    if top_down:
                        sc.pl.violin(sex_ct_adata, top_down, groupby='Tissue_Type', 
                                    save=f"_{ct}_{sex}_down_genes.pdf")
                        if os.path.exists(f"figures/violin_{ct}_{sex}_down_genes.pdf"):
                            os.rename(f"figures/violin_{ct}_{sex}_down_genes.pdf", 
                                    f"{sex_dir}/violin_down_genes.pdf")
    
    # Save summary
    summary_df.to_csv(os.path.join(celltype_dir, "sex_tumor_celltype_summary.csv"))
    
    # Create summary heatmap
    total_data = summary_df.xs('Total', level='Direction', axis=1)
    if not total_data.empty and not total_data.isna().all().all():
        # Filter rows with data
        total_data = total_data.loc[total_data.sum(axis=1) > 0]
        if len(total_data) > 0:
            plt.figure(figsize=(10, 8))
            sns.heatmap(total_data, annot=True, fmt=".0f", cmap="YlOrRd")
            plt.title("Number of Tumor-Normal DEGs by Cell Type and Sex")
            plt.tight_layout()
            plt.savefig(os.path.join(celltype_dir, "sex_tumor_celltype_heatmap.pdf"))
            plt.close()
    
    return summary_df

def main():
    parser = argparse.ArgumentParser(description='Sex-Tumor Interaction Analysis')
    parser.add_argument('-i', '--input', required=True, help='Input h5ad file')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    parser.add_argument('-m', '--min_cells', type=int, default=40, help='Minimum cells required')
    parser.add_argument('-c', '--cell_types', action='store_true', help='Perform cell type-specific analysis')
    parser.add_argument('-f', '--fold_change', type=float, default=1.0, 
                       help='Minimum difference in log2FC between sexes to consider sex-specific')
    args = parser.parse_args()
    
    # Load data
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells with {adata.n_vars} genes")
    
    # Run global analysis
    analyze_sex_tumor_interaction(adata, args.output_dir, args.min_cells, args.fold_change)
    
    # Run cell type-specific analysis if requested
    if args.cell_types:
        analyze_sex_tumor_celltype_interactions(adata, args.output_dir, args.min_cells)

if __name__ == "__main__":
    main()