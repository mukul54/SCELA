#!/usr/bin/env python3
"""
Sex-Marker Volcano Plots
Generates volcano plots for differential expression analysis, highlighting key
sex-specific marker genes mentioned in the report.
"""

import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

def generate_volcano_plots(adata, output_dir):
    """Generate volcano plots highlighting marker genes"""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Define the marker genes of interest from the report
    marker_genes = {
        'SPP1': 'Secreted Phosphoprotein 1/Osteopontin',
        'KRT19': 'Keratin 19',
        'AGR2': 'Anterior Gradient 2',
        'CDC20B': 'Cell Division Cycle 20B',
        'MT1G': 'Metallothionein 1G',
        'MT1X': 'Metallothionein 1X',
        'MT2A': 'Metallothionein 2A',
        'COL11A1': 'Collagen Type XI Alpha 1',
        'COL10A1': 'Collagen Type X Alpha 1',
        'SERPING1': 'Serpin Family G Member 1',
        'GBP1': 'Guanylate Binding Protein 1',
        'GBP4': 'Guanylate Binding Protein 4',
        'FCGR2B': 'Fc Fragment of IgG Receptor IIb',
        'FOS': 'Fos Proto-Oncogene',
        'SFTPC': 'Surfactant Protein C',
        'MALAT1': 'Metastasis Associated Lung Adenocarcinoma Transcript 1'
    }
    
    # Ensure tissue type classification like in the original script
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
    
    # Create volcano plots for different comparisons
    
    # 1. Female vs Male (overall)
    print("Creating Female vs Male volcano plot...")
    run_de_analysis_and_create_volcano(
        adata, 
        'Sex', 
        group='Female', 
        reference='Male', 
        highlight_genes=marker_genes, 
        title='Female vs Male (All Cells)', 
        output_file=os.path.join(output_dir, 'volcano_female_vs_male.pdf')
    )
    
    # 2. Tumor vs Normal (overall)
    print("Creating Tumor vs Normal volcano plot...")
    run_de_analysis_and_create_volcano(
        adata, 
        'Tissue_Type', 
        group='Tumor', 
        reference='Normal', 
        highlight_genes=marker_genes, 
        title='Tumor vs Normal (All Cells)', 
        output_file=os.path.join(output_dir, 'volcano_tumor_vs_normal.pdf')
    )
    
    # 3. Female Tumor vs Male Tumor
    print("Creating Female Tumor vs Male Tumor volcano plot...")
    tumor_only = adata[adata.obs['Tissue_Type'] == 'Tumor'].copy()
    if tumor_only.n_obs > 0:
        run_de_analysis_and_create_volcano(
            tumor_only, 
            'Sex', 
            group='Female', 
            reference='Male', 
            highlight_genes=marker_genes, 
            title='Female vs Male (Tumor Only)', 
            output_file=os.path.join(output_dir, 'volcano_female_vs_male_tumor.pdf')
        )
    
    # 4. Female Normal vs Male Normal
    print("Creating Female Normal vs Male Normal volcano plot...")
    normal_only = adata[adata.obs['Tissue_Type'] == 'Normal'].copy()
    if normal_only.n_obs > 0:
        run_de_analysis_and_create_volcano(
            normal_only, 
            'Sex', 
            group='Female', 
            reference='Male', 
            highlight_genes=marker_genes, 
            title='Female vs Male (Normal Only)', 
            output_file=os.path.join(output_dir, 'volcano_female_vs_male_normal.pdf')
        )
    
    # 5. Cell type-specific volcano plots if possible
    if 'predicted_cell_type' in adata.obs.columns:
        print("Creating cell type-specific volcano plots...")
        cell_types = adata.obs['predicted_cell_type'].unique()
        
        for cell_type in cell_types:
            # Only analyze cell types with sufficient cells
            ct_data = adata[adata.obs['predicted_cell_type'] == cell_type].copy()
            if ct_data.n_obs < 100:  # Skip if too few cells
                continue
                
            # Check if both sexes present
            if len(ct_data.obs['Sex'].unique()) < 2:
                continue
                
            print(f"  Processing {cell_type}...")
            run_de_analysis_and_create_volcano(
                ct_data, 
                'Sex', 
                group='Female', 
                reference='Male', 
                highlight_genes=marker_genes, 
                title=f'Female vs Male in {cell_type}', 
                output_file=os.path.join(output_dir, f'volcano_female_vs_male_{cell_type.replace(" ", "_")}.pdf')
            )
    
    print(f"All volcano plots saved to {output_dir}")
    return


def run_de_analysis_and_create_volcano(adata, groupby, group, reference, highlight_genes, title, output_file, min_log2fc=0.5, max_padj=0.05):
    """Run differential expression and create a volcano plot highlighting marker genes"""
    # Check if we have enough cells
    if np.sum(adata.obs[groupby] == group) < 30 or np.sum(adata.obs[groupby] == reference) < 30:
        print(f"  Skipping {group} vs {reference} - insufficient cells")
        return
        
    # Run differential expression
    try:
        sc.tl.rank_genes_groups(adata, groupby, groups=[group], reference=reference, method='wilcoxon')
        de_results = sc.get.rank_genes_groups_df(adata, group=group)
    except Exception as e:
        print(f"  Error running DE analysis for {group} vs {reference}: {str(e)}")
        return
    
    # Save the results
    output_dir = os.path.dirname(output_file)
    os.makedirs(output_dir, exist_ok=True)
    de_results.to_csv(os.path.join(output_dir, f"{os.path.basename(output_file).replace('.pdf', '.csv')}"), index=False)
    
    # Filter significant DEGs
    significant_degs = de_results[
        (de_results['pvals_adj'] < max_padj) & 
        (abs(de_results['logfoldchanges']) > min_log2fc)
    ]
    
    # Create volcano plot
    plt.figure(figsize=(12, 10))
    epsilon = 1e-300  # To avoid log10(0)
    
    # Create safe log10 p-values for plotting
    plot_df = de_results.copy()
    plot_df['safe_log10_pval'] = -np.log10(np.maximum(plot_df['pvals'].values, epsilon))
    
    # Plot all genes in gray
    plt.scatter(plot_df['logfoldchanges'], plot_df['safe_log10_pval'], 
              alpha=0.4, s=15, color='gray', label="Not significant")
    
    # Highlight significant genes
    sig_up = significant_degs[significant_degs['logfoldchanges'] > 0].copy()
    sig_down = significant_degs[significant_degs['logfoldchanges'] < 0].copy()
    
    sig_up['safe_log10_pval'] = -np.log10(np.maximum(sig_up['pvals'].values, epsilon))
    sig_down['safe_log10_pval'] = -np.log10(np.maximum(sig_down['pvals'].values, epsilon))
    
    plt.scatter(sig_up['logfoldchanges'], sig_up['safe_log10_pval'], 
              alpha=0.7, s=25, color='red', label=f'Up in {group} ({len(sig_up)})')
    plt.scatter(sig_down['logfoldchanges'], sig_down['safe_log10_pval'], 
              alpha=0.7, s=25, color='blue', label=f'Down in {group} ({len(sig_down)})')
    
    # Find marker genes in the results and highlight them
    marker_de_genes = []
    for gene in highlight_genes:
        if gene in de_results['names'].values:
            gene_data = de_results[de_results['names'] == gene].iloc[0]
            marker_de_genes.append({
                'gene': gene,
                'logfc': gene_data['logfoldchanges'],
                'pval': gene_data['pvals'],
                'padj': gene_data['pvals_adj'],
                'log10_pval': -np.log10(max(gene_data['pvals'], epsilon)),
                'significant': gene_data['pvals_adj'] < max_padj and abs(gene_data['logfoldchanges']) > min_log2fc
            })
    
    # Plot marker genes with special markers
    for i, gene_data in enumerate(marker_de_genes):
        color = 'green' if gene_data['significant'] else 'purple'
        plt.scatter(gene_data['logfc'], gene_data['log10_pval'], s=150, 
                  color=color, marker='*', edgecolor='black', linewidth=1,
                  alpha=0.8, zorder=10)
        
        # Add labels with position adjustments to avoid overlap
        offset_x = 0.2 if i % 2 == 0 else -0.2
        offset_y = 0.2 if i < len(marker_de_genes)/2 else -0.2
        
        # Adjust text alignment based on position
        ha = 'left' if offset_x > 0 else 'right'
        
        plt.text(gene_data['logfc'] + offset_x, gene_data['log10_pval'] + offset_y,
               gene_data['gene'], fontsize=10, fontweight='bold',
               ha=ha, va='center', bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.2'))
    
    # Add reference lines
    plt.axhline(-np.log10(max_padj), linestyle='--', color='gray', alpha=0.7)
    plt.axvline(-min_log2fc, linestyle='--', color='gray', alpha=0.7)
    plt.axvline(min_log2fc, linestyle='--', color='gray', alpha=0.7)
    
    # Add count of highlighted markers
    if marker_de_genes:
        sig_markers = sum(1 for g in marker_de_genes if g['significant'])
        plt.text(0.02, 0.02, f"Marker genes: {len(marker_de_genes)} total, {sig_markers} significant",
               transform=plt.gca().transAxes, fontsize=10, bbox=dict(facecolor='white', alpha=0.8))
    
    # Format plot
    plt.xlabel('Log2 Fold Change', fontsize=14)
    plt.ylabel('-log10(p-value)', fontsize=14)
    plt.title(title, fontsize=16, pad=20)
    plt.legend(loc='upper left', fontsize=10)
    plt.grid(alpha=0.2)
    
    # Add annotation for significance thresholds
    plt.text(plt.xlim()[1]*0.98, -np.log10(max_padj)*1.1, f'FDR {max_padj}', 
           ha='right', fontsize=9, style='italic', color='gray')
    plt.text(min_log2fc*1.1, plt.ylim()[1]*0.98, f'log2FC {min_log2fc}', 
           va='top', fontsize=9, style='italic', color='gray')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"  Saved {os.path.basename(output_file)}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Generate volcano plots highlighting sex-specific marker genes')
    parser.add_argument('-i', '--input', required=True, help='Input h5ad file')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for plots')
    args = parser.parse_args()
    
    # Load data
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells with {adata.n_vars} genes")
    
    # Generate volcano plots
    generate_volcano_plots(adata, args.output_dir)

if __name__ == "__main__":
    main()
