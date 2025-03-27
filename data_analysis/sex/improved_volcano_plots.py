#!/usr/bin/env python3
"""
Improved Volcano Plots
Generates publication-quality volcano plots for differential expression analysis,
with non-overlapping labels for marker genes.
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text  # For non-overlapping text labels

def create_improved_volcano(de_results_file, output_file, marker_genes, title, min_log2fc=0.5, max_padj=0.05):
    """Create a volcano plot with non-overlapping text from saved DE results"""
    # Check if the file exists
    if not os.path.exists(de_results_file):
        print(f"Error: Results file not found: {de_results_file}")
        return
    
    # Load the DE results
    de_results = pd.read_csv(de_results_file)
    
    # Filter significant DEGs
    significant_degs = de_results[
        (de_results['pvals_adj'] < max_padj) & 
        (abs(de_results['logfoldchanges']) > min_log2fc)
    ]
    
    # Create volcano plot
    plt.figure(figsize=(12, 10))
    epsilon = 1e-300  # To avoid log10(0)
    
    # Create safe log10 p-values for plotting
    de_results['safe_log10_pval'] = -np.log10(np.maximum(de_results['pvals'].values, epsilon))
    
    # Plot all genes in gray
    plt.scatter(de_results['logfoldchanges'], de_results['safe_log10_pval'], 
              alpha=0.4, s=15, color='gray', label="Not significant")
    
    # Highlight significant genes
    sig_up = significant_degs[significant_degs['logfoldchanges'] > 0].copy()
    sig_down = significant_degs[significant_degs['logfoldchanges'] < 0].copy()
    
    sig_up['safe_log10_pval'] = -np.log10(np.maximum(sig_up['pvals'].values, epsilon))
    sig_down['safe_log10_pval'] = -np.log10(np.maximum(sig_down['pvals'].values, epsilon))
    
    plt.scatter(sig_up['logfoldchanges'], sig_up['safe_log10_pval'], 
              alpha=0.7, s=25, color='red', label=f'Up-regulated ({len(sig_up)})')
    plt.scatter(sig_down['logfoldchanges'], sig_down['safe_log10_pval'], 
              alpha=0.7, s=25, color='blue', label=f'Down-regulated ({len(sig_down)})')
    
    # Find marker genes in the results
    texts = []  # For adjust_text
    for gene, description in marker_genes.items():
        if gene in de_results['names'].values:
            gene_data = de_results[de_results['names'] == gene].iloc[0]
            is_significant = gene_data['pvals_adj'] < max_padj and abs(gene_data['logfoldchanges']) > min_log2fc
            
            # Determine marker color
            color = 'green' if is_significant else 'purple'
            
            # Plot marker
            plt.scatter(gene_data['logfoldchanges'], gene_data['safe_log10_pval'], s=150, 
                      color=color, marker='*', edgecolor='black', linewidth=1,
                      alpha=0.8, zorder=10)
            
            # Add to texts for adjustment
            text = plt.text(gene_data['logfoldchanges'], gene_data['safe_log10_pval'],
                         gene, fontsize=10, fontweight='bold', ha='center',
                         va='center', bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.2'))
            texts.append(text)
    
    # Adjust text positions to avoid overlap
    if texts:
        adjust_text(texts, 
                   arrowprops=dict(arrowstyle='-', color='black', lw=0.5),
                   expand_points=(1.5, 1.5),
                   force_text=(0.5, 0.5))
    
    # Add reference lines
    plt.axhline(-np.log10(max_padj), linestyle='--', color='gray', alpha=0.7)
    plt.axvline(-min_log2fc, linestyle='--', color='gray', alpha=0.7)
    plt.axvline(min_log2fc, linestyle='--', color='gray', alpha=0.7)
    
    # Format plot
    plt.xlabel('Log2 Fold Change', fontsize=14)
    plt.ylabel('-log10(p-value)', fontsize=14)
    plt.title(title, fontsize=16, pad=20)
    plt.legend(loc='upper left', fontsize=10)
    plt.grid(alpha=0.2)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved {os.path.basename(output_file)}")
    plt.close()


def generate_volcano_plots(results_dir, output_dir):
    """Generate improved volcano plots from existing DE results"""
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
    
    # Find CSV files with DE results
    csv_files = []
    for root, dirs, files in os.walk(results_dir):
        for file in files:
            if file.endswith('.csv') and 'vs' in file:
                csv_files.append(os.path.join(root, file))
    
    print(f"Found {len(csv_files)} DE result files")
    
    # Process each file
    for csv_file in csv_files:
        try:
            # Extract comparison info from filename
            filename = os.path.basename(csv_file)
            base_name = os.path.splitext(filename)[0]
            
            # Create appropriate title from filename
            title_parts = base_name.replace('volcano_', '').split('_')
            title = ' '.join(title_parts).replace('vs', 'vs.').title()
            
            # Create output file
            output_file = os.path.join(output_dir, f"improved_{base_name}.pdf")
            
            print(f"Processing {filename}...")
            create_improved_volcano(
                csv_file, 
                output_file, 
                marker_genes, 
                title
            )
        except Exception as e:
            print(f"  Error processing {csv_file}: {str(e)}")
    
    print(f"All improved volcano plots saved to {output_dir}")
    return


def main():
    parser = argparse.ArgumentParser(description='Generate improved volcano plots from DE results')
    parser.add_argument('-i', '--input_dir', required=True, help='Directory containing DE result CSV files')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for plots')
    args = parser.parse_args()
    
    # Generate improved volcano plots
    generate_volcano_plots(args.input_dir, args.output_dir)

if __name__ == "__main__":
    main()
