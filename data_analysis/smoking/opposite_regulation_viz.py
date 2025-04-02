#!/usr/bin/env python3
"""
Visualization of genes with opposite regulation patterns across smoking categories.
Creates slope graphs showing the top 5 genes with opposite regulation for Never vs Ex,
Never vs Current, and Ex vs Current smoking comparisons.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import os

# Define output directory
output_dir = "/home/mukul.ranjan/Documents/cancer/SCELA/smoking_tumor_analysis_results"

# Load data for all three comparisons
never_ex_data = pd.read_csv(os.path.join(output_dir, "never_vs_ex_opposite_regulation.csv"))
never_cur_data = pd.read_csv(os.path.join(output_dir, "never_vs_cur_opposite_regulation.csv"))
ex_cur_data = pd.read_csv(os.path.join(output_dir, "ex_vs_cur_opposite_regulation.csv"))

print(f"Loaded data: Never vs Ex: {len(never_ex_data)} genes")
print(f"Loaded data: Never vs Current: {len(never_cur_data)} genes")
print(f"Loaded data: Ex vs Current: {len(ex_cur_data)} genes")

# Create slope graph function for consistent visualization
def create_slope_graph(data1, col1, data2, col2, labels, title, output_file, top_n=5):
    # Sort by magnitude of difference
    data = data1.copy()
    data['magnitude'] = abs(data[col1] - data[col2])
    top_data = data.sort_values('magnitude', ascending=False).head(top_n)
    
    # Create figure and axes
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create positions for genes
    genes = top_data['Gene'].tolist()
    positions = np.arange(len(genes))
    
    # Define x-positions
    x_positions = [1, 2]
    
    # Plot lines for each gene
    for i, gene in enumerate(genes):
        fc1 = top_data.iloc[i][col1]
        fc2 = top_data.iloc[i][col2]
        
        # Color based on direction (red: down-to-up, blue: up-to-down)
        color = 'red' if fc1 < 0 and fc2 > 0 else 'blue'
        
        # Plot line
        ax.plot(x_positions, [fc1, fc2], color=color, alpha=0.8, linewidth=2)
        
        # Add gene labels
        ax.text(x_positions[0]-0.1, fc1, gene, ha='right', va='center', fontsize=12)
        ax.text(x_positions[1]+0.1, fc2, gene, ha='left', va='center', fontsize=12)
    
    # Add labels for x-axis
    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, fontsize=14)
    
    # Add horizontal line at y=0
    ax.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    
    # Add labels and title
    ax.set_xlabel('Smoking Status', fontsize=14)
    ax.set_ylabel('Log2 Fold Change', fontsize=14)
    ax.set_title(title, fontsize=16)
    
    # Create legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='blue', lw=2, label='Up to Down'),
        Line2D([0], [0], color='red', lw=2, label='Down to Up')
    ]
    ax.legend(handles=legend_elements, loc='best', fontsize=12)
    
    # Remove spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Save figure with specified DPI
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f"Saved figure to {output_file}")

# Create Never vs Ex slope graph (top 5)
create_slope_graph(
    never_ex_data, 'Never_log2FC', never_ex_data, 'Ex_log2FC',
    ['Never', 'Ex'],
    'Top 5 Genes with Opposite Regulation: Never vs Ex Smokers',
    os.path.join(output_dir, 'never_ex_slope_graph.pdf')
)

# Create Never vs Current slope graph (top 5)
create_slope_graph(
    never_cur_data, 'Never_log2FC', never_cur_data, 'Cur_log2FC',
    ['Never', 'Current'],
    'Top 5 Genes with Opposite Regulation: Never vs Current Smokers',
    os.path.join(output_dir, 'never_cur_slope_graph.pdf')
)

# Create Ex vs Current slope graph (top 5)
create_slope_graph(
    ex_cur_data, 'Ex_log2FC', ex_cur_data, 'Cur_log2FC',
    ['Ex', 'Current'],
    'Top 5 Genes with Opposite Regulation: Ex vs Current Smokers',
    os.path.join(output_dir, 'ex_cur_slope_graph.pdf')
)

# Create a consolidated heatmap of top 5 genes from all comparisons
def create_comparison_heatmap(output_file):
    # Get top 5 genes from each comparison
    never_ex_data['magnitude'] = abs(never_ex_data['Never_log2FC'] - never_ex_data['Ex_log2FC'])
    never_cur_data['magnitude'] = abs(never_cur_data['Never_log2FC'] - never_cur_data['Cur_log2FC'])
    ex_cur_data['magnitude'] = abs(ex_cur_data['Ex_log2FC'] - ex_cur_data['Cur_log2FC'])
    
    top_never_ex = never_ex_data.sort_values('magnitude', ascending=False).head(5)[['Gene', 'Never_log2FC', 'Ex_log2FC']]
    top_never_cur = never_cur_data.sort_values('magnitude', ascending=False).head(5)[['Gene', 'Never_log2FC', 'Cur_log2FC']]
    top_ex_cur = ex_cur_data.sort_values('magnitude', ascending=False).head(5)[['Gene', 'Ex_log2FC', 'Cur_log2FC']]
    
    # Combine all genes
    all_genes = pd.concat([
        top_never_ex['Gene'], 
        top_never_cur['Gene'], 
        top_ex_cur['Gene']
    ]).drop_duplicates().tolist()
    
    # Create data for the heatmap
    heatmap_data = []
    for gene in all_genes:
        row = {'Gene': gene, 'Never': 0, 'Ex': 0, 'Current': 0}
        
        # Get values from Never vs Ex
        if gene in top_never_ex['Gene'].values:
            gene_row = top_never_ex[top_never_ex['Gene'] == gene].iloc[0]
            row['Never'] = gene_row['Never_log2FC']
            row['Ex'] = gene_row['Ex_log2FC']
        
        # Get values from Never vs Cur
        if gene in top_never_cur['Gene'].values:
            gene_row = top_never_cur[top_never_cur['Gene'] == gene].iloc[0]
            if 'Never' in row and row['Never'] == 0:  # Only override if not set
                row['Never'] = gene_row['Never_log2FC']
            row['Current'] = gene_row['Cur_log2FC']
        
        # Get values from Ex vs Cur
        if gene in top_ex_cur['Gene'].values:
            gene_row = top_ex_cur[top_ex_cur['Gene'] == gene].iloc[0]
            if 'Ex' in row and row['Ex'] == 0:  # Only override if not set
                row['Ex'] = gene_row['Ex_log2FC']
            if 'Current' in row and row['Current'] == 0:  # Only override if not set
                row['Current'] = gene_row['Cur_log2FC']
            
        heatmap_data.append(row)
    
    # Convert to DataFrame and pivot for heatmap
    heatmap_df = pd.DataFrame(heatmap_data)
    heatmap_pivot = heatmap_df.set_index('Gene')[['Never', 'Ex', 'Current']]
    
    # Create the heatmap
    plt.figure(figsize=(8, 10))
    sns.heatmap(heatmap_pivot, cmap="RdBu_r", center=0, 
                annot=True, fmt=".2f", linewidths=0.5, 
                cbar_kws={'label': 'Log2 Fold Change'})
    plt.title('Top Genes with Opposite Regulation Across Smoking Categories', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f"Saved comparison heatmap to {output_file}")

# Create a combined heatmap of top genes
create_comparison_heatmap(os.path.join(output_dir, 'top_smoking_opposite_regulation_heatmap.pdf'))

print("All visualizations complete.")