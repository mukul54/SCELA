#!/usr/bin/env python3
"""
Cell-Specific Sex Marker Violin Plots
Generates violin plots for specific marker genes ONLY in their relevant cell types,
showing expression patterns across sex and tumor/normal status.
"""

import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

def generate_cell_specific_violins(adata, output_dir, min_cells=10):
    """Generate violin plots for markers in their relevant cell types"""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Define marker gene to cell type mapping based on the report
    marker_to_cell_types = {
        'SPP1': ['Pulmonary_alveolar_type_II_cells', 'Mast_cells'],
        'KRT19': ['T_memory_cells'],
        'AGR2': ['T_memory_cells', 'Endothelial_cells'],
        'CDC20B': ['Ependymal_cells'],
        'MT1G': ['T_memory_cells'],
        'MT1X': ['T_memory_cells'],
        'MT2A': ['T_memory_cells'],
        'COL11A1': ['Hepatic_stellate_cells'],
        'COL10A1': ['Hepatic_stellate_cells'],
        'SERPING1': ['Dendritic_cells'],
        'GBP1': ['Dendritic_cells'],
        'GBP4': ['Dendritic_cells'],
        'FCGR2B': ['Dendritic_cells'],
        'FOS': ['Dendritic_cells'],
        'SFTPC': ['Pulmonary_alveolar_type_II_cells', 'Myoepithelial_cells'],
        'MALAT1': ['Myoepithelial_cells']
    }
    
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
    
    # Check if we have predicted cell types
    if 'predicted_cell_type' not in adata.obs.columns:
        print("Error: 'predicted_cell_type' column not found in the dataset.")
        return
    
    # Create a combined grouping variable - convert categorical to string first
    adata.obs['Sex_Tissue'] = adata.obs['Sex'].astype(str) + '_' + adata.obs['Tissue_Type'].astype(str)
    
    # Define order and colors for the combined groups
    sex_tissue_order = ['Male_Normal', 'Male_Tumor', 'Female_Normal', 'Female_Tumor']
    sex_tissue_colors = ['#1f77b4', '#1f77b4', '#d62728', '#d62728']  # Blue for male, red for female
    sex_tissue_alpha = [0.5, 0.9, 0.5, 0.9]  # Lower alpha for normal, higher for tumor
    
    # Create palette
    sex_tissue_palette = {}
    for i, group in enumerate(sex_tissue_order):
        r, g, b = [int(sex_tissue_colors[i][j:j+2], 16)/255 for j in (1, 3, 5)]
        sex_tissue_palette[group] = (r, g, b, sex_tissue_alpha[i])
    
    # Process each marker and its relevant cell types
    for gene, cell_types in marker_to_cell_types.items():
        print(f"Processing {gene} in relevant cell types...")
        
        # Skip if gene not in dataset
        if gene not in adata.var_names:
            print(f"  Gene {gene} not found in dataset, skipping")
            continue
        
        # Create a subfolder for this gene
        gene_dir = os.path.join(output_dir, gene)
        os.makedirs(gene_dir, exist_ok=True)
        
        # Process each relevant cell type
        for cell_type in cell_types:
            # Handle possible formatting differences in cell type names
            cell_type_clean = cell_type.replace('_', ' ')
            matching_cell_types = [ct for ct in adata.obs['predicted_cell_type'].unique() 
                                  if ct.lower().replace('_', ' ') == cell_type_clean.lower()]
            
            if not matching_cell_types:
                print(f"  Cell type {cell_type} not found in dataset, skipping")
                continue
                
            # Use the actual cell type name from the dataset
            actual_cell_type = matching_cell_types[0]
            
            # Filter to this cell type
            ct_mask = adata.obs['predicted_cell_type'] == actual_cell_type
            if np.sum(ct_mask) < min_cells:
                print(f"  Insufficient cells for {actual_cell_type}, skipping")
                continue
                
            # Filter data to current cell type
            cell_data = adata[ct_mask]
            
            # Further check if we have enough cells in each group
            counts_by_group = cell_data.obs['Sex_Tissue'].value_counts()
            if any(count < min_cells for count in counts_by_group):
                print(f"  Some groups have fewer than {min_cells} cells in {actual_cell_type}, proceeding with caution")
            
            print(f"  Creating plot for {gene} in {actual_cell_type} (n={np.sum(ct_mask)})")
            
            # Extract data for plotting
            plot_data = pd.DataFrame({
                'Gene_Expression': cell_data[:, gene].X.toarray().flatten(),
                'Sex_Tissue': cell_data.obs['Sex_Tissue'],
                'Sex': cell_data.obs['Sex'],
                'Tissue_Type': cell_data.obs['Tissue_Type'],
                'Cell_Type': actual_cell_type
            })
            
            # 1. CREATE VIOLIN PLOT
            plt.figure(figsize=(10, 6))
            
            # Create violins for each sex-tissue combination
            ax = sns.violinplot(data=plot_data, x='Sex_Tissue', y='Gene_Expression', 
                             order=sex_tissue_order, palette=sex_tissue_palette,
                             scale='width', inner=None)
            
            # Add data points with jitter (small number for clarity)
            sns.stripplot(data=plot_data, x='Sex_Tissue', y='Gene_Expression',
                        order=sex_tissue_order, palette=sex_tissue_palette,
                        size=3, alpha=0.4, jitter=True)
            
            # Calculate and add mean markers
            means = plot_data.groupby('Sex_Tissue')['Gene_Expression'].mean()
            for i, group in enumerate(sex_tissue_order):
                if group in means:
                    plt.plot(i, means[group], 'ko', markersize=10, alpha=0.8)
            
            # Add titles and labels
            plt.title(f"{gene} in {actual_cell_type}\nExpression by Sex and Tissue Type", fontsize=14)
            plt.xlabel('')
            plt.ylabel('Expression')
            
            # Format x-axis labels
            plt.xticks(range(len(sex_tissue_order)), 
                     [f"{x.split('_')[0]} {x.split('_')[1]}" for x in sex_tissue_order],
                     rotation=45)
            
            # Calculate and add statistics
            stat_text = f"Cells: {len(plot_data)}\n"
            for group in sex_tissue_order:
                if group in counts_by_group:
                    stat_text += f"{group}: {counts_by_group[group]} cells\n"
            
            plt.figtext(0.02, 0.02, stat_text, fontsize=8)
            
            plt.tight_layout()
            plt.savefig(os.path.join(gene_dir, f"{gene}_{actual_cell_type.replace(' ', '_')}_violin.pdf"), dpi=300)
            plt.close()
            
            # 2. CREATE BOX PLOT VERSION
            plt.figure(figsize=(10, 6))
            
            # Create box plots
            sns.boxplot(data=plot_data, x='Sex', y='Gene_Expression', hue='Tissue_Type',
                      palette={'Normal': '#66c2a5', 'Tumor': '#fc8d62'})
            
            # Add data points with jitter
            sns.stripplot(data=plot_data, x='Sex', y='Gene_Expression', hue='Tissue_Type',
                        palette={'Normal': '#66c2a5', 'Tumor': '#fc8d62'},
                        dodge=True, alpha=0.3, size=3)
            
            plt.title(f"{gene} in {actual_cell_type}\nExpression by Sex and Tissue Type", fontsize=14)
            plt.xlabel('Sex')
            plt.ylabel('Expression')
            plt.legend(title='Tissue Type')
            
            plt.tight_layout()
            plt.savefig(os.path.join(gene_dir, f"{gene}_{actual_cell_type.replace(' ', '_')}_boxplot.pdf"), dpi=300)
            plt.close()
            
            # 3. CREATE A GROUPED BAR PLOT SHOWING MEAN EXPRESSION
            plt.figure(figsize=(10, 6))
            
            # Calculate means and standard errors
            summary_stats = plot_data.groupby(['Sex', 'Tissue_Type'])['Gene_Expression'].agg(['mean', 'sem', 'count']).reset_index()
            
            # Plot
            ax = sns.barplot(data=summary_stats, x='Sex', y='mean', hue='Tissue_Type',
                          palette={'Normal': '#66c2a5', 'Tumor': '#fc8d62'})
            
            # Add error bars
            for i, row in enumerate(summary_stats.itertuples()):
                plt.errorbar(i, row.mean, yerr=row.sem, color='black', linestyle='none', marker='_', markersize=20)
                
            # Add sample sizes
            # Make sure we only iterate through valid indices
            # We should have 4 patches (2 sexes × 2 tissue types)
            for i, p in enumerate(ax.patches):
                if i < len(summary_stats):
                    row = summary_stats.iloc[i]
                    height = p.get_height()
                    ax.text(p.get_x() + p.get_width()/2., height + 0.1,
                           f'n={row["count"]}', ha='center', fontsize=8)
            
            plt.title(f"{gene} in {actual_cell_type}\nMean Expression by Sex and Tissue Type", fontsize=14)
            plt.xlabel('Sex')
            plt.ylabel('Mean Expression')
            plt.legend(title='Tissue Type')
            
            plt.tight_layout()
            plt.savefig(os.path.join(gene_dir, f"{gene}_{actual_cell_type.replace(' ', '_')}_barplot.pdf"), dpi=300)
            plt.close()
        
    print(f"All cell-specific marker plots saved to {output_dir}")
    return

def main():
    parser = argparse.ArgumentParser(description='Generate cell type-specific violin plots for sex marker genes')
    parser.add_argument('-i', '--input', required=True, help='Input h5ad file')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    parser.add_argument('-m', '--min_cells', type=int, default=10, help='Minimum cells required for plotting')
    args = parser.parse_args()
    
    # Load data
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells with {adata.n_vars} genes")
    
    # Generate cell-specific marker plots
    generate_cell_specific_violins(adata, args.output_dir, args.min_cells)

if __name__ == "__main__":
    main()
