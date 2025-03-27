#!/usr/bin/env python3
"""
Sex-Tumor Marker Individual Plots
Generates individual violin plots for key sex-specific marker genes identified in the analysis.
"""

import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

def generate_marker_violin_plots(adata, output_dir, min_cells=10):
    """Generate individual violin plots for key sex-specific marker genes"""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Define the marker genes of interest from the report
    primary_markers = {
        'SPP1': 'Secreted Phosphoprotein 1/Osteopontin',
        'KRT19': 'Keratin 19',
        'AGR2': 'Anterior Gradient 2',
        'CDC20B': 'Cell Division Cycle 20B',
        'MT1G': 'Metallothionein 1G',
        'MT1X': 'Metallothionein 1X',
        'MT2A': 'Metallothionein 2A',
        'COL11A1': 'Collagen Type XI Alpha 1',
        'COL10A1': 'Collagen Type X Alpha 1'
    }
    
    # Additional genes with opposite regulation between sexes
    opposite_regulation_markers = {
        'SERPING1': 'Serpin Family G Member 1',
        'GBP1': 'Guanylate Binding Protein 1',
        'GBP4': 'Guanylate Binding Protein 4',
        'FCGR2B': 'Fc Fragment of IgG Receptor IIb',
        'FOS': 'Fos Proto-Oncogene',
        'SFTPC': 'Surfactant Protein C',
        'MALAT1': 'Metastasis Associated Lung Adenocarcinoma Transcript 1'
    }
    
    # Check which markers are in our dataset
    all_markers = {**primary_markers, **opposite_regulation_markers}
    available_markers = [gene for gene in all_markers if gene in adata.var_names]
    missing_markers = [gene for gene in all_markers if gene not in adata.var_names]
    
    if missing_markers:
        print(f"Warning: The following markers are not in the dataset: {', '.join(missing_markers)}")
    
    if not available_markers:
        print("Error: None of the specified marker genes are in the dataset.")
        return
    
    print(f"Generating violin plots for {len(available_markers)} marker genes...")
    
    # Ensure we have the required columns
    if 'Tissue_Type' not in adata.obs.columns:
        print("Error: 'Tissue_Type' column not found in the dataset.")
        return
    
    if 'Sex' not in adata.obs.columns:
        print("Error: 'Sex' column not found in the dataset.")
        return

    # Create a combined grouping variable - convert categorical to string first
    adata.obs['Sex_Tissue'] = adata.obs['Sex'].astype(str) + '_' + adata.obs['Tissue_Type'].astype(str)
    
    # Define colors and order for the combined groups
    sex_tissue_order = ['Male_Normal', 'Male_Tumor', 'Female_Normal', 'Female_Tumor']
    sex_tissue_colors = ['#1f77b4', '#1f77b4', '#d62728', '#d62728']  # Blue for male, red for female
    sex_tissue_alpha = [0.5, 0.9, 0.5, 0.9]  # Lower alpha for normal, higher for tumor
    
    # Create color palette with alpha variations
    sex_tissue_palette = {}
    for i, group in enumerate(sex_tissue_order):
        sex_tissue_palette[group] = sex_tissue_colors[i]
    
    # Create directory for individual plots
    print("Creating individual gene violin plots...")
    
    # Create plots for each gene
    for gene in available_markers:
        if gene not in adata.var_names:
            continue
            
        # Skip genes with very low expression
        cells_with_expression = adata[:, gene].X.toarray().flatten() > 0
        if np.sum(cells_with_expression) < min_cells:
            print(f"Skipping {gene} - insufficient cells with expression")
            continue
            
        # Create data for plotting
        data = pd.DataFrame({
            'Gene_Expression': adata[:, gene].X.toarray().flatten(),
            'Sex_Tissue': adata.obs['Sex_Tissue'],
            'Sex': adata.obs['Sex'],
            'Tissue_Type': adata.obs['Tissue_Type']
        })
        
        plt.figure(figsize=(10, 6))
        
        # Draw violin plot
        ax = sns.violinplot(data=data, x='Sex_Tissue', y='Gene_Expression', 
                         order=sex_tissue_order, palette=sex_tissue_palette,
                         scale='width', inner=None)
        
        # Add stripplot with jitter
        sns.stripplot(data=data, x='Sex_Tissue', y='Gene_Expression',
                    order=sex_tissue_order, palette=sex_tissue_palette,
                    size=3, alpha=0.4, jitter=True)
        
        # Add mean expression points
        means = data.groupby('Sex_Tissue')['Gene_Expression'].mean()
        for j, group in enumerate(sex_tissue_order):
            if group in means:
                plt.plot(j, means[group], 'ko', markersize=10, alpha=0.8)
        
        # Set title and labels
        plt.title(f"{gene} - {all_markers.get(gene, '')}", fontsize=14)
        plt.xlabel('Sex and Tissue Type')
        plt.ylabel('Expression')
        
        # Adjust x-tick labels
        plt.xticks(range(len(sex_tissue_order)), 
                 [f"{x.split('_')[0]} {x.split('_')[1]}" for x in sex_tissue_order],
                 rotation=45)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{gene}_violin.pdf"), dpi=300)
        print(f"  Saved {gene}_violin.pdf")
        plt.close()
        
        # Also create cell type-specific plot for this gene
        if np.sum(cells_with_expression) >= min_cells * 3:  # Require more cells for cell type analysis
            # Identify key cell types based on the report
            key_cell_types = {
                'T_memory_cells': ['KRT19', 'AGR2', 'MT1G'],
                'Hepatic_stellate_cells': ['COL11A1', 'COL10A1'],
                'Endothelial_cells': ['AGR2'],
                'Ependymal_cells': ['CDC20B'],
                'Pulmonary_alveolar_type_II_cells': ['SPP1', 'SFTPC'],
                'Myoepithelial_cells': ['SFTPC', 'MALAT1'],
                'Dendritic_cells': ['SERPING1', 'GBP1', 'GBP4', 'FCGR2B', 'FOS'],
                'Mast_cells': ['SPP1']
            }
            
            # Find cell types for this gene
            gene_cell_types = []
            for ct, ct_genes in key_cell_types.items():
                if gene in ct_genes:
                    gene_cell_types.append(ct)
            
            # If gene is a key marker for specific cell types, create cell type plot
            if gene_cell_types:
                ct_data = []
                for ct in gene_cell_types:
                    # Check if cell type exists
                    if ct not in adata.obs['predicted_cell_type'].unique():
                        continue
                    
                    # Filter to cell type
                    ct_mask = adata.obs['predicted_cell_type'] == ct
                    if np.sum(ct_mask) < min_cells:
                        continue
                        
                    # Get data
                    temp_df = pd.DataFrame({
                        'Expression': adata[ct_mask, gene].X.toarray().flatten(),
                        'Sex': adata[ct_mask].obs['Sex'],
                        'Tissue_Type': adata[ct_mask].obs['Tissue_Type'],
                        'Cell_Type': ct
                    })
                    ct_data.append(temp_df)
                
                # If we have cell type data, create plot
                if ct_data:
                    combined_df = pd.concat(ct_data)
                    
                    # Create facet grid by cell type
                    g = sns.FacetGrid(combined_df, col="Cell_Type", height=4, 
                                     aspect=1.2, sharex=True, sharey=False)
                    
                    # Map violin plot to each facet
                    g.map_dataframe(lambda data, **kws: sns.violinplot(
                        x='Sex', y='Expression', hue='Tissue_Type', 
                        data=data, palette="Set1", **kws))
                    
                    # Add title
                    g.fig.suptitle(f"{gene} expression by cell type", y=1.05, fontsize=14)
                    g.set_axis_labels("Sex", "Expression")
                    g.add_legend()
                    
                    # Save plot
                    plt.tight_layout()
                    plt.savefig(os.path.join(output_dir, f"{gene}_by_celltype.pdf"), dpi=300)
                    print(f"  Saved {gene}_by_celltype.pdf")
                    plt.close()
    
    print(f"All plots saved to {output_dir}")
    return

def main():
    parser = argparse.ArgumentParser(description='Generate individual violin plots for sex-specific marker genes')
    parser.add_argument('-i', '--input', required=True, help='Input h5ad file')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    parser.add_argument('-m', '--min_cells', type=int, default=10, help='Minimum cells required for plotting')
    args = parser.parse_args()
    
    # Load data
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells with {adata.n_vars} genes")
    
    # Generate marker plots
    generate_marker_violin_plots(adata, args.output_dir, args.min_cells)

if __name__ == "__main__":
    main()
