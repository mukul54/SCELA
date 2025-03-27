#!/usr/bin/env python3
"""
Simple Marker Violin Plots
Generates violin plots for key sex-specific marker genes mentioned in the report.
Follows the structure of the original sex_tumor_interaction.py script.
"""

import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

def generate_marker_violins(adata, output_dir, min_cells=10):
    """Generate violin plots for key sex-specific marker genes"""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Define the marker genes of interest from the report
    markers = [
        'SPP1',    # Secreted Phosphoprotein 1/Osteopontin
        'KRT19',   # Keratin 19
        'AGR2',    # Anterior Gradient 2
        'CDC20B',  # Cell Division Cycle 20B
        'MT1G',    # Metallothionein 1G
        'MT1X',    # Metallothionein 1X
        'MT2A',    # Metallothionein 2A
        'COL11A1', # Collagen Type XI Alpha 1
        'COL10A1', # Collagen Type X Alpha 1
        'SERPING1', # Serpin Family G Member 1
        'GBP1',    # Guanylate Binding Protein 1
        'GBP4',    # Guanylate Binding Protein 4
        'FCGR2B',  # Fc Fragment of IgG Receptor IIb
        'FOS',     # Fos Proto-Oncogene
        'SFTPC',   # Surfactant Protein C
        'MALAT1'   # Metastasis Associated Lung Adenocarcinoma Transcript 1
    ]
    
    # Check which markers are available in the dataset
    available_markers = [gene for gene in markers if gene in adata.var_names]
    if not available_markers:
        print("Error: None of the specified marker genes are in the dataset.")
        return
    
    print(f"Generating violin plots for {len(available_markers)} marker genes...")
    
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
    
    # Process each gene
    for gene in available_markers:
        print(f"Processing {gene}...")
        
        # Skip if gene not found
        if gene not in adata.var_names:
            print(f"  Gene {gene} not found in dataset, skipping")
            continue
        
        # 1. Male vs Female overall
        plt.figure(figsize=(10, 6))
        sc.pl.violin(adata, [gene], groupby='Sex', show=False)
        plt.title(f"{gene} - Expression by Sex", fontsize=14)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{gene}_by_sex.pdf"), dpi=300)
        plt.close()
        
        # 2. Tumor vs Normal overall
        plt.figure(figsize=(10, 6))
        sc.pl.violin(adata, [gene], groupby='Tissue_Type', show=False)
        plt.title(f"{gene} - Expression by Tissue Type", fontsize=14)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{gene}_by_tissue.pdf"), dpi=300)
        plt.close()
        
        # 3. Combined Sex-Tissue plot
        # Create combo groups manually to avoid categorical issues
        adata.obs['Sex_Tissue'] = adata.obs['Sex'].astype(str) + "_" + adata.obs['Tissue_Type'].astype(str)
        
        # Prepare data for better plotting
        data = pd.DataFrame({
            'Gene_Expression': adata[:, gene].X.toarray().flatten(),
            'Sex_Tissue': adata.obs['Sex_Tissue'],
            'Sex': adata.obs['Sex'],
            'Tissue_Type': adata.obs['Tissue_Type']
        })
        
        # Define plotting order and colors
        order = ['Male_Normal', 'Male_Tumor', 'Female_Normal', 'Female_Tumor']
        colors = ['#1f77b4', '#1f77b4', '#d62728', '#d62728']  # Blue for male, red for female
        alpha = [0.5, 0.9, 0.5, 0.9]  # Lower alpha for Normal
        palette = dict(zip(order, colors))
        
        plt.figure(figsize=(10, 6))
        
        # Create violin plot
        sns.violinplot(data=data, x='Sex_Tissue', y='Gene_Expression', 
                     order=order, palette=palette,
                     scale='width', inner=None)
        
        # Add individual points
        sns.stripplot(data=data, x='Sex_Tissue', y='Gene_Expression',
                    order=order, palette=palette,
                    size=3, alpha=0.4, jitter=True)
        
        # Mark means
        means = data.groupby('Sex_Tissue')['Gene_Expression'].mean()
        for i, group in enumerate(order):
            if group in means:
                plt.plot(i, means[group], 'ko', markersize=10, alpha=0.8)
        
        # Format plot
        plt.title(f"{gene} - Expression by Sex and Tissue Type", fontsize=14)
        plt.xlabel('')
        plt.ylabel('Expression')
        plt.xticks(range(len(order)), 
                 [f"{x.split('_')[0]} {x.split('_')[1]}" for x in order],
                 rotation=45)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{gene}_combo.pdf"), dpi=300)
        plt.close()
        
        # 4. Cell-type specific plots (if predicted_cell_type is available)
        if 'predicted_cell_type' in adata.obs.columns:
            # Define key cell types for each marker based on the report
            key_cell_types = {
                'SPP1': ['Pulmonary_alveolar_type_II_cells', 'Mast_cells'],
                'KRT19': ['T_memory_cells'],
                'AGR2': ['T_memory_cells', 'Endothelial_cells'],
                'CDC20B': ['Ependymal_cells'],
                'MT1G': ['T_memory_cells'],
                'COL11A1': ['Hepatic_stellate_cells'],
                'COL10A1': ['Hepatic_stellate_cells'],
                'SERPING1': ['Dendritic_cells'],
                'GBP1': ['Dendritic_cells'],
                'SFTPC': ['Pulmonary_alveolar_type_II_cells', 'Myoepithelial_cells'],
                'MALAT1': ['Myoepithelial_cells']
            }
            
            # Only plot for cell types mentioned in the report for this gene
            if gene in key_cell_types:
                for cell_type in key_cell_types[gene]:
                    # Check if this cell type exists in our data
                    if cell_type not in adata.obs['predicted_cell_type'].unique():
                        continue
                    
                    # Filter to this cell type
                    ct_mask = adata.obs['predicted_cell_type'] == cell_type
                    if np.sum(ct_mask) < min_cells:
                        continue
                    
                    ct_adata = adata[ct_mask].copy()
                    
                    # Create violin plot for this cell type
                    plt.figure(figsize=(10, 6))
                    sc.pl.violin(ct_adata, [gene], groupby='Sex_Tissue', show=False)
                    plt.title(f"{gene} in {cell_type.replace('_', ' ')}", fontsize=14)
                    plt.tight_layout()
                    plt.savefig(os.path.join(output_dir, f"{gene}_{cell_type}.pdf"), dpi=300)
                    plt.close()
    
    print(f"All violin plots saved to {output_dir}")
    return

def main():
    parser = argparse.ArgumentParser(description='Generate violin plots for sex-specific marker genes')
    parser.add_argument('-i', '--input', required=True, help='Input h5ad file')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    parser.add_argument('-m', '--min_cells', type=int, default=10, help='Minimum cells required for plotting')
    args = parser.parse_args()
    
    # Load data
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells with {adata.n_vars} genes")
    
    # Generate marker violin plots
    generate_marker_violins(adata, args.output_dir, args.min_cells)

if __name__ == "__main__":
    main()
