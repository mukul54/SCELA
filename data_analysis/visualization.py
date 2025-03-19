"""
Functions for generating visualizations for cell type analysis
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from .config import VIS_DIR

def generate_visualizations(adata):
    """
    Create comprehensive cell type distribution visualizations
    """
    print("\n=== Generating Visualizations ===")
    
    # Make sure the visualization directory exists
    os.makedirs(VIS_DIR, exist_ok=True)
    
    # Cell type column to use
    cell_type_column = 'cell_type' if 'cell_type' in adata.obs.columns else 'predicted_cell_type'
    
    # Ensure cell_type_column is categorical
    if not pd.api.types.is_categorical_dtype(adata.obs[cell_type_column]):
        adata.obs[cell_type_column] = adata.obs[cell_type_column].astype('category')
        
    # Create a display-friendly version of the column name for plot titles
    cell_type_display = 'Cell Type' if cell_type_column == 'cell_type' else 'Predicted Cell Type'
    
    # 1. Distribution by sample
    if 'sample_id' in adata.obs.columns:
        print("Generating sample distribution plot...")
        plt.figure(figsize=(14, 8))
        
        # Generate cross-tabulation
        sample_dist = pd.crosstab(
            adata.obs[cell_type_column], 
            adata.obs['sample_id'],
            normalize='columns'
        ) * 100  # Convert to percentage
        
        # Plot stacked bar chart
        sample_dist.plot(kind='bar', stacked=True, colormap='viridis')
        plt.title('Cell Type Distribution Across Samples', fontsize=14, pad=20)
        plt.xlabel('Sample ID', fontsize=12)
        plt.ylabel('Percentage', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.legend(title=cell_type_display, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        
        # Save figure
        sample_dist_file = f"{VIS_DIR}/sample_distribution.pdf"
        plt.savefig(sample_dist_file, bbox_inches='tight', dpi=150)
        plt.close()
        print(f"Saved sample distribution to {sample_dist_file}")
    
    # 2. Cell type UMAP with better aesthetics
    if 'X_umap' in adata.obsm:
        print("Generating enhanced UMAP visualization...")
        
        # Format category names for display
        if cell_type_column in adata.obs:
            # Ensure the column is categorical
            if not pd.api.types.is_categorical_dtype(adata.obs[cell_type_column]):
                adata.obs[cell_type_column] = adata.obs[cell_type_column].astype('category')
                
            # Create formatted version of categories
            categories = adata.obs[cell_type_column].cat.categories
            formatted_categories = [cat.replace('_', ' ').replace('-', ' ').title() for cat in categories]
            
            # Create a temporary formatted column
            temp_col = f"{cell_type_column}_formatted"
            category_map = dict(zip(categories, formatted_categories))
            adata.obs[temp_col] = adata.obs[cell_type_column].map(category_map)
            adata.obs[temp_col] = adata.obs[temp_col].astype('category')
            
            # Create publication-quality UMAP with better cell type separation
            plt.figure(figsize=(12, 10))
            # Set figure title font properties before plotting
            plt.rcParams['axes.titlesize'] = 14
            
            sc.pl.umap(
                adata,
                color=temp_col,
                legend_loc='on data',
                legend_fontsize=10,
                legend_fontoutline=2,
                title='Cell Types in Lung Cancer Dataset',
                frameon=False,
                s=50,  # Larger point size
                alpha=0.7,
                palette='tab20',  # Better color palette for distinction
                show=False
            )
            
            # Reset title size to default after plotting
            plt.rcParams['axes.titlesize'] = plt.rcParamsDefault['axes.titlesize']
            
            # Get the current figure and axes
            fig = plt.gcf()
            ax = plt.gca()
            
            # Add a light grid to help with spatial reference
            ax.grid(True, linestyle='--', alpha=0.3, linewidth=0.5)
            
            # Clean up temporary column
            del adata.obs[temp_col]
            
            # Add annotations
            plt.xlabel("UMAP Dimension 1", fontsize=12, labelpad=10)
            plt.ylabel("UMAP Dimension 2", fontsize=12, labelpad=10)
            
            # Save high-quality figure
            enhanced_umap_file = f"{VIS_DIR}/enhanced_cell_types_umap.pdf"
            plt.savefig(enhanced_umap_file, bbox_inches='tight', dpi=300)
            plt.close()
            print(f"Saved enhanced UMAP to {enhanced_umap_file}")
    
    # 3. Create cell count summary table
    print("Generating cell count summary...")
    cell_counts = adata.obs[cell_type_column].value_counts().reset_index()
    cell_counts.columns = ['Cell Type', 'Count']
    cell_counts['Percentage'] = (cell_counts['Count'] / cell_counts['Count'].sum() * 100).round(2)
    
    # Generate a simple visualization of the counts
    plt.figure(figsize=(10, 6))
    sns.barplot(x='Cell Type', y='Count', data=cell_counts)
    plt.title('Cell Type Counts', fontsize=14)
    plt.xlabel('Cell Type', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    # Save the figure
    count_summary_file = f"{VIS_DIR}/cell_type_counts.pdf"
    plt.savefig(count_summary_file, bbox_inches='tight', dpi=150)
    plt.close()
    
    # Save the summary table as CSV
    count_table_file = f"{VIS_DIR}/cell_type_counts.csv"
    cell_counts.to_csv(count_table_file, index=False)
    print(f"Saved cell count summary to {count_summary_file} and {count_table_file}")
    
    return adata
