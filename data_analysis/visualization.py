"""
Functions for generating visualizations for cell type analysis
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from .config import VIS_DIR, DYNAMIC_FIGURE_SIZE

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
        
        # Determine dynamic figure size if enabled
        n_cell_types = len(adata.obs[cell_type_column].unique())
        n_samples = len(adata.obs['sample_id'].unique())
        
        if DYNAMIC_FIGURE_SIZE:
            # Adjust figure width based on number of samples
            fig_width = max(14, min(20, n_samples * 1.5))
            # Adjust figure height based on number of cell types (for legend)
            fig_height = max(8, min(16, n_cell_types * 0.5))
            plt.figure(figsize=(fig_width, fig_height))
        else:
            plt.figure(figsize=(16, 10))
        
        # Generate cross-tabulation - each column (sample) should sum to 100%
        sample_dist = pd.crosstab(
            adata.obs[cell_type_column], 
            adata.obs['sample_id'],
            normalize='columns'
        ) * 100  # Convert to percentage
        
        # Verify percentages sum to 100 for each sample
        for col in sample_dist.columns:
            if not (99.0 <= sample_dist[col].sum() <= 101.0):
                print(f"Warning: Percentages for sample {col} sum to {sample_dist[col].sum():.2f}%, not 100%")
                # Normalize to exactly 100% to fix any floating point issues
                sample_dist[col] = sample_dist[col] / sample_dist[col].sum() * 100
        
        # Plot a properly stacked bar chart where each sample adds up to exactly 100%
        ax = sample_dist.T.plot(kind='bar', stacked=True, figsize=(16, 10), colormap='viridis')
        
        # Set proper title and labels
        plt.title('Cell Type Distribution Within Each Sample', fontsize=16, pad=20)
        plt.xlabel('Sample ID', fontsize=14, labelpad=10)
        plt.ylabel('Percentage of Cells (%)', fontsize=14, labelpad=10)
        
        # Improve tick readability
        plt.xticks(rotation=45, ha='right', fontsize=10)
        plt.yticks(fontsize=10)
        
        # Add horizontal grid lines
        plt.grid(axis='y', linestyle='--', alpha=0.6)
        
        # Ensure y-axis goes from 0 to 100%
        plt.ylim(0, 100)
        
        # Add a horizontal line at 100% for visual reference
        plt.axhline(y=100, color='black', linestyle='-', alpha=0.3, linewidth=1)
        
        # Add percentage annotations at important thresholds
        for y in [25, 50, 75]:
            plt.axhline(y=y, color='black', linestyle=':', alpha=0.2, linewidth=0.5)
            plt.text(len(sample_dist.columns), y, f'{y}%', va='center', ha='left', fontsize=8, alpha=0.7)
        
        # Determine optimal number of legend columns based on number of categories
        # More categories = more columns, but not too many to make it unreadable
        n_categories = len(adata.obs[cell_type_column].unique())
        n_legend_cols = max(1, min(n_categories // 5, 5))
        
        # Position legend outside plot to avoid overlap
        # The bbox_to_anchor ensures it's placed to the right of the plot
        plt.legend(title=cell_type_display, bbox_to_anchor=(1.02, 1), 
                  loc='upper left', ncol=n_legend_cols)
        
        # Adjust figure size to accommodate the legend if needed
        if n_categories > 15:
            # Get the current figure
            fig = plt.gcf()
            # Get the current size
            fig_size = fig.get_size_inches()
            # Add extra width for the legend
            fig.set_size_inches(fig_size[0] + 3, fig_size[1])
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
    # Dynamic figure size based on number of cell types
    if DYNAMIC_FIGURE_SIZE:
        fig_width = max(10, min(20, len(cell_counts) * 0.8))
        plt.figure(figsize=(fig_width, 8))
    else:
        plt.figure(figsize=(10, 6))
        
    # Sort by count for better visualization
    cell_counts_sorted = cell_counts.sort_values('Count', ascending=False)
    
    # Use a horizontal bar chart for better readability with many cell types
    if len(cell_counts) > 10:
        ax = sns.barplot(y='Cell Type', x='Count', data=cell_counts_sorted)
        # Add count labels to the bars
        for i, p in enumerate(ax.patches):
            width = p.get_width()
            ax.text(width + width*0.02, p.get_y() + p.get_height()/2,
                    f'{width:.0f}', ha='left', va='center')
        plt.xlabel('Count', fontsize=12)
        plt.ylabel('Cell Type', fontsize=12)
    else:
        ax = sns.barplot(x='Cell Type', y='Count', data=cell_counts_sorted)
        # Add count labels to the bars
        for i, p in enumerate(ax.patches):
            height = p.get_height()
            ax.text(p.get_x() + p.get_width()/2, height + height*0.02,
                    f'{height:.0f}', ha='center')
        plt.xlabel('Cell Type', fontsize=12)
        plt.ylabel('Count', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        
    plt.title('Cell Type Counts', fontsize=14)
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
