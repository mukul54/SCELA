"""
Functions for demographic comparisons of cell type distributions
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from .config import VIS_DIR

def compare_demographics(adata):
    """
    Compare cell type composition between different demographic groups
    """
    print("\n=== Comparing Demographics ===")
    
    # Determine which cell type column to use
    cell_type_column = 'cell_type' if 'cell_type' in adata.obs.columns else 'predicted_cell_type'
    
    # Ensure cell_type_column is categorical
    if not pd.api.types.is_categorical_dtype(adata.obs[cell_type_column]):
        adata.obs[cell_type_column] = adata.obs[cell_type_column].astype('category')
        
    # Create a display-friendly version of the column name for plot titles
    cell_type_display = 'Cell Type' if cell_type_column == 'cell_type' else 'Predicted Cell Type'
    if cell_type_column not in adata.obs.columns:
        print(f"Error: No {cell_type_column} column found in the data")
        return adata
    
    # Create directory for visualizations
    os.makedirs(VIS_DIR, exist_ok=True)
    
    # 1. Male vs. Female comparison
    if 'Sex' in adata.obs.columns:
        # Clean data - handle missing values
        if adata.obs['Sex'].isna().any():
            print(f"Warning: {adata.obs['Sex'].isna().sum()} cells have missing sex data. These will be excluded.")
        
        valid_sex = adata.obs['Sex'].dropna()
        sex_groups = valid_sex.value_counts()
        print(f"Sex groups: {dict(sex_groups)}")
        
        # Only proceed if we have multiple groups
        if len(sex_groups) > 1:
            # Get cell type proportions per sex
            cell_counts = pd.crosstab(
                adata.obs[cell_type_column], 
                adata.obs['Sex'],
                normalize='columns'
            ) * 100  # Convert to percentage
            
            # Plot proportions
            plt.figure(figsize=(14, 9))
            cell_counts.plot(kind='bar', stacked=False)
            plt.title(f'{cell_type_display} Distribution: Male vs. Female', fontsize=14, pad=20)
            plt.ylabel('Percentage of Cells', fontsize=12, labelpad=10)
            plt.xticks(rotation=45, ha='right', fontsize=10)
            plt.yticks(fontsize=10)
            plt.legend(title='Sex', fontsize=10)
            plt.grid(axis='y', linestyle='--', alpha=0.7)
            plt.tight_layout(pad=3.0)
            
            # Save the figure
            sex_comparison_file = f"{VIS_DIR}/sex_comparison.pdf"
            plt.savefig(sex_comparison_file, bbox_inches='tight', dpi=150)
            plt.close()
            print(f"Saved sex comparison to {sex_comparison_file}")
    
    # 2. Smoking status comparison
    if 'Smoking' in adata.obs.columns:
        # Clean data
        if adata.obs['Smoking'].isna().any():
            print(f"Warning: {adata.obs['Smoking'].isna().sum()} cells have missing smoking data. These will be excluded.")
        
        valid_smoking = adata.obs['Smoking'].dropna()
        smoking_groups = valid_smoking.value_counts()
        print(f"Smoking groups: {dict(smoking_groups)}")
        
        # Only proceed if we have multiple groups
        if len(smoking_groups) > 1:
            # Get cell type proportions per smoking status
            smoking_counts = pd.crosstab(
                adata.obs[cell_type_column], 
                adata.obs['Smoking'],
                normalize='columns'
            ) * 100  # Convert to percentage
            
            # Plot proportions
            plt.figure(figsize=(14, 9))
            smoking_counts.plot(kind='bar', stacked=False)
            plt.title(f'{cell_type_display} Distribution by Smoking Status', fontsize=14, pad=20)
            plt.ylabel('Percentage of Cells', fontsize=12, labelpad=10)
            plt.xticks(rotation=45, ha='right', fontsize=10)
            plt.yticks(fontsize=10)
            plt.legend(title='Smoking Status', fontsize=10)
            plt.grid(axis='y', linestyle='--', alpha=0.7)
            plt.tight_layout(pad=3.0)
            
            # Save the figure
            smoking_comparison_file = f"{VIS_DIR}/smoking_comparison.pdf"
            plt.savefig(smoking_comparison_file, bbox_inches='tight', dpi=150)
            plt.close()
            print(f"Saved smoking comparison to {smoking_comparison_file}")
    
    # 3. Tumor vs. Normal comparison
    if 'Sample_Origin' in adata.obs.columns or 'tissue_type' in adata.obs.columns:
        # Determine which column to use
        origin_column = 'Sample_Origin' if 'Sample_Origin' in adata.obs.columns else 'tissue_type'
        
        # Clean data
        if adata.obs[origin_column].isna().any():
            print(f"Warning: {adata.obs[origin_column].isna().sum()} cells have missing tissue origin data. These will be excluded.")
        
        valid_origin = adata.obs[origin_column].dropna()
        origin_groups = valid_origin.value_counts()
        print(f"Tissue origin groups: {dict(origin_groups)}")
        
        # Only proceed if we have multiple groups
        if len(origin_groups) > 1:
            # Get cell type proportions per tissue origin
            origin_counts = pd.crosstab(
                adata.obs[cell_type_column], 
                adata.obs[origin_column],
                normalize='columns'
            ) * 100  # Convert to percentage
            
            # Plot proportions
            plt.figure(figsize=(14, 9))
            origin_counts.plot(kind='bar', stacked=False)
            plt.title(f'{cell_type_display} Distribution: Tumor vs. Normal Tissue', fontsize=14, pad=20)
            plt.ylabel('Percentage of Cells', fontsize=12, labelpad=10)
            plt.xticks(rotation=45, ha='right', fontsize=10)
            plt.yticks(fontsize=10)
            plt.legend(title='Tissue Origin', fontsize=10)
            plt.grid(axis='y', linestyle='--', alpha=0.7)
            plt.tight_layout(pad=3.0)
            
            # Save the figure
            origin_comparison_file = f"{VIS_DIR}/tumor_normal_comparison.pdf"
            plt.savefig(origin_comparison_file, bbox_inches='tight', dpi=150)
            plt.close()
            print(f"Saved tumor vs. normal comparison to {origin_comparison_file}")
    
    return adata
