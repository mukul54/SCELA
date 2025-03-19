"""
Functions for demographic comparisons of cell type distributions
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from .config import VIS_DIR, DYNAMIC_FIGURE_SIZE

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
    
    # Get all possible demographic columns
    demographic_columns = {
        'sex': next((col for col in adata.obs.columns if col.lower() in ['sex', 'gender']), None),
        'smoking': next((col for col in adata.obs.columns if 'smok' in col.lower()), None),
        'age': next((col for col in adata.obs.columns if col.lower() in ['age']), None)
    }
    
    print(f"Found demographic columns: {demographic_columns}")
    
    # Ensure all demographics are in proper format
    for demo_type, column in demographic_columns.items():
        if column and column in adata.obs.columns:
            # Handle missing values
            if adata.obs[column].isna().any():
                missing_count = adata.obs[column].isna().sum()
                print(f"Warning: {missing_count} cells have missing {demo_type} data. These will be excluded from analysis.")
            
            # Convert to categorical if not already
            if not pd.api.types.is_categorical_dtype(adata.obs[column]):
                adata.obs[column] = adata.obs[column].astype('category')
    
    # 1. Sex/Gender comparison
    sex_col = demographic_columns['sex']
    if sex_col and sex_col in adata.obs.columns:
        valid_sex = adata.obs[sex_col].dropna()
        sex_groups = valid_sex.value_counts()
        print(f"Sex/Gender groups: {dict(sex_groups)}")
        
        # Only proceed if we have multiple groups
        if len(sex_groups) > 1:
            # Get cell type proportions per sex/gender
            cell_counts = pd.crosstab(
                adata.obs[cell_type_column], 
                adata.obs[sex_col],
                normalize='columns'
            ) * 100  # Convert to percentage
            
            # Dynamic figure sizing if enabled
            if DYNAMIC_FIGURE_SIZE:
                n_cell_types = len(adata.obs[cell_type_column].unique())
                fig_width = max(12, min(18, len(sex_groups) * 2.5))
                fig_height = max(8, min(16, n_cell_types * 0.6))
                plt.figure(figsize=(fig_width, fig_height))
            else:
                plt.figure(figsize=(14, 9))
            
            # Plot proportions
            cell_counts.plot(kind='bar', stacked=False)
            plt.title(f'{cell_type_display} Distribution by {sex_col}', fontsize=14, pad=20)
            plt.ylabel('Percentage of Cells', fontsize=12, labelpad=10)
            plt.xticks(rotation=45, ha='right', fontsize=10)
            plt.yticks(fontsize=10)
            
            # Determine optimal number of legend columns based on number of groups
            n_legend_cols = max(1, min(len(sex_groups), 5))
            plt.legend(title=sex_col, fontsize=10, ncol=n_legend_cols, bbox_to_anchor=(1.02, 1), loc='upper left')
            
            # Adjust figure size if needed to accommodate legend
            if len(sex_groups) > 5:
                fig = plt.gcf()
                fig_size = fig.get_size_inches()
                fig.set_size_inches(fig_size[0] + 2, fig_size[1])
            
            plt.grid(axis='y', linestyle='--', alpha=0.7)
            plt.tight_layout(pad=3.0)
            
            # Save the figure
            filename = sex_col.lower().replace('/', '_').replace(' ', '_')
            sex_comparison_file = f"{VIS_DIR}/{filename}_comparison.pdf"
            plt.savefig(sex_comparison_file, bbox_inches='tight', dpi=150)
            plt.close()
            print(f"Saved {sex_col} comparison to {sex_comparison_file}")
            
            # Create an improved heatmap visualization without numerical annotations
            plt.figure(figsize=(12, max(8, n_cell_types * 0.5)))
            
            # Use better colormap and vertical colorbar
            cbar_kws = {
                'orientation': 'vertical',
                'label': 'Percentage (%)',
                'shrink': 0.8,
                'pad': 0.02,
                'aspect': 20
            }
            
            # Create heatmap without annotations
            ax = sns.heatmap(cell_counts, annot=False, cmap='viridis', 
                       cbar_kws=cbar_kws, robust=True)
                       
            # Improve axis labels
            plt.xlabel(sex_col, fontsize=12, labelpad=10)
            plt.ylabel(cell_type_display, fontsize=12, labelpad=10)
            
            # Format the ticks
            plt.yticks(rotation=0, fontsize=10)
            plt.xticks(rotation=45, ha='right', fontsize=10)
            
            plt.title(f'{cell_type_display} Distribution by {sex_col}', fontsize=14)
            plt.tight_layout()
            
            heatmap_file = f"{VIS_DIR}/{filename}_heatmap.pdf"
            plt.savefig(heatmap_file, bbox_inches='tight', dpi=150)
            plt.close()
            print(f"Saved improved {sex_col} heatmap to {heatmap_file}")
    
    # 2. Smoking status comparison
    smoking_col = demographic_columns['smoking']
    if smoking_col and smoking_col in adata.obs.columns:
        valid_smoking = adata.obs[smoking_col].dropna()
        smoking_groups = valid_smoking.value_counts()
        print(f"Smoking groups: {dict(smoking_groups)}")
        
        # Only proceed if we have multiple groups
        if len(smoking_groups) > 1:
            # Get cell type proportions per smoking status
            smoking_counts = pd.crosstab(
                adata.obs[cell_type_column], 
                adata.obs[smoking_col],
                normalize='columns'
            ) * 100  # Convert to percentage
            
            # Dynamic figure sizing if enabled
            if DYNAMIC_FIGURE_SIZE:
                n_cell_types = len(adata.obs[cell_type_column].unique())
                fig_width = max(12, min(18, len(smoking_groups) * 2.5))
                fig_height = max(8, min(16, n_cell_types * 0.6))
                plt.figure(figsize=(fig_width, fig_height))
            else:
                plt.figure(figsize=(14, 9))
            
            # Plot proportions
            smoking_counts.plot(kind='bar', stacked=False)
            plt.title(f'{cell_type_display} Distribution by {smoking_col}', fontsize=14, pad=20)
            plt.ylabel('Percentage of Cells', fontsize=12, labelpad=10)
            plt.xticks(rotation=45, ha='right', fontsize=10)
            plt.yticks(fontsize=10)
            
            # Determine optimal number of legend columns based on number of groups
            n_legend_cols = max(1, min(len(smoking_groups), 5))
            plt.legend(title=smoking_col, fontsize=10, ncol=n_legend_cols, bbox_to_anchor=(1.02, 1), loc='upper left')
            
            # Adjust figure size if needed to accommodate legend
            if len(smoking_groups) > 5:
                fig = plt.gcf()
                fig_size = fig.get_size_inches()
                fig.set_size_inches(fig_size[0] + 2, fig_size[1])
            
            plt.grid(axis='y', linestyle='--', alpha=0.7)
            plt.tight_layout(pad=3.0)
            
            # Save the figure
            filename = smoking_col.lower().replace('/', '_').replace(' ', '_')
            smoking_comparison_file = f"{VIS_DIR}/{filename}_comparison.pdf"
            plt.savefig(smoking_comparison_file, bbox_inches='tight', dpi=150)
            plt.close()
            print(f"Saved {smoking_col} comparison to {smoking_comparison_file}")
            
            # Create an improved heatmap visualization without numerical annotations
            plt.figure(figsize=(12, max(8, n_cell_types * 0.5)))
            
            # Use better colormap and vertical colorbar
            cbar_kws = {
                'orientation': 'vertical',
                'label': 'Percentage (%)',
                'shrink': 0.8,
                'pad': 0.02,
                'aspect': 20
            }
            
            # Create heatmap without annotations
            ax = sns.heatmap(smoking_counts, annot=False, cmap='viridis', 
                       cbar_kws=cbar_kws, robust=True)
                       
            # Improve axis labels
            plt.xlabel(smoking_col, fontsize=12, labelpad=10)
            plt.ylabel(cell_type_display, fontsize=12, labelpad=10)
            
            # Format the ticks
            plt.yticks(rotation=0, fontsize=10)
            plt.xticks(rotation=45, ha='right', fontsize=10)
            
            plt.title(f'{cell_type_display} Distribution by {smoking_col}', fontsize=14)
            plt.tight_layout()
            
            heatmap_file = f"{VIS_DIR}/{filename}_heatmap.pdf"
            plt.savefig(heatmap_file, bbox_inches='tight', dpi=150)
            plt.close()
            print(f"Saved improved {smoking_col} heatmap to {heatmap_file}")
    
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
            
            # Dynamic figure sizing if enabled
            if DYNAMIC_FIGURE_SIZE:
                n_cell_types = len(adata.obs[cell_type_column].unique())
                fig_width = max(12, min(18, len(origin_groups) * 2.5))
                fig_height = max(8, min(16, n_cell_types * 0.6))
                plt.figure(figsize=(fig_width, fig_height))
            else:
                plt.figure(figsize=(14, 9))
                
            # Plot proportions
            origin_counts.plot(kind='bar', stacked=False)
            plt.title(f'{cell_type_display} Distribution: Tumor vs. Normal Tissue', fontsize=14, pad=20)
            plt.ylabel('Percentage of Cells', fontsize=12, labelpad=10)
            plt.xticks(rotation=45, ha='right', fontsize=10)
            plt.yticks(fontsize=10)
            
            # Determine optimal number of legend columns based on number of groups
            n_legend_cols = max(1, min(len(origin_groups), 5))
            plt.legend(title='Tissue Origin', fontsize=10, ncol=n_legend_cols, bbox_to_anchor=(1.02, 1), loc='upper left')
            
            # Adjust figure size if needed to accommodate legend
            if len(origin_groups) > 5:
                fig = plt.gcf()
                fig_size = fig.get_size_inches()
                fig.set_size_inches(fig_size[0] + 2, fig_size[1])
            
            plt.grid(axis='y', linestyle='--', alpha=0.7)
            plt.tight_layout(pad=3.0)
            
            # Save the figure
            origin_comparison_file = f"{VIS_DIR}/tumor_normal_comparison.pdf"
            plt.savefig(origin_comparison_file, bbox_inches='tight', dpi=150)
            plt.close()
            print(f"Saved tumor vs. normal comparison to {origin_comparison_file}")
            
            # Create an improved heatmap visualization without numerical annotations
            plt.figure(figsize=(12, max(8, n_cell_types * 0.5)))
            
            # Use better colormap and vertical colorbar
            cbar_kws = {
                'orientation': 'vertical',
                'label': 'Percentage (%)',
                'shrink': 0.8,
                'pad': 0.02,
                'aspect': 20
            }
            
            # Create heatmap without annotations (annot=False)
            ax = sns.heatmap(origin_counts, annot=False, cmap='viridis', 
                       cbar_kws=cbar_kws, robust=True)
                       
            # Improve axis labels
            plt.xlabel('Sample Origin', fontsize=12, labelpad=10)
            plt.ylabel(cell_type_display, fontsize=12, labelpad=10)
            
            # Format the ticks
            plt.yticks(rotation=0, fontsize=10)
            plt.xticks(rotation=45, ha='right', fontsize=10)
            
            plt.title(f'{cell_type_display} Distribution: Tumor vs. Normal (Heatmap)', fontsize=14)
            plt.tight_layout()
            
            heatmap_file = f"{VIS_DIR}/tumor_normal_heatmap.pdf"
            plt.savefig(heatmap_file, bbox_inches='tight', dpi=150)
            plt.close()
            print(f"Saved improved tumor/normal heatmap to {heatmap_file}")
    
    # 4. Age analysis
    age_col = demographic_columns['age']
    if age_col and age_col in adata.obs.columns:
        # Check if we have sufficient non-missing age data
        if adata.obs[age_col].isna().sum() / len(adata.obs) > 0.5:
            print(f"Warning: More than 50% of cells have missing age data. Skipping age analysis.")
        else:
            print("Analyzing cell type distribution by age...")
            
            # Create age groups for better analysis
            valid_age = adata.obs[adata.obs[age_col].notna()]
            
            # Determine age bins based on data distribution
            age_min = valid_age[age_col].min()
            age_max = valid_age[age_col].max()
            
            # Create age groups with approximately equal bin sizes
            n_bins = min(5, max(3, (age_max - age_min) // 10))
            age_bins = np.linspace(age_min, age_max, n_bins + 1)
            age_labels = [f'{int(age_bins[i])}-{int(age_bins[i+1])}' for i in range(len(age_bins)-1)]
            
            # Add age group column
            adata.obs['age_group'] = pd.cut(adata.obs[age_col], bins=age_bins, labels=age_labels, include_lowest=True)
            
            # Count by age group
            age_groups = adata.obs['age_group'].value_counts()
            print(f"Age groups: {dict(age_groups)}")
            
            if len(age_groups) > 1:
                # Calculate cell type distribution by age group
                age_counts = pd.crosstab(
                    adata.obs[cell_type_column],
                    adata.obs['age_group'],
                    normalize='columns'
                ) * 100
                
                # Plot age group distribution
                if DYNAMIC_FIGURE_SIZE:
                    n_cell_types = len(adata.obs[cell_type_column].unique())
                    fig_width = max(12, min(18, len(age_groups) * 2.5))
                    fig_height = max(8, min(16, n_cell_types * 0.6))
                    plt.figure(figsize=(fig_width, fig_height))
                else:
                    plt.figure(figsize=(14, 9))
                
                age_counts.plot(kind='bar', stacked=False)
                plt.title(f'{cell_type_display} Distribution by Age Group', fontsize=14, pad=20)
                plt.ylabel('Percentage of Cells', fontsize=12, labelpad=10)
                plt.xticks(rotation=45, ha='right', fontsize=10)
                plt.yticks(fontsize=10)
                
                # Determine optimal number of legend columns based on number of groups
                n_legend_cols = max(1, min(len(age_groups), 5))
                plt.legend(title='Age Group', fontsize=10, ncol=n_legend_cols, bbox_to_anchor=(1.02, 1), loc='upper left')
                
                # Adjust figure size if needed to accommodate legend
                if len(age_groups) > 5:
                    fig = plt.gcf()
                    fig_size = fig.get_size_inches()
                    fig.set_size_inches(fig_size[0] + 2, fig_size[1])
                
                plt.grid(axis='y', linestyle='--', alpha=0.7)
                plt.tight_layout(pad=3.0)
                
                # Save the figure
                age_comparison_file = f"{VIS_DIR}/age_group_comparison.pdf"
                plt.savefig(age_comparison_file, bbox_inches='tight', dpi=150)
                plt.close()
                print(f"Saved age group comparison to {age_comparison_file}")
                
                # Create an improved heatmap visualization without numerical annotations
                plt.figure(figsize=(12, max(8, n_cell_types * 0.5)))
                
                # Use better colormap and vertical colorbar
                cbar_kws = {
                    'orientation': 'vertical',
                    'label': 'Percentage (%)',
                    'shrink': 0.8,
                    'pad': 0.02,
                    'aspect': 20
                }
                
                # Create heatmap without annotations
                ax = sns.heatmap(age_counts, annot=False, cmap='viridis', 
                           cbar_kws=cbar_kws, robust=True)
                           
                # Improve axis labels
                plt.xlabel('Age Group', fontsize=12, labelpad=10)
                plt.ylabel(cell_type_display, fontsize=12, labelpad=10)
                
                # Format the ticks
                plt.yticks(rotation=0, fontsize=10)
                plt.xticks(rotation=45, ha='right', fontsize=10)
                
                plt.title(f'{cell_type_display} Distribution by Age Group', fontsize=14)
                plt.tight_layout()
                
                heatmap_file = f"{VIS_DIR}/age_group_heatmap.pdf"
                plt.savefig(heatmap_file, bbox_inches='tight', dpi=150)
                plt.close()
                print(f"Saved improved age group heatmap to {heatmap_file}")
    
    return adata
