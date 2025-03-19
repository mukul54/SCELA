#!/usr/bin/env python
# Cell Type Annotation for GSE131907 Lung Cancer Dataset

import os
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse
import warnings
warnings.filterwarnings('ignore')

# File Paths
BASE_DATA_PATH = "/l/users/mukul.ranjan/cancer/data/GSE131907"
ANNDATA_FILE = f"{BASE_DATA_PATH}/GSE131907_anndata.h5ad"
METADATA_FILE = f"{BASE_DATA_PATH}/GSE131907_metadata.csv"
RESULTS_DIR = f"{BASE_DATA_PATH}/analysis_results"
VIS_DIR = f"{RESULTS_DIR}/visualizations"
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(VIS_DIR, exist_ok=True)

# Set plotting defaults
sc.settings.set_figure_params(dpi=100, frameon=False)
plt.rcParams['figure.figsize'] = (8, 6)

# Configure scanpy to save figures directly (without prefixing)
import matplotlib as mpl
sc.settings.figdir = './'  # Set base figure directory to current directory

# Create a function to safely save plots
def save_fig(fig, filename):
    """Save a figure to a specific path, creating directories if needed"""
    path = os.path.dirname(filename)
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
    fig.savefig(filename, bbox_inches='tight', dpi=150)
    
# Configure figure appearance to reduce text overlap
plt.rcParams['figure.autolayout'] = True  # Better automatic layout
plt.rcParams['axes.labelpad'] = 10  # Add padding to axis labels
plt.rcParams['xtick.major.pad'] = 5  # Padding for x-tick labels
plt.rcParams['ytick.major.pad'] = 5  # Padding for y-tick labels

def load_data():
    """
    Load the pre-processed AnnData object and metadata
    """
    print("\n=== Loading Data ===")
    
    if not os.path.exists(ANNDATA_FILE):
        print(f"Error: AnnData file not found at {ANNDATA_FILE}")
        print("Please run the data processing script first.")
        return None
    
    print(f"Loading AnnData from {ANNDATA_FILE}...")
    adata = sc.read_h5ad(ANNDATA_FILE)
    
    # Load and fix metadata
    print("Loading and integrating metadata...")
    metadata = pd.read_csv(METADATA_FILE)
    
    # Create a mapping from sample ID to metadata
    sample_to_metadata = {}
    for _, row in metadata.iterrows():
        sample_id = row['Samples']  # Using 'Samples' column, not 'sample_id'
        sample_to_metadata[sample_id] = row.to_dict()
    
    # Add metadata to AnnData object
    for cell_id in adata.obs_names:
        # Extract sample ID from cell barcode
        match = None
        for sample in sample_to_metadata.keys():
            if sample in cell_id:
                match = sample
                break
        
        if match:
            meta = sample_to_metadata[match]
            # Add metadata fields to observation
            for key, value in meta.items():
                if key != 'Samples':  # Skip the sample ID field
                    adata.obs.loc[cell_id, key] = value
    
    # Create categorical variables for grouping
    if 'Sex' in adata.obs.columns:
        adata.obs['Sex'] = adata.obs['Sex'].astype('category')
    
    if 'Smoking' in adata.obs.columns:
        adata.obs['Smoking'] = adata.obs['Smoking'].astype('category')
    
    # Create tissue type variable (Tumor vs Normal)
    adata.obs['tissue_type'] = adata.obs['sample_id'].apply(
        lambda x: 'Tumor' if ('LUNG_T' in str(x) or 'EBUS' in str(x) or 
                              'BRONCHO' in str(x) or 'EFFUSION' in str(x)) 
                  else 'Normal' if 'LUNG_N' in str(x) 
                  else 'Other'
    )
    adata.obs['tissue_type'] = adata.obs['tissue_type'].astype('category')
    
    print(f"Dataset dimensions: {adata.shape[0]} cells × {adata.shape[1]} genes")
    print(f"Available metadata: {', '.join(adata.obs.columns)}")
    
    return adata

def identify_marker_genes(adata):
    """
    Identify marker genes for each cluster
    """
    print("\n=== Identifying Marker Genes ===")
    
    # Make sure we have clustering information
    if 'leiden' not in adata.obs.columns:
        print("Performing clustering...")
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.leiden(adata, resolution=0.5)
    
    # Calculate UMAP if not already done
    if 'X_umap' not in adata.obsm:
        print("Computing UMAP embedding...")
        sc.tl.umap(adata)
        
    # Apply log normalization if not already done (to prevent warnings)
    if not adata.uns.get('log1p', {}).get('base', None):
        print("Logarithmizing the data before finding marker genes...")
        sc.pp.log1p(adata)
    
    # Find markers for each cluster
    print("Finding marker genes for each cluster...")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    
    # Extract markers and save to file
    print("Extracting top markers...")
    markers_df = pd.DataFrame()
    
    for cluster in adata.obs['leiden'].unique():
        cluster_markers = sc.get.rank_genes_groups_df(adata, group=cluster)
        top_markers = cluster_markers.head(10)  # Top 10 markers per cluster
        markers_df = pd.concat([markers_df, top_markers])
    
    # Save markers to file
    markers_file = f"{RESULTS_DIR}/cluster_markers.csv"
    markers_df.to_csv(markers_file)
    print(f"Saved cluster markers to {markers_file}")
    
    # Create marker gene heatmap
    plt.figure(figsize=(14, 10))  # Increased figure size
    ax = sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, show_gene_labels=True, show=False,
                                        gene_symbols=None)  # Use gene IDs directly if symbols not available
    
    # Adjust gene labels to prevent overlap
    plt.tight_layout(pad=3.0)  # Increase padding
    
    # Get the axes and adjust y-axis label properties
    for ax in plt.gcf().get_axes():
        for tick in ax.get_yticklabels():
            tick.set_fontsize(8)  # Smaller font for y-axis labels
    
    # Save the figure manually to the correct path
    marker_heatmap_file = f"{VIS_DIR}/marker_genes_heatmap.pdf"
    os.makedirs(VIS_DIR, exist_ok=True)  # Ensure directory exists
    plt.savefig(marker_heatmap_file, bbox_inches='tight', dpi=150)
    plt.close()
    print(f"Saved marker heatmap to {marker_heatmap_file}")
    
    return adata

def annotate_cell_types(adata):
    """
    Annotate cell types based on marker genes and existing annotations
    """
    print("\n=== Annotating Cell Types ===")
    
    # Use pre-existing cell annotations if available
    if 'cell_type' in adata.obs.columns:
        print("Using pre-existing cell type annotations")
        cell_types = adata.obs['cell_type'].value_counts()
        print(f"Found {len(cell_types)} cell types")
        for cell_type, count in cell_types.items():
            print(f"  - {cell_type}: {count} cells")
    else:
        # Define known cell type markers
        # These are common markers for lung cancer/tissue cell types
        cell_markers = {
            'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19'],
            'Endothelial': ['PECAM1', 'VWF', 'CDH5'],
            'Fibroblasts': ['COL1A1', 'DCN', 'LUM', 'COL3A1'],
            'T-cells': ['CD3D', 'CD3E', 'CD3G', 'CD8A', 'CD4'],
            'B-cells': ['CD79A', 'CD79B', 'MS4A1', 'CD19'],
            'Myeloid': ['LYZ', 'CD68', 'CSF1R', 'MARCO'],
            'NK cells': ['NCAM1', 'NKG7', 'KLRD1'],
            'Alveolar cells': ['SFTPC', 'SFTPB', 'SFTPD'],
        }
        
        # Calculate gene expression scores for each cell type
        print("Calculating cell type scores based on marker genes...")
        for cell_type, markers in cell_markers.items():
            # Filter for genes that exist in our dataset
            available_markers = [m for m in markers if m in adata.var_names]
            if available_markers:
                sc.tl.score_genes(adata, available_markers, score_name=f"{cell_type}_score")
        
        # Assign cell types based on highest score
        print("Assigning cell types based on marker scores...")
        score_columns = [col for col in adata.obs.columns if col.endswith('_score')]
        
        if score_columns:
            adata.obs['predicted_cell_type'] = adata.obs[score_columns].idxmax(axis=1)
            adata.obs['predicted_cell_type'] = adata.obs['predicted_cell_type'].str.replace('_score', '')
            
            # Ensure the predicted_cell_type column is categorical
            adata.obs['predicted_cell_type'] = adata.obs['predicted_cell_type'].astype('category')
            
            # Plot cell type distribution
            plt.figure(figsize=(10, 6))
            adata.obs['predicted_cell_type'].value_counts().plot(kind='bar')
            plt.title('Predicted Cell Type Distribution')
            plt.xlabel('Cell Type')
            plt.ylabel('Number of Cells')
            plt.tight_layout()
            plt.savefig(f"{VIS_DIR}/predicted_cell_type_distribution.pdf")
            plt.close()
    
    # Plot UMAP with cell types
    cell_type_column = 'cell_type' if 'cell_type' in adata.obs.columns else 'predicted_cell_type'
    if cell_type_column in adata.obs.columns:
        # Plot UMAP with cell types - with improved label spacing
        plt.figure(figsize=(12, 10))  # Larger figure
        # Format category names for display by replacing underscores with spaces
        if cell_type_column in adata.obs:
            # Ensure the column is categorical before accessing .cat
            if not pd.api.types.is_categorical_dtype(adata.obs[cell_type_column]):
                adata.obs[cell_type_column] = adata.obs[cell_type_column].astype('category')
                
            # Get current categories and create a mapping to nicely formatted versions
            categories = adata.obs[cell_type_column].cat.categories
            category_map = {cat: cat.replace('_', ' ').replace('-', ' ').title() for cat in categories}
            
            # Create a temporary column with formatted names for plotting
            temp_col = f"{cell_type_column}_formatted"
            adata.obs[temp_col] = adata.obs[cell_type_column].map(category_map)
            
            # Make sure the new column is categorical with the right categories
            adata.obs[temp_col] = adata.obs[temp_col].astype('category')
            adata.obs[temp_col] = adata.obs[temp_col].cat.reorder_categories([category_map[cat] for cat in categories])
            
            # Plot with the formatted names
            sc.pl.umap(adata, color=temp_col, 
                     legend_loc='right',  # Move legend to the right instead of on data 
                     legend_fontsize=8,   # Smaller legend font
                     title='Cell Types', 
                     s=30,  # Smaller point size
                     show=False)
            
            # Delete temporary column after plotting
            del adata.obs[temp_col]
        else:
            # Fallback: use the original column
            sc.pl.umap(adata, color=cell_type_column, 
                     legend_loc='right',
                     legend_fontsize=8,
                     title='Cell Types', 
                     s=30,
                     show=False)
        
        # Get current axes and adjust
        ax = plt.gca()
        ax.set_title('Cell Types UMAP', fontsize=14, pad=20)  # Add more padding to title
        
        cell_types_umap_file = f"{VIS_DIR}/cell_types_umap.pdf"
        plt.savefig(cell_types_umap_file, bbox_inches='tight', dpi=150)
        plt.close()
        print(f"Saved cell types UMAP to {cell_types_umap_file}")
    
    return adata

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
    
    # 1. Male vs Female
    if 'Sex' in adata.obs.columns:
        print("Comparing cell types: Male vs. Female")
        male_vs_female = pd.crosstab(
            adata.obs[cell_type_column], 
            adata.obs['Sex'], 
            normalize='columns'
        ) * 100  # Convert to percentage
        
        plt.figure(figsize=(14, 9))
        male_vs_female.plot(kind='bar', stacked=False)
        plt.title(f'{cell_type_display} Distribution: Male vs. Female', fontsize=14, pad=20)
        plt.ylabel('Percentage of Cells', fontsize=12, labelpad=10)
        plt.xticks(rotation=45, ha='right', fontsize=10)
        plt.yticks(fontsize=10)
        plt.legend(fontsize=10, title_fontsize=11)
        plt.subplots_adjust(bottom=0.2)  # More space for x-labels
        plt.tight_layout()
        plt.savefig(f"{VIS_DIR}/male_vs_female_cell_types.pdf", dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved male vs. female comparison to {VIS_DIR}/male_vs_female_cell_types.pdf")
    
    # 2. Smoking status
    if 'Smoking' in adata.obs.columns:
        print("Comparing cell types by smoking status")
        smoking_groups = pd.crosstab(
            adata.obs[cell_type_column], 
            adata.obs['Smoking'], 
            normalize='columns'
        ) * 100
        
        plt.figure(figsize=(14, 9))
        smoking_groups.plot(kind='bar', stacked=False)
        plt.title(f'{cell_type_display} Distribution by Smoking Status', fontsize=14, pad=20)
        plt.ylabel('Percentage of Cells', fontsize=12, labelpad=10)
        plt.xticks(rotation=45, ha='right', fontsize=10)
        plt.yticks(fontsize=10)
        plt.legend(fontsize=10, title_fontsize=11)
        plt.subplots_adjust(bottom=0.2)  # More space for x-labels
        plt.tight_layout()
        plt.savefig(f"{VIS_DIR}/smoking_status_cell_types.pdf", dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved smoking status comparison to {VIS_DIR}/smoking_status_cell_types.pdf")
    
    # 3. Tumor vs. Normal
    if 'tissue_type' in adata.obs.columns:
        print("Comparing cell types: Tumor vs. Normal tissue")
        tumor_vs_normal = pd.crosstab(
            adata.obs[cell_type_column], 
            adata.obs['tissue_type'], 
            normalize='columns'
        ) * 100
        
        plt.figure(figsize=(14, 9))
        tumor_vs_normal.plot(kind='bar', stacked=False)
        plt.title(f'{cell_type_display} Distribution: Tumor vs. Normal Tissue', fontsize=14, pad=20)
        plt.ylabel('Percentage of Cells', fontsize=12, labelpad=10)
        plt.xticks(rotation=45, ha='right', fontsize=10)
        plt.yticks(fontsize=10)
        plt.legend(fontsize=10, title_fontsize=11)
        plt.subplots_adjust(bottom=0.2)  # More space for x-labels
        plt.tight_layout()
        plt.savefig(f"{VIS_DIR}/tumor_vs_normal_cell_types.pdf", dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved tumor vs. normal comparison to {VIS_DIR}/tumor_vs_normal_cell_types.pdf")
    
    return adata

def generate_visualizations(adata):
    """
    Create comprehensive cell type distribution visualizations
    """
    print("\n=== Generating Visualizations ===")
    
    # Determine which cell type column to use
    cell_type_column = 'cell_type' if 'cell_type' in adata.obs.columns else 'predicted_cell_type'
    
    # Ensure cell_type_column is categorical
    if not pd.api.types.is_categorical_dtype(adata.obs[cell_type_column]):
        adata.obs[cell_type_column] = adata.obs[cell_type_column].astype('category')
        
    # Create a display-friendly version of the column name for plot titles
    cell_type_display = 'Cell Type' if cell_type_column == 'cell_type' else 'Predicted Cell Type'
    if cell_type_column not in adata.obs.columns:
        print(f"Error: No {cell_type_column} column found in the data")
        return
    
    # 1. Overall cell type distribution
    plt.figure(figsize=(12, 6))
    adata.obs[cell_type_column].value_counts().plot(kind='bar')
    plt.title('Overall Cell Type Distribution')
    plt.ylabel('Number of Cells')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{VIS_DIR}/overall_cell_type_distribution.pdf")
    plt.close()
    
    # 2. Cell type distribution by patient
    if 'Patient id' in adata.obs.columns:
        plt.figure(figsize=(14, 8))
        patient_cell_types = pd.crosstab(
            adata.obs['Patient id'], 
            adata.obs[cell_type_column], 
            normalize='index'
        ) * 100
        
        patient_cell_types.plot(kind='bar', stacked=True, colormap='tab20')
        plt.title('Cell Type Distribution by Patient')
        plt.ylabel('Percentage of Cells')
        plt.xlabel('Patient ID')
        plt.xticks(rotation=45, ha='right')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(f"{VIS_DIR}/patient_cell_type_distribution.pdf")
        plt.close()
    
    # 3. Heatmap of gene expression by cell type
    print("Generating gene expression heatmap by cell type...")
    n_groups = min(len(adata.obs[cell_type_column].unique()), 8)  # Limit to top 8 cell types
    sc.tl.dendrogram(adata, groupby=cell_type_column)
    # Select top marker genes for each cell type to plot in the heatmap
    # We'll use the matrix plot's native features to improve visualization
    print("Selecting top genes for heatmap visualization...")
    
    # Take top N genes from previously computed markers
    marker_genes_per_cluster = {}
    genes_to_show = []
    if 'rank_genes_groups' in adata.uns and adata.uns['rank_genes_groups']:
        try:
            # Create dictionary of marker genes grouped by cell type for better visualization
            for cluster in adata.obs[cell_type_column].cat.categories:
                # Get top 5 marker genes for each cell type
                markers = sc.get.rank_genes_groups_df(adata, group=cluster)
                if not markers.empty:
                    top_markers = markers.sort_values('scores', ascending=False).head(5)['names'].tolist()
                    marker_genes_per_cluster[cluster] = top_markers
                    genes_to_show.extend(top_markers)
        except Exception as e:
            print(f"Could not extract markers by cell type: {e}")
            # Fall back to using top genes by expression
            genes_to_show = adata.var_names[:20].tolist()
    else:
        # Fall back to using top genes by expression
        genes_to_show = adata.var_names[:20].tolist()
    
    # Remove duplicates while preserving order
    genes_to_show = list(dict.fromkeys(genes_to_show))
    
    # Limit to a reasonable number to prevent overcrowding
    genes_to_show = genes_to_show[:25]
    
    print(f"Creating expression heatmap for {len(genes_to_show)} genes across {len(adata.obs[cell_type_column].cat.categories)} cell types...")
    
    # Create a temporary column with formatted cell type names for plotting
    temp_col = f"{cell_type_column}_display"
    if cell_type_column in adata.obs:
        # Ensure the column is categorical before accessing .cat
        if not pd.api.types.is_categorical_dtype(adata.obs[cell_type_column]):
            adata.obs[cell_type_column] = adata.obs[cell_type_column].astype('category')
            
        # Format category names by replacing underscores with spaces and capitalizing properly
        categories = adata.obs[cell_type_column].cat.categories
        category_map = {cat: cat.replace('_', ' ').replace('-', ' ').title() for cat in categories}
        
        # Create formatted column
        adata.obs[temp_col] = adata.obs[cell_type_column].map(category_map)
        adata.obs[temp_col] = adata.obs[temp_col].astype('category')
        adata.obs[temp_col] = adata.obs[temp_col].cat.reorder_categories([category_map[cat] for cat in categories])
        plot_groupby = temp_col
    else:
        plot_groupby = cell_type_column
    
    # Format gene names if needed
    formatted_genes = [gene.replace('_', ' ') for gene in genes_to_show]
    
    # Use the MatrixPlot class directly for more customization options
    matrix_plot = sc.pl.matrixplot(
        adata, 
        genes_to_show, 
        groupby=plot_groupby,
        dendrogram=True, 
        standard_scale='var', 
        cmap='viridis',
        colorbar_title='Scaled\nexpression',
        figsize=(16, 12),
        return_fig=True  # Return the figure object for customization
    )
    
    # Clean up temporary column
    if temp_col in adata.obs.columns:
        del adata.obs[temp_col]
    
    # Apply styling improvements
    matrix_plot.style(edge_color='black', edge_lw=0.5)
    
    # Add better axis labels
    matrix_plot.add_totals().style(edge_color='black')
    plt.gcf().text(0.5, 0.02, 'Cell Types', ha='center', fontsize=14)
    plt.gcf().text(0.02, 0.5, 'Marker Genes', va='center', rotation='vertical', fontsize=14)
    
    # Get the figure and customize it further
    fig = matrix_plot.get_figure()
    
    # Adjust gene labels (make them smaller and rotated)
    ax_dict = matrix_plot.get_axes()
    if 'mainplot_ax' in ax_dict:
        main_ax = ax_dict['mainplot_ax']
        
        # Adjust x-tick labels (genes)
        if not matrix_plot._swap_axes:
            for tick in main_ax.get_xticklabels():
                tick.set_fontsize(9)
                tick.set_rotation(45)
                tick.set_ha('right')
        
        # Adjust y-tick labels (cell types)
        for tick in main_ax.get_yticklabels():
            tick.set_fontsize(10)
    
    # Save the customized figure
    heatmap_file = f"{VIS_DIR}/cell_type_expression_heatmap.pdf"
    fig.savefig(heatmap_file, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"Saved expression heatmap to {heatmap_file}")
    
    # Create a simplified version without dendrogram for easier interpretation
    print("Creating simplified expression heatmap...")
    
    # Create a temporary column with formatted cell type names for the simplified plot
    temp_col = f"{cell_type_column}_display"
    if cell_type_column in adata.obs:
        # Ensure the column is categorical before accessing .cat
        if not pd.api.types.is_categorical_dtype(adata.obs[cell_type_column]):
            adata.obs[cell_type_column] = adata.obs[cell_type_column].astype('category')
            
        # Format category names by replacing underscores with spaces and capitalizing properly
        categories = adata.obs[cell_type_column].cat.categories
        category_map = {cat: cat.replace('_', ' ').replace('-', ' ').title() for cat in categories}
        
        # Create formatted column
        adata.obs[temp_col] = adata.obs[cell_type_column].map(category_map)
        adata.obs[temp_col] = adata.obs[temp_col].astype('category')
        adata.obs[temp_col] = adata.obs[temp_col].cat.reorder_categories([category_map[cat] for cat in categories])
        simple_plot_groupby = temp_col
    else:
        simple_plot_groupby = cell_type_column
    
    sc.pl.matrixplot(
        adata, 
        genes_to_show[:15],  # Limit to fewer genes for better clarity
        groupby=simple_plot_groupby,
        standard_scale='var', 
        cmap='viridis',
        figsize=(12, 10),
        colorbar_title='Scaled\nexpression',
        title=f"Gene Expression by {cell_type_display}",
        show=False
    )
    
    # Clean up temporary column
    if temp_col in adata.obs.columns:
        del adata.obs[temp_col]
    plt.tight_layout(pad=3.0)
    simple_heatmap_file = f"{VIS_DIR}/simple_expression_heatmap.pdf"
    plt.savefig(simple_heatmap_file, bbox_inches='tight', dpi=150)
    plt.close()
    print(f"Saved simplified heatmap to {simple_heatmap_file}")
    
    print(f"All visualizations saved to {VIS_DIR}")

def main():
    """
    Main function to run the entire analysis pipeline
    """
    print("\n=== GSE131907 Lung Cancer Cell Type Annotation Pipeline ===\n")
    
    # Set pandas display options and warning handling
    pd.set_option('display.max_columns', 100)
    pd.set_option('display.max_rows', 100)
    
    # Create output directories if they don't exist
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(VIS_DIR, exist_ok=True)
    
    # Load the data
    adata = load_data()
    if adata is None:
        return
    
    try:
        # Identify marker genes for each cluster
        adata = identify_marker_genes(adata)
        
        # Make sure columns are properly typed before annotation
        for col in adata.obs.columns:
            if col.endswith('_type') and not pd.api.types.is_categorical_dtype(adata.obs[col]):
                print(f"Converting {col} to categorical type")
                adata.obs[col] = adata.obs[col].astype('category')
        
        # Annotate cell types
        adata = annotate_cell_types(adata)
        
        # Ensure cell type columns are categorical after annotation
        for col in ['cell_type', 'predicted_cell_type']:
            if col in adata.obs.columns and not pd.api.types.is_categorical_dtype(adata.obs[col]):
                print(f"Converting {col} to categorical type after annotation")
                adata.obs[col] = adata.obs[col].astype('category')
        
        # Compare demographics
        adata = compare_demographics(adata)
        
        # Generate visualizations
        generate_visualizations(adata)
        
        # Visualize marker genes
        visualize_marker_genes(adata)
        
        # Save the annotated data
        print("\nSaving annotated data...")
        adata.write(f"{RESULTS_DIR}/GSE131907_annotated.h5ad")
        
        print("\n=== Analysis Pipeline Completed Successfully! ===")
    except Exception as e:
        print(f"\nError during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
