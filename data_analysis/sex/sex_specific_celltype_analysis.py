#!/usr/bin/env python3
"""
Sex-Specific Cell Type DEG Analysis
-----------------------------------
This script analyzes sex-specific differential expression genes (DEGs) across different
cell types to identify markers responsible for higher cancer risk in males vs females.

It provides:
1. Global statistics on DEG counts and patterns by sex and cell type
2. Identification of key sex-specific markers with highest fold changes and significance
3. Comparative analysis of cell type specific patterns between sexes
4. Visualization of results through heatmaps and volcano plots
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import glob
from scipy import stats
import argparse

# Set up global styling for plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set(font_scale=1.2)
sns.set_style("white")


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Cell type-specific sex difference analysis in cancer')
    parser.add_argument('--base_dir', type=str, default='/home/mukul.ranjan/Documents/cancer/SCELA/sex_tumor_interaction',
                      help='Base directory containing sex-specific data')
    parser.add_argument('--output_dir', type=str, default='sex_celltype_results',
                      help='Directory to save analysis results')
    parser.add_argument('--cell_types', type=str, nargs='+', 
                      help='Specific cell types to analyze (default: all)')
    parser.add_argument('--top_genes', type=int, default=20,
                      help='Number of top genes to analyze in detail')
    parser.add_argument('--min_fc', type=float, default=1.5,
                      help='Minimum fold change threshold (log2)')
    parser.add_argument('--max_pval', type=float, default=0.01,
                      help='Maximum adjusted p-value threshold')
    return parser.parse_args()


def load_global_degs(base_dir):
    """Load and analyze the global DEGs for males and females."""
    male_degs = pd.read_csv(os.path.join(base_dir, 'Male/tumor_vs_normal_significant.csv'))
    female_degs = pd.read_csv(os.path.join(base_dir, 'Female/tumor_vs_normal_significant.csv'))
    
    male_specific = pd.read_csv(os.path.join(base_dir, 'male_specific_tumor_degs.csv'))
    female_specific = pd.read_csv(os.path.join(base_dir, 'female_specific_tumor_degs.csv'))
    shared = pd.read_csv(os.path.join(base_dir, 'shared_tumor_degs.csv'))
    
    # Create summary statistics
    summary = {
        'Male': {
            'total_degs': len(male_degs),
            'up_regulated': len(male_degs[male_degs['logfoldchanges'] > 0]),
            'down_regulated': len(male_degs[male_degs['logfoldchanges'] < 0]),
            'highly_significant': len(male_degs[male_degs['pvals'] < 0.001]),
            'male_specific': len(male_specific),
            'median_fold_change': male_degs['logfoldchanges'].median(),
            'max_fold_change': male_degs['logfoldchanges'].max(),
            'min_fold_change': male_degs['logfoldchanges'].min()
        },
        'Female': {
            'total_degs': len(female_degs),
            'up_regulated': len(female_degs[female_degs['logfoldchanges'] > 0]),
            'down_regulated': len(female_degs[female_degs['logfoldchanges'] < 0]),
            'highly_significant': len(female_degs[female_degs['pvals'] < 0.001]),
            'female_specific': len(female_specific),
            'median_fold_change': female_degs['logfoldchanges'].median(),
            'max_fold_change': female_degs['logfoldchanges'].max(),
            'min_fold_change': female_degs['logfoldchanges'].min()
        },
        'Shared': {
            'shared_degs': len(shared)
        }
    }
    
    return {
        'male_degs': male_degs,
        'female_degs': female_degs,
        'male_specific': male_specific,
        'female_specific': female_specific,
        'shared': shared,
        'summary': summary
    }


def get_cell_types(base_dir):
    """Get all available cell types from the directory structure."""
    cell_types_dir = os.path.join(base_dir, 'cell_types')
    cell_types = [d for d in os.listdir(cell_types_dir) 
                 if os.path.isdir(os.path.join(cell_types_dir, d)) and 
                 os.path.exists(os.path.join(cell_types_dir, d, 'Male')) and
                 os.path.exists(os.path.join(cell_types_dir, d, 'Female'))]
    return cell_types


def load_cell_specific_degs(base_dir, cell_type, logger=None):
    """Load and analyze DEGs for a specific cell type in males and females."""
    cell_type_path = os.path.join(base_dir, 'cell_types', cell_type)
    
    # Check if both files exist
    male_file = os.path.join(cell_type_path, 'Male', 'tumor_vs_normal_significant.csv')
    female_file = os.path.join(cell_type_path, 'Female', 'tumor_vs_normal_significant.csv')
    
    if not (os.path.exists(male_file) and os.path.exists(female_file)):
        if logger:
            logger.warning(f"Missing DEG files for {cell_type}")
        return None
    
    male_degs = pd.read_csv(male_file)
    female_degs = pd.read_csv(female_file)
    
    if logger:
        logger.info(f"\n=== CELL TYPE: {cell_type} ===")
        logger.info(f"Loaded {len(male_degs)} male DEGs and {len(female_degs)} female DEGs")
        
        # Log distribution of fold changes
        m_fc_quantiles = male_degs['logfoldchanges'].quantile([0, 0.25, 0.5, 0.75, 1]).to_dict()
        f_fc_quantiles = female_degs['logfoldchanges'].quantile([0, 0.25, 0.5, 0.75, 1]).to_dict()
        
        logger.info("Male log2FC distribution:")
        logger.info(f"  Min: {m_fc_quantiles[0]:.2f}, Q1: {m_fc_quantiles[0.25]:.2f}, Median: {m_fc_quantiles[0.5]:.2f}, Q3: {m_fc_quantiles[0.75]:.2f}, Max: {m_fc_quantiles[1]:.2f}")
        
        logger.info("Female log2FC distribution:")
        logger.info(f"  Min: {f_fc_quantiles[0]:.2f}, Q1: {f_fc_quantiles[0.25]:.2f}, Median: {f_fc_quantiles[0.5]:.2f}, Q3: {f_fc_quantiles[0.75]:.2f}, Max: {f_fc_quantiles[1]:.2f}")
    
    # Add gender column
    male_degs['sex'] = 'Male'
    female_degs['sex'] = 'Female'
    
    # Create male-specific, female-specific, and shared gene sets
    male_genes = set(male_degs['names'])
    female_genes = set(female_degs['names'])
    
    male_specific_genes = male_genes - female_genes
    female_specific_genes = female_genes - male_genes
    shared_genes = male_genes.intersection(female_genes)
    
    if logger:
        logger.info(f"Gene set analysis:")
        logger.info(f"  Male-specific genes: {len(male_specific_genes)}")
        logger.info(f"  Female-specific genes: {len(female_specific_genes)}")
        logger.info(f"  Shared genes: {len(shared_genes)}")
    
    # Filter DEGs to get the specific sets
    male_specific = male_degs[male_degs['names'].isin(male_specific_genes)]
    female_specific = female_degs[female_degs['names'].isin(female_specific_genes)]
    
    # Get the shared DEGs
    male_shared = male_degs[male_degs['names'].isin(shared_genes)]
    female_shared = female_degs[female_degs['names'].isin(shared_genes)]
    
    # Merge shared DEGs to calculate fold change differences
    shared_degs = male_shared.merge(female_shared, on='names', suffixes=('_male', '_female'))
    shared_degs['fc_diff'] = shared_degs['logfoldchanges_male'] - shared_degs['logfoldchanges_female']
    shared_degs['abs_fc_diff'] = abs(shared_degs['fc_diff'])
    shared_degs['direction_switch'] = (shared_degs['logfoldchanges_male'] * 
                                       shared_degs['logfoldchanges_female']) < 0
    
    # Log top genes with opposite regulation
    if logger and len(shared_degs) > 0:
        direction_switch_count = len(shared_degs[shared_degs['direction_switch']])
        large_diff_count = len(shared_degs[shared_degs['abs_fc_diff'] > 2])
        
        logger.info(f"Shared gene analysis:")
        logger.info(f"  Genes with opposite regulation: {direction_switch_count} ({direction_switch_count/len(shared_degs)*100:.1f}%)")
        logger.info(f"  Genes with large sex difference (abs diff > 2): {large_diff_count} ({large_diff_count/len(shared_degs)*100:.1f}%)")
        
        if direction_switch_count > 0:
            top_switches = shared_degs[shared_degs['direction_switch']].sort_values('abs_fc_diff', ascending=False).head(5)
            logger.info("  Top 5 genes with opposite regulation:")
            for idx, row in top_switches.iterrows():
                logger.info(f"    {row['names']}: Male log2FC = {row['logfoldchanges_male']:.2f}, Female log2FC = {row['logfoldchanges_female']:.2f}")
    
    # Create summary statistics
    summary = {
        'Male': {
            'total_degs': len(male_degs),
            'up_regulated': len(male_degs[male_degs['logfoldchanges'] > 0]),
            'down_regulated': len(male_degs[male_degs['logfoldchanges'] < 0]),
            'highly_significant': len(male_degs[male_degs['pvals'] < 0.001]),
            'male_specific': len(male_specific),
            'median_fold_change': male_degs['logfoldchanges'].median() if len(male_degs) > 0 else 0,
            'max_fold_change': male_degs['logfoldchanges'].max() if len(male_degs) > 0 else 0
        },
        'Female': {
            'total_degs': len(female_degs),
            'up_regulated': len(female_degs[female_degs['logfoldchanges'] > 0]),
            'down_regulated': len(female_degs[female_degs['logfoldchanges'] < 0]),
            'highly_significant': len(female_degs[female_degs['pvals'] < 0.001]),
            'female_specific': len(female_specific),
            'median_fold_change': female_degs['logfoldchanges'].median() if len(female_degs) > 0 else 0,
            'max_fold_change': female_degs['logfoldchanges'].max() if len(female_degs) > 0 else 0
        },
        'Shared': {
            'shared_degs': len(shared_degs),
            'direction_switches': len(shared_degs[shared_degs['direction_switch']]),
            'large_diff_genes': len(shared_degs[shared_degs['abs_fc_diff'] > 2])
        }
    }
    
    # Log regulation patterns
    if logger:
        logger.info(f"Regulation patterns:")
        logger.info(f"  Male: {summary['Male']['up_regulated']} up-regulated, {summary['Male']['down_regulated']} down-regulated (ratio: {summary['Male']['up_regulated']/max(1, summary['Male']['down_regulated']):.2f})")
        logger.info(f"  Female: {summary['Female']['up_regulated']} up-regulated, {summary['Female']['down_regulated']} down-regulated (ratio: {summary['Female']['up_regulated']/max(1, summary['Female']['down_regulated']):.2f})")
        
        # Log top genes for both sexes
        if len(male_degs) > 0:
            top_male_genes = male_degs.sort_values('logfoldchanges', ascending=False).head(5)
            logger.info("Top 5 up-regulated genes in males:")
            for idx, row in top_male_genes.iterrows():
                logger.info(f"  {row['names']}: log2FC = {row['logfoldchanges']:.2f}, p-val = {row['pvals']:.2e}")
        
        if len(female_degs) > 0:
            top_female_genes = female_degs.sort_values('logfoldchanges', ascending=False).head(5)
            logger.info("Top 5 up-regulated genes in females:")
            for idx, row in top_female_genes.iterrows():
                logger.info(f"  {row['names']}: log2FC = {row['logfoldchanges']:.2f}, p-val = {row['pvals']:.2e}")
    
    return {
        'male_degs': male_degs,
        'female_degs': female_degs,
        'male_specific': male_specific,
        'female_specific': female_specific,
        'shared_degs': shared_degs,
        'summary': summary
    }


def get_top_sex_markers(cell_data, top_n=20, min_fc=1.5, max_pval=0.01):
    """Identify top markers specific to each sex based on fold change and significance."""
    if cell_data is None:
        return None
    
    # Filter by fold change and p-value
    male_filtered = cell_data['male_specific'][
        (cell_data['male_specific']['logfoldchanges'].abs() >= min_fc) & 
        (cell_data['male_specific']['pvals_adj'] <= max_pval)
    ]
    
    female_filtered = cell_data['female_specific'][
        (cell_data['female_specific']['logfoldchanges'].abs() >= min_fc) & 
        (cell_data['female_specific']['pvals_adj'] <= max_pval)
    ]
    
    # Sort by fold change
    male_top_fc = male_filtered.sort_values('logfoldchanges', ascending=False).head(top_n)
    female_top_fc = female_filtered.sort_values('logfoldchanges', ascending=False).head(top_n)
    
    # Sort by significance
    male_top_sig = male_filtered.sort_values('pvals').head(top_n)
    female_top_sig = female_filtered.sort_values('pvals').head(top_n)
    
    # Get genes with opposite regulation in shared DEGs
    direction_switch_genes = cell_data['shared_degs'][cell_data['shared_degs']['direction_switch']]
    
    # Sort by magnitude of fold change difference
    direction_switch_genes = direction_switch_genes.sort_values('abs_fc_diff', ascending=False)
    
    return {
        'male_top_fc': male_top_fc,
        'female_top_fc': female_top_fc,
        'male_top_sig': male_top_sig,
        'female_top_sig': female_top_sig,
        'direction_switch': direction_switch_genes.head(top_n)
    }


def plot_sex_celltype_heatmap(cell_summaries, output_dir, logger=None):
    """Create a heatmap visualization of DEG patterns across cell types by sex."""
    if logger:
        logger.info("\n=== CELL TYPE HEATMAP ANALYSIS ===\n")
        logger.info("Generating heatmap visualizations of sex differences across cell types")
        logger.info("Analysis includes metrics: total DEGs, up-regulated, down-regulated, male-specific, and female-specific genes")
        logger.info("For total DEGs: log2(female/male) ratio is shown - positive values indicate female bias, negative values indicate male bias")
        logger.info("For other metrics: normalized difference between sexes is shown (female proportion - male proportion)")
    
    # Prepare data for heatmap
    cell_types = list(cell_summaries.keys())
    
    if logger:
        logger.info(f"Analyzing {len(cell_types)} cell types for sex differences")
    
    # Create dataframes for different metrics
    metrics = ['total_degs', 'up_regulated', 'down_regulated', 'male_specific', 'female_specific']
    heatmap_data = {}
    
    # Store detailed values for logging
    detailed_data = []
    
    for metric in metrics:
        data = []
        for cell_type in cell_types:
            if cell_summaries[cell_type] is not None:
                male_value = cell_summaries[cell_type]['Male'].get(metric, 0)
                female_value = cell_summaries[cell_type]['Female'].get(metric, 0)
                
                # Calculate ratio or difference based on metric
                if metric == 'total_degs':
                    # Use log2 ratio for visualization
                    ratio = np.log2(female_value / male_value) if male_value > 0 else np.nan
                    
                    # Store detailed value for logging
                    detailed_data.append({
                        'cell_type': cell_type,
                        'metric': metric,
                        'male_value': male_value,
                        'female_value': female_value,
                        'ratio_or_diff': ratio,
                        'description': f"log2(female/male) ratio"
                    })
                else:
                    # Use normalized difference
                    total_male = cell_summaries[cell_type]['Male']['total_degs']
                    total_female = cell_summaries[cell_type]['Female']['total_degs']
                    
                    male_norm = male_value / total_male if total_male > 0 else 0
                    female_norm = female_value / total_female if total_female > 0 else 0
                    
                    ratio = female_norm - male_norm
                    
                    # Store detailed value for logging
                    detailed_data.append({
                        'cell_type': cell_type,
                        'metric': metric,
                        'male_value': male_value,
                        'female_value': female_value,
                        'male_proportion': male_norm,
                        'female_proportion': female_norm,
                        'ratio_or_diff': ratio,
                        'description': f"female proportion - male proportion"
                    })
                
                data.append({
                    'cell_type': cell_type,
                    'metric': metric,
                    'value': ratio
                })
        
        if data:
            heatmap_data[metric] = pd.DataFrame(data)
            
    # Log the most extreme sex differences by cell type
    if logger and detailed_data:
        detail_df = pd.DataFrame(detailed_data)
        
        # Log the most female-biased cell types (by total DEGs)
        if 'total_degs' in detail_df['metric'].values:
            total_degs_data = detail_df[detail_df['metric'] == 'total_degs'].sort_values('ratio_or_diff', ascending=False)
            
            logger.info("\nTop female-biased cell types (by total DEG count):")
            for i, row in total_degs_data.head(5).iterrows():
                if not np.isnan(row['ratio_or_diff']):
                    logger.info(f"  {row['cell_type']}: {row['female_value']} female DEGs vs {row['male_value']} male DEGs (log2 ratio: {row['ratio_or_diff']:.2f})")
            
            logger.info("\nTop male-biased cell types (by total DEG count):")
            for i, row in total_degs_data.tail(5).sort_values('ratio_or_diff').iterrows():
                if not np.isnan(row['ratio_or_diff']):
                    logger.info(f"  {row['cell_type']}: {row['male_value']} male DEGs vs {row['female_value']} female DEGs (log2 ratio: {row['ratio_or_diff']:.2f})")
    
    # Create and save heatmap for each metric
    for metric, df in heatmap_data.items():
        if not df.empty:
            plt.figure(figsize=(10, len(cell_types)/2 + 2))
            
            pivot_data = df.pivot(index='cell_type', columns='metric', values='value')
            
            # Determine colormap based on metric
            if metric == 'total_degs':
                cmap = 'RdBu_r'  # Red-Blue diverging colormap
                title = 'Log2 Ratio of Female:Male DEGs by Cell Type'
                center = 0
            else:
                cmap = 'RdBu_r'
                title = f'Female-Male Normalized Difference in {metric.replace("_", " ").title()}'
                center = 0
            
            # Create the heatmap
            sns.heatmap(pivot_data, cmap=cmap, center=center, 
                      linewidths=0.5, cbar_kws={'label': 'Female:Male Ratio (log2)'})
            
            plt.title(title)
            plt.tight_layout()
            plot_path = os.path.join(output_dir, f'sex_celltype_heatmap_{metric}.pdf')
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            if logger:
                logger.info(f"\nHeatmap visualization for {metric.replace('_', ' ')} saved to: {plot_path}")
                logger.info(f"  Plot description: {title}")
                logger.info(f"  Color scale: {cmap} (diverging colormap centered at {center})")
                
                # Log extreme values for this metric
                if not pivot_data.empty:
                    max_val = pivot_data.max().max()
                    min_val = pivot_data.min().min()
                    max_cell = pivot_data.stack().idxmax()
                    min_cell = pivot_data.stack().idxmin()
                    
                    if not np.isnan(max_val):
                        logger.info(f"  Strongest female bias: {max_cell[0]} ({max_val:.2f})")
                    if not np.isnan(min_val):
                        logger.info(f"  Strongest male bias: {min_cell[0]} ({min_val:.2f})")


def plot_celltype_volcano(cell_data, cell_type, output_dir, logger=None):
    """Create volcano plots highlighting sex-specific markers."""
    if cell_data is None:
        return
    
    if logger:
        logger.info(f"\n--- GENERATING VOLCANO PLOT FOR {cell_type} ---")
        logger.info(f"Plot type: Volcano plot showing differential gene expression between tumor and normal samples")
        logger.info(f"X-axis: Log2 Fold Change (tumor vs normal)")
        logger.info(f"Y-axis: -Log10 P-value (statistical significance)")
        logger.info(f"Color coding: Gray = all genes, Blue = male-specific genes, Red = female-specific genes")
        logger.info(f"Thresholds: Horizontal line at p-value = 0.05, Vertical lines at log2FC = ±1.5")
    
    # Combine male and female data
    combined = pd.concat([cell_data['male_degs'], cell_data['female_degs']])
    
    # Log statistical summary of the data being plotted
    if logger:
        # Calculate statistics for the volcano plot data
        male_sig = len(cell_data['male_degs'][cell_data['male_degs']['pvals'] < 0.05])
        female_sig = len(cell_data['female_degs'][cell_data['female_degs']['pvals'] < 0.05])
        
        male_high_fc = len(cell_data['male_degs'][cell_data['male_degs']['logfoldchanges'].abs() > 1.5])
        female_high_fc = len(cell_data['female_degs'][cell_data['female_degs']['logfoldchanges'].abs() > 1.5])
        
        male_sig_high_fc = len(cell_data['male_degs'][(cell_data['male_degs']['pvals'] < 0.05) & 
                                                     (cell_data['male_degs']['logfoldchanges'].abs() > 1.5)])
        female_sig_high_fc = len(cell_data['female_degs'][(cell_data['female_degs']['pvals'] < 0.05) & 
                                                         (cell_data['female_degs']['logfoldchanges'].abs() > 1.5)])
        
        logger.info(f"Volcano plot statistics:")
        logger.info(f"  Total genes plotted: {len(combined)}")
        logger.info(f"  Male genes: {len(cell_data['male_degs'])} total, {male_sig} significant (p<0.05), {male_high_fc} with |log2FC|>1.5, {male_sig_high_fc} both")
        logger.info(f"  Female genes: {len(cell_data['female_degs'])} total, {female_sig} significant (p<0.05), {female_high_fc} with |log2FC|>1.5, {female_sig_high_fc} both")
    
    # Create the volcano plot
    plt.figure(figsize=(12, 8))
    
    # Basic scatter of all points (gray)
    plt.scatter(combined['logfoldchanges'], -np.log10(combined['pvals']), 
               alpha=0.3, s=10, color='gray')
    
    # Overlay male-specific (blue)
    if len(cell_data['male_specific']) > 0:
        plt.scatter(cell_data['male_specific']['logfoldchanges'], 
                   -np.log10(cell_data['male_specific']['pvals']), 
                   alpha=0.7, s=20, color='blue', label='Male-Specific')
    
    # Overlay female-specific (red)
    if len(cell_data['female_specific']) > 0:
        plt.scatter(cell_data['female_specific']['logfoldchanges'], 
                   -np.log10(cell_data['female_specific']['pvals']), 
                   alpha=0.7, s=20, color='red', label='Female-Specific')
    
    # Label top genes
    top_genes = []
    
    if 'male_top_fc' in cell_data and len(cell_data['male_top_fc']) > 0:
        top_male = cell_data['male_top_fc'].nlargest(5, 'logfoldchanges')
        top_genes.extend([(row['names'], row['logfoldchanges'], row['pvals']) 
                         for _, row in top_male.iterrows()])
    
    if 'female_top_fc' in cell_data and len(cell_data['female_top_fc']) > 0:
        top_female = cell_data['female_top_fc'].nlargest(5, 'logfoldchanges')
        top_genes.extend([(row['names'], row['logfoldchanges'], row['pvals']) 
                         for _, row in top_female.iterrows()])
    
    # Log top gene information
    if logger and top_genes:
        logger.info(f"Top labeled genes in volcano plot:")
        for gene, fc, pval in top_genes:
            logger.info(f"  {gene}: log2FC = {fc:.2f}, p-value = {pval:.2e}, -log10(p) = {-np.log10(pval):.2f}")
    
    for gene, fc, pval in top_genes:
        plt.annotate(gene, (fc, -np.log10(pval)), 
                    xytext=(5, 5), textcoords='offset points', 
                    fontsize=8, fontweight='bold')
    
    # Add lines for thresholds
    plt.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.3)
    plt.axvline(-1.5, color='gray', linestyle='--', alpha=0.3)
    plt.axvline(1.5, color='gray', linestyle='--', alpha=0.3)
    
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 P-value')
    plt.title(f'Sex-Specific DEGs in {cell_type.replace("_", " ")}')
    plt.legend()
    plt.grid(False)
    
    # Save the plot
    plot_path = os.path.join(output_dir, f'{cell_type}_sex_volcano.pdf')
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    if logger:
        logger.info(f"Volcano plot saved to: {plot_path}")


def export_top_markers(cell_data_dict, output_dir, top_n=20, logger=None):
    """Export tables of top markers for each cell type."""
    if logger:
        logger.info("\n=== EXPORTING TOP SEX-SPECIFIC MARKERS ===\n")
        logger.info(f"Analyzing and exporting top {top_n} markers for each cell type, grouped by sex")
    
    all_cell_markers = {
        'male': [],
        'female': []
    }
    
    # Track markers for summary statistics
    cell_type_stats = {
        'male': {},
        'female': {}
    }
    
    for cell_type, cell_data in cell_data_dict.items():
        if cell_data is None or 'male_top_fc' not in cell_data:
            if logger:
                logger.warning(f"Skipping {cell_type}: missing required data")
            continue
        
        # Extract top markers by fold change
        for sex, df_key in [('male', 'male_top_fc'), ('female', 'female_top_fc')]:
            if len(cell_data[df_key]) > 0:
                # Add cell type and get top genes
                top_genes = cell_data[df_key].head(top_n).copy()
                top_genes['cell_type'] = cell_type
                
                # Record statistics for this cell type
                cell_type_stats[sex][cell_type] = len(top_genes)
                
                # Add to the overall list
                all_cell_markers[sex].append(top_genes)
                
                if logger:
                    # Log the top 3 genes for this cell type
                    if len(top_genes) > 0:
                        logger.info(f"Top {sex} markers for {cell_type}:")
                        for idx, row in top_genes.head(3).iterrows():
                            logger.info(f"  {row['names']}: log2FC = {row['logfoldchanges']:.2f}, p-value = {row['pvals']:.2e}")
    
    # Combine and save
    export_files = {}
    for sex in ['male', 'female']:
        if all_cell_markers[sex]:
            combined = pd.concat(all_cell_markers[sex])
            combined.sort_values(['cell_type', 'logfoldchanges'], ascending=[True, False], inplace=True)
            
            output_file = os.path.join(output_dir, f'{sex}_specific_celltype_markers.csv')
            combined.to_csv(output_file, index=False)
            export_files[sex] = output_file
            
            if logger:
                # Log summary statistics about exported markers
                total_markers = len(combined)
                total_cell_types = len(cell_type_stats[sex])
                avg_markers_per_celltype = total_markers / total_cell_types if total_cell_types > 0 else 0
                
                logger.info(f"\nExported {total_markers} {sex}-specific markers across {total_cell_types} cell types to {output_file}")
                logger.info(f"Average of {avg_markers_per_celltype:.1f} markers per cell type")
                
                # Cell types with most markers
                if cell_type_stats[sex]:
                    top_cells = pd.Series(cell_type_stats[sex]).sort_values(ascending=False)
                    logger.info(f"\nCell types with most {sex}-specific markers:")
                    for cell_type, count in top_cells.head(3).items():
                        logger.info(f"  {cell_type}: {count} markers")
                
                # Top markers by fold change
                logger.info(f"\nTop {sex}-specific markers across all cell types (by fold change):")
                for idx, row in combined.sort_values('logfoldchanges', ascending=False).head(5).iterrows():
                    logger.info(f"  {row['names']} ({row['cell_type']}): log2FC = {row['logfoldchanges']:.2f}, p-value = {row['pvals']:.2e}")
    
    # Create a table of top genes with direction switches
    direction_switch_genes = []
    switch_counts = {}
    
    if logger:
        logger.info("\n--- ANALYZING GENES WITH SEX-DEPENDENT DIRECTION SWITCHES ---")
        logger.info("Identifying genes that change expression direction between sexes")
        logger.info("These are genes that are up-regulated in one sex but down-regulated in the other")
    
    for cell_type, cell_data in cell_data_dict.items():
        if cell_data is None or 'direction_switch' not in cell_data:
            continue
            
        if len(cell_data['direction_switch']) > 0:
            switch_genes = cell_data['direction_switch'].head(top_n).copy()
            switch_genes['cell_type'] = cell_type
            direction_switch_genes.append(switch_genes)
            switch_counts[cell_type] = len(switch_genes)
            
            if logger and len(switch_genes) > 0:
                logger.info(f"\nFound {len(switch_genes)} direction-switching genes in {cell_type}")
                for idx, row in switch_genes.head(3).iterrows():
                    # Direction switch genes use the index as the gene name
                    # Example from logs: "JCHAIN: Male log2FC = -0.94, Female log2FC = 1.30"
                    gene_name = row.name if isinstance(row.name, str) else idx
                    
                    # Log fold changes using the actual column names from the dataframe
                    male_fc = row.get('male_logfoldchanges', row.get('logfoldchanges_male', 0))
                    female_fc = row.get('female_logfoldchanges', row.get('logfoldchanges_female', 0))
                    
                    if pd.notna(male_fc) and pd.notna(female_fc):
                        logger.info(f"  {gene_name}: Male log2FC = {male_fc:.2f}, Female log2FC = {female_fc:.2f}, difference = {abs(male_fc - female_fc):.2f}")
                    elif 'abs_fc_diff' in row:
                        logger.info(f"  {gene_name}: fold change difference = {row['abs_fc_diff']:.2f}")
                    else:
                        logger.info(f"  {gene_name}: direction switching gene")
    
    if direction_switch_genes:
        combined_switches = pd.concat(direction_switch_genes)
        combined_switches.sort_values(['abs_fc_diff'], ascending=False, inplace=True)
        
        output_file = os.path.join(output_dir, 'direction_switch_genes.csv')
        combined_switches.to_csv(output_file, index=False)
        
        if logger:
            # Log summary statistics about direction-switching genes
            total_switches = len(combined_switches)
            total_cell_types = len(switch_counts)
            
            logger.info(f"\nExported {total_switches} direction-switching genes across {total_cell_types} cell types to {output_file}")
            
            # Cell types with most direction switches
            if switch_counts:
                top_cells = pd.Series(switch_counts).sort_values(ascending=False)
                logger.info(f"\nCell types with most direction-switching genes:")
                for cell_type, count in top_cells.head(3).items():
                    logger.info(f"  {cell_type}: {count} genes")
            
            # Top direction-switching genes
            logger.info(f"\nTop direction-switching genes across all cell types (by absolute FC difference):")
            for idx, row in combined_switches.head(5).iterrows():
                # Direction switch genes use the index as the gene name
                gene_name = row.name if isinstance(row.name, str) else idx
                
                # Log the cell type and fold changes
                male_fc = row.get('male_logfoldchanges', row.get('logfoldchanges_male', 0))
                female_fc = row.get('female_logfoldchanges', row.get('logfoldchanges_female', 0))
                
                if pd.notna(male_fc) and pd.notna(female_fc):
                    logger.info(f"  {gene_name} ({row['cell_type']}): Male log2FC = {male_fc:.2f}, Female log2FC = {female_fc:.2f}, diff = {abs(male_fc - female_fc):.2f}")
                elif 'abs_fc_diff' in row:
                    logger.info(f"  {gene_name} ({row['cell_type']}): fold change difference = {row['abs_fc_diff']:.2f}")
                else:
                    logger.info(f"  {gene_name} ({row['cell_type']}): direction switching gene")


def analyze_pathway_enrichment(gene_lists, output_dir):
    """
    Simple enrichment analysis of gene lists against known pathways.
    In a real implementation, this would use a proper pathway analysis package.
    """
    try:
        import gseapy as gp
        
        for name, genes in gene_lists.items():
            if len(genes) < 10:  # Skip if too few genes
                continue
                
            # This is a placeholder - in real implementation would do proper enrichment
            # Example using Enrichr via gseapy
            try:
                enr = gp.enrichr(gene_list=genes, 
                                gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2021'],
                                outdir=None,
                                cutoff=0.05)
                
                # Save results
                if hasattr(enr, 'results'):
                    enr.results.to_csv(os.path.join(output_dir, f'{name}_pathway_enrichment.csv'))
            except Exception as e:
                print(f"Enrichment analysis failed for {name}: {str(e)}")
                
    except ImportError:
        print("gseapy not installed. Skipping pathway enrichment analysis.")
        print("Install with: pip install gseapy")


def setup_logging(output_dir):
    """Set up logging to save analysis results for LLM report generation."""
    log_file = os.path.join(output_dir, 'sex_analysis_report.log')
    
    # Configure logging
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger()

def main():
    """Main function to run the analysis."""
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Set up logging
    logger = setup_logging(args.output_dir)
    
    # Load global DEGs
    logger.info("Loading global DEG data...")
    global_deg_data = load_global_degs(args.base_dir)
    
    # Log global statistics for LLM report
    logger.info(f"=== GLOBAL DEG STATISTICS ===")
    logger.info(f"Male total DEGs: {global_deg_data['summary']['Male']['total_degs']}")
    logger.info(f"Female total DEGs: {global_deg_data['summary']['Female']['total_degs']}")
    logger.info(f"Male up-regulated: {global_deg_data['summary']['Male']['up_regulated']}")
    logger.info(f"Male down-regulated: {global_deg_data['summary']['Male']['down_regulated']}")
    logger.info(f"Female up-regulated: {global_deg_data['summary']['Female']['up_regulated']}")
    logger.info(f"Female down-regulated: {global_deg_data['summary']['Female']['down_regulated']}")
    logger.info(f"Male-specific genes: {global_deg_data['summary']['Male']['male_specific']}")
    logger.info(f"Female-specific genes: {global_deg_data['summary']['Female']['female_specific']}")
    logger.info(f"Shared genes: {global_deg_data['summary']['Shared']['shared_degs']}")
    
    # Get cell types to analyze
    if args.cell_types:
        cell_types = args.cell_types
    else:
        cell_types = get_cell_types(args.base_dir)
    
    logger.info(f"Analyzing {len(cell_types)} cell types: {', '.join(cell_types)}")
    
    # Initialize dictionaries to store results
    cell_data_dict = {}
    cell_summary_dict = {}
    
    # Analyze each cell type
    for cell_type in cell_types:
        logger.info(f"Processing {cell_type}...")
        
        # Load cell-specific DEGs
        cell_data = load_cell_specific_degs(args.base_dir, cell_type, logger)
        
        if cell_data is not None:
            cell_data_dict[cell_type] = cell_data
            cell_summary_dict[cell_type] = cell_data['summary']
            
            # Get top markers
            top_markers = get_top_sex_markers(cell_data, 
                                             top_n=args.top_genes, 
                                             min_fc=args.min_fc, 
                                             max_pval=args.max_pval)
            
            if top_markers:
                cell_data_dict[cell_type].update(top_markers)
            
            # Create volcano plot
            plot_celltype_volcano(cell_data_dict[cell_type], cell_type, args.output_dir, logger)
    
    # Create overall heatmap visualization
    logger.info("Creating cell type heatmap...")
    plot_sex_celltype_heatmap(cell_summary_dict, args.output_dir, logger)
    
    # Export tables of top markers
    logger.info("Exporting top markers...")
    export_top_markers(cell_data_dict, args.output_dir, top_n=args.top_genes, logger=logger)
    
    # Prepare gene lists for pathway analysis
    gene_lists = {}
    
    # Add global DEGs
    gene_lists['male_specific_global'] = global_deg_data['male_specific']['gene'].tolist()
    gene_lists['female_specific_global'] = global_deg_data['female_specific']['gene'].tolist()
    
    # Add cell-specific top genes
    for cell_type, cell_data in cell_data_dict.items():
        if cell_data is None or 'male_top_fc' not in cell_data or 'female_top_fc' not in cell_data:
            continue
            
        if len(cell_data['male_top_fc']) > 0:
            gene_lists[f'male_specific_{cell_type}'] = cell_data['male_top_fc']['names'].tolist()
            
        if len(cell_data['female_top_fc']) > 0:
            gene_lists[f'female_specific_{cell_type}'] = cell_data['female_top_fc']['names'].tolist()
    
    # Run pathway enrichment analysis
    logger.info("Running pathway enrichment analysis...")
    analyze_pathway_enrichment(gene_lists, args.output_dir)
    
    # Create a summary table of all cell types
    summary_rows = []
    for cell_type, summary in cell_summary_dict.items():
        if summary is not None:
            male_stats = summary['Male']
            female_stats = summary['Female']
            shared_stats = summary['Shared']
            
            row = {
                'cell_type': cell_type,
                'male_total_degs': male_stats['total_degs'],
                'female_total_degs': female_stats['total_degs'],
                'male_up': male_stats['up_regulated'],
                'male_down': male_stats['down_regulated'],
                'female_up': female_stats['up_regulated'],
                'female_down': female_stats['down_regulated'],
                'male_specific': male_stats['male_specific'],
                'female_specific': female_stats['female_specific'],
                'shared_degs': shared_stats.get('shared_degs', 0),
                'direction_switches': shared_stats.get('direction_switches', 0),
                'large_diff_genes': shared_stats.get('large_diff_genes', 0),
                'log2_ratio_f_to_m': np.log2(female_stats['total_degs'] / male_stats['total_degs']) 
                                    if male_stats['total_degs'] > 0 else np.nan
            }
            summary_rows.append(row)
    
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_df.sort_values('log2_ratio_f_to_m', ascending=False, inplace=True)
        summary_df.to_csv(os.path.join(args.output_dir, 'cell_type_summary.csv'), index=False)
    
    # Log cell type-specific findings for LLM report
    logger.info("\n=== CELL TYPE-SPECIFIC FINDINGS ===\n")
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        
        # Log top male-biased cell types
        male_biased = summary_df.sort_values('log2_ratio_f_to_m').head(5)
        logger.info("Top 5 male-biased cell types:")
        for idx, row in male_biased.iterrows():
            logger.info(f"  {row['cell_type']}: {row['male_total_degs']} male DEGs vs {row['female_total_degs']} female DEGs (log2 ratio: {row['log2_ratio_f_to_m']})")
        
        # Log top female-biased cell types
        female_biased = summary_df.sort_values('log2_ratio_f_to_m', ascending=False).head(5)
        logger.info("\nTop 5 female-biased cell types:")
        for idx, row in female_biased.iterrows():
            logger.info(f"  {row['cell_type']}: {row['female_total_degs']} female DEGs vs {row['male_total_degs']} male DEGs (log2 ratio: {row['log2_ratio_f_to_m']})")
    
    # Log marker gene highlights
    logger.info("\n=== TOP MARKER GENES ===\n")
    for sex in ['male', 'female']:
        if os.path.exists(os.path.join(args.output_dir, f'{sex}_specific_celltype_markers.csv')):
            markers_df = pd.read_csv(os.path.join(args.output_dir, f'{sex}_specific_celltype_markers.csv'))
            top_markers = markers_df.sort_values('logfoldchanges', ascending=False).head(20)
            
            logger.info(f"Top 20 {sex.capitalize()}-specific marker genes across all cell types:")
            for idx, row in top_markers.iterrows():
                logger.info(f"  {row['names']} in {row['cell_type']}: log2FC = {row['logfoldchanges']:.2f}, p-value = {row['pvals']:.2e}")
            logger.info("")
    
    # Log direction switch genes
    if os.path.exists(os.path.join(args.output_dir, 'direction_switch_genes.csv')):
        switch_df = pd.read_csv(os.path.join(args.output_dir, 'direction_switch_genes.csv'))
        top_switches = switch_df.sort_values('abs_fc_diff', ascending=False).head(10)
        
        logger.info("\nTop 10 genes with opposite regulation patterns between sexes:")
        for idx, row in top_switches.iterrows():
            logger.info(f"  {row['names']} in {row['cell_type']}: Male log2FC = {row['logfoldchanges_male']:.2f}, Female log2FC = {row['logfoldchanges_female']:.2f}")
    
    logger.info(f"\nAnalysis complete! Results saved to {args.output_dir}")


if __name__ == "__main__":
    main()
