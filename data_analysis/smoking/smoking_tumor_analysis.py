#!/usr/bin/env python3
"""
Smoking-Tumor Interaction Analysis
This script analyzes differential gene expression patterns between smoking categories
(Never smoker, Ex-smoker, Current smoker) across different cell types in lung cancer data.
It identifies genes with significant differential expression and highlights those with
opposite regulation patterns between smoking groups.
"""

import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import logging
import glob
import re
from scipy import stats
import math

# Configure logging
log_dir = "/home/mukul.ranjan/Documents/cancer/SCELA/smoking_tumor_analysis_results"
os.makedirs(log_dir, exist_ok=True)
log_file = os.path.join(log_dir, "smoking_analysis_report.log")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger()

# Define paths
base_dir = "/home/mukul.ranjan/Documents/cancer/SCELA"
smoking_dir = os.path.join(base_dir, "smoking_tumor_interaction")
output_dir = os.path.join(base_dir, "smoking_tumor_analysis_results")
os.makedirs(output_dir, exist_ok=True)

# Smoking categories
smoking_categories = ["Never", "Ex", "Cur"]

def load_smoking_comparison_data():
    """Load and return the smoking group comparison data."""
    comparison_file = os.path.join(smoking_dir, "smoking_group_deg_comparison.csv")
    return pd.read_csv(comparison_file, index_col=0)

def load_smoking_dependent_genes():
    """Load and return the smoking-dependent tumor genes."""
    dependent_file = os.path.join(smoking_dir, "smoking_dependent_tumor_genes.csv")
    return pd.read_csv(dependent_file)

def get_cell_type_list():
    """Get list of all cell types with data."""
    cell_types_dir = os.path.join(smoking_dir, "cell_types")
    return [d for d in os.listdir(cell_types_dir) 
            if os.path.isdir(os.path.join(cell_types_dir, d)) and 
            not d.startswith(".")]

def load_celltype_data(cell_type, smoking_status):
    """Load gene expression data for a specific cell type and smoking status."""
    cell_dir = os.path.join(smoking_dir, "cell_types", cell_type, smoking_status)
    
    # Check if directory exists
    if not os.path.exists(cell_dir):
        return None
    
    # Load significant DEGs
    sig_file = os.path.join(cell_dir, "tumor_vs_normal_significant.csv")
    if os.path.exists(sig_file):
        return pd.read_csv(sig_file)
    return None

def analyze_cell_type_deg_counts():
    """Analyze the counts of differentially expressed genes for each cell type and smoking status."""
    cell_types = get_cell_type_list()
    results = []
    
    logger.info("\n=== ANALYZING CELL TYPE DEG COUNTS ===\n")
    
    for cell_type in cell_types:
        logger.info(f"Processing {cell_type}...")
        
        for smoking in smoking_categories:
            data = load_celltype_data(cell_type, smoking)
            
            if data is not None:
                # Count up and down-regulated genes
                up_genes = sum(data["logfoldchanges"] > 0)
                down_genes = sum(data["logfoldchanges"] < 0)
                total_genes = len(data)
                
                results.append({
                    "Cell_Type": cell_type,
                    "Smoking_Status": smoking,
                    "Up_Regulated": up_genes,
                    "Down_Regulated": down_genes,
                    "Total_DEGs": total_genes
                })
                
                logger.info(f"  {smoking}: {total_genes} DEGs ({up_genes} up, {down_genes} down)")
            else:
                logger.info(f"  No data for {smoking}")
    
    # Create DataFrame and save to CSV
    df = pd.DataFrame(results)
    output_file = os.path.join(output_dir, "smoking_celltype_deg_counts.csv")
    df.to_csv(output_file, index=False)
    logger.info(f"Saved cell type DEG counts to {output_file}")
    
    return df

def identify_opposite_regulation_genes():
    """
    Identify genes with opposite regulation patterns between smoking categories.
    This function analyzes the smoking group comparison data to find genes that:
    1. Are upregulated in one smoking category but downregulated in another
    2. Show significant expression differences between smoking categories
    """
    logger.info("\n=== IDENTIFYING GENES WITH OPPOSITE REGULATION PATTERNS ===\n")
    
    # Load the smoking comparison data
    comp_data = load_smoking_comparison_data()
    
    # Convert padj columns to numeric, handling any non-numeric values
    for smoking in smoking_categories:
        padj_col = f"{smoking}_padj"
        comp_data[padj_col] = pd.to_numeric(comp_data[padj_col], errors='coerce')
    
    # Define significance thresholds
    log2fc_threshold = 0.5
    padj_threshold = 0.05
    
    # Lists to store different comparison results
    never_vs_ex = []
    never_vs_cur = []
    ex_vs_cur = []
    all_pattern_switches = []
    
    # Analyze Never vs Ex smokers
    logger.info("Analyzing Never vs Ex smokers...")
    for idx, row in comp_data.iterrows():
        gene = idx
        never_fc = row.get('Never_log2FC', 0)
        ex_fc = row.get('Ex_log2FC', 0)
        never_padj = row.get('Never_padj', 1)
        ex_padj = row.get('Ex_padj', 1)
        
        # Check if gene shows opposite regulation
        if (never_fc > log2fc_threshold and ex_fc < -log2fc_threshold and 
            never_padj < padj_threshold and ex_padj < padj_threshold):
            never_vs_ex.append({
                'Gene': gene,
                'Never_log2FC': never_fc,
                'Ex_log2FC': ex_fc,
                'Never_padj': never_padj,
                'Ex_padj': ex_padj,
                'Comparison': 'Never_vs_Ex',
                'Pattern': 'Up in Never, Down in Ex'
            })
            all_pattern_switches.append({
                'Gene': gene,
                'Never_log2FC': never_fc,
                'Ex_log2FC': ex_fc,
                'Cur_log2FC': row.get('Cur_log2FC', 0),
                'Never_padj': never_padj,
                'Ex_padj': ex_padj,
                'Cur_padj': row.get('Cur_padj', 1),
                'Pattern_Type': 'Never_vs_Ex',
                'Description': 'Up in Never, Down in Ex'
            })
        elif (never_fc < -log2fc_threshold and ex_fc > log2fc_threshold and 
              never_padj < padj_threshold and ex_padj < padj_threshold):
            never_vs_ex.append({
                'Gene': gene,
                'Never_log2FC': never_fc,
                'Ex_log2FC': ex_fc,
                'Never_padj': never_padj,
                'Ex_padj': ex_padj,
                'Comparison': 'Never_vs_Ex',
                'Pattern': 'Down in Never, Up in Ex'
            })
            all_pattern_switches.append({
                'Gene': gene,
                'Never_log2FC': never_fc,
                'Ex_log2FC': ex_fc,
                'Cur_log2FC': row.get('Cur_log2FC', 0),
                'Never_padj': never_padj,
                'Ex_padj': ex_padj,
                'Cur_padj': row.get('Cur_padj', 1),
                'Pattern_Type': 'Never_vs_Ex',
                'Description': 'Down in Never, Up in Ex'
            })
    
    # Analyze Never vs Current smokers
    logger.info("Analyzing Never vs Current smokers...")
    for idx, row in comp_data.iterrows():
        gene = idx
        never_fc = row.get('Never_log2FC', 0)
        cur_fc = row.get('Cur_log2FC', 0)
        never_padj = row.get('Never_padj', 1)
        cur_padj = row.get('Cur_padj', 1)
        
        # Check if gene shows opposite regulation
        if (never_fc > log2fc_threshold and cur_fc < -log2fc_threshold and 
            never_padj < padj_threshold and cur_padj < padj_threshold):
            never_vs_cur.append({
                'Gene': gene,
                'Never_log2FC': never_fc,
                'Cur_log2FC': cur_fc,
                'Never_padj': never_padj,
                'Cur_padj': cur_padj,
                'Comparison': 'Never_vs_Cur',
                'Pattern': 'Up in Never, Down in Current'
            })
            all_pattern_switches.append({
                'Gene': gene,
                'Never_log2FC': never_fc,
                'Ex_log2FC': row.get('Ex_log2FC', 0),
                'Cur_log2FC': cur_fc,
                'Never_padj': never_padj,
                'Ex_padj': row.get('Ex_padj', 1),
                'Cur_padj': cur_padj,
                'Pattern_Type': 'Never_vs_Cur',
                'Description': 'Up in Never, Down in Current'
            })
        elif (never_fc < -log2fc_threshold and cur_fc > log2fc_threshold and 
              never_padj < padj_threshold and cur_padj < padj_threshold):
            never_vs_cur.append({
                'Gene': gene,
                'Never_log2FC': never_fc,
                'Cur_log2FC': cur_fc,
                'Never_padj': never_padj,
                'Cur_padj': cur_padj,
                'Comparison': 'Never_vs_Cur',
                'Pattern': 'Down in Never, Up in Current'
            })
            all_pattern_switches.append({
                'Gene': gene,
                'Never_log2FC': never_fc,
                'Ex_log2FC': row.get('Ex_log2FC', 0),
                'Cur_log2FC': cur_fc,
                'Never_padj': never_padj,
                'Ex_padj': row.get('Ex_padj', 1),
                'Cur_padj': cur_padj,
                'Pattern_Type': 'Never_vs_Cur',
                'Description': 'Down in Never, Up in Current'
            })
    
    # Analyze Ex vs Current smokers
    logger.info("Analyzing Ex vs Current smokers...")
    for idx, row in comp_data.iterrows():
        gene = idx
        ex_fc = row.get('Ex_log2FC', 0)
        cur_fc = row.get('Cur_log2FC', 0)
        ex_padj = row.get('Ex_padj', 1)
        cur_padj = row.get('Cur_padj', 1)
        
        # Check if gene shows opposite regulation
        if (ex_fc > log2fc_threshold and cur_fc < -log2fc_threshold and 
            ex_padj < padj_threshold and cur_padj < padj_threshold):
            ex_vs_cur.append({
                'Gene': gene,
                'Ex_log2FC': ex_fc,
                'Cur_log2FC': cur_fc,
                'Ex_padj': ex_padj,
                'Cur_padj': cur_padj,
                'Comparison': 'Ex_vs_Cur',
                'Pattern': 'Up in Ex, Down in Current'
            })
            all_pattern_switches.append({
                'Gene': gene,
                'Never_log2FC': row.get('Never_log2FC', 0),
                'Ex_log2FC': ex_fc,
                'Cur_log2FC': cur_fc,
                'Never_padj': row.get('Never_padj', 1),
                'Ex_padj': ex_padj,
                'Cur_padj': cur_padj,
                'Pattern_Type': 'Ex_vs_Cur',
                'Description': 'Up in Ex, Down in Current'
            })
        elif (ex_fc < -log2fc_threshold and cur_fc > log2fc_threshold and 
              ex_padj < padj_threshold and cur_padj < padj_threshold):
            ex_vs_cur.append({
                'Gene': gene,
                'Ex_log2FC': ex_fc,
                'Cur_log2FC': cur_fc,
                'Ex_padj': ex_padj,
                'Cur_padj': cur_padj,
                'Comparison': 'Ex_vs_Cur',
                'Pattern': 'Down in Ex, Up in Current'
            })
            all_pattern_switches.append({
                'Gene': gene,
                'Never_log2FC': row.get('Never_log2FC', 0),
                'Ex_log2FC': ex_fc,
                'Cur_log2FC': cur_fc,
                'Never_padj': row.get('Never_padj', 1),
                'Ex_padj': ex_padj,
                'Cur_padj': cur_padj,
                'Pattern_Type': 'Ex_vs_Cur',
                'Description': 'Down in Ex, Up in Current'
            })
    
    # Log summary
    logger.info(f"Found {len(never_vs_ex)} genes with opposite regulation between Never and Ex smokers")
    logger.info(f"Found {len(never_vs_cur)} genes with opposite regulation between Never and Current smokers")
    logger.info(f"Found {len(ex_vs_cur)} genes with opposite regulation between Ex and Current smokers")
    logger.info(f"Total of {len(all_pattern_switches)} opposite regulation patterns identified")
    
    # Save results to CSV files
    if never_vs_ex:
        never_ex_df = pd.DataFrame(never_vs_ex)
        never_ex_file = os.path.join(output_dir, "never_vs_ex_opposite_regulation.csv")
        never_ex_df.to_csv(never_ex_file, index=False)
        logger.info(f"Saved Never vs Ex results to {never_ex_file}")
    
    if never_vs_cur:
        never_cur_df = pd.DataFrame(never_vs_cur)
        never_cur_file = os.path.join(output_dir, "never_vs_cur_opposite_regulation.csv")
        never_cur_df.to_csv(never_cur_file, index=False)
        logger.info(f"Saved Never vs Current results to {never_cur_file}")
    
    if ex_vs_cur:
        ex_cur_df = pd.DataFrame(ex_vs_cur)
        ex_cur_file = os.path.join(output_dir, "ex_vs_cur_opposite_regulation.csv")
        ex_cur_df.to_csv(ex_cur_file, index=False)
        logger.info(f"Saved Ex vs Current results to {ex_cur_file}")
    
    if all_pattern_switches:
        all_switches_df = pd.DataFrame(all_pattern_switches)
        all_switches_file = os.path.join(output_dir, "all_smoking_direction_switch_genes.csv")
        all_switches_df.to_csv(all_switches_file, index=False)
        logger.info(f"Saved all direction switch results to {all_switches_file}")
        
        # Sort by the magnitude of the direction change
        all_switches_df['magnitude'] = abs(all_switches_df['Never_log2FC'] - all_switches_df['Cur_log2FC'])
        top_switches = all_switches_df.sort_values('magnitude', ascending=False).head(20)
        
        # Log top genes with opposite regulation
        logger.info("\nTop genes with the largest opposite regulation effects between smoking categories:")
        for _, row in top_switches.iterrows():
            gene = row['Gene']
            pattern = row['Pattern_Type']
            never_fc = row['Never_log2FC']
            ex_fc = row['Ex_log2FC']
            cur_fc = row['Cur_log2FC']
            logger.info(f"  {gene} ({pattern}): Never log2FC = {never_fc:.2f}, Ex log2FC = {ex_fc:.2f}, Cur log2FC = {cur_fc:.2f}")
    
    return all_pattern_switches

def analyze_celltype_opposite_regulation():
    """
    Analyze cell type-specific genes with opposite regulation patterns between smoking categories.
    This function identifies genes within each cell type that show opposite regulation patterns
    between different smoking categories.
    """
    logger.info("\n=== ANALYZING CELL TYPE-SPECIFIC OPPOSITE REGULATION PATTERNS ===\n")
    
    cell_types = get_cell_type_list()
    all_celltype_switches = []
    
    for cell_type in cell_types:
        logger.info(f"Processing {cell_type}...")
        
        # Create dictionaries to store gene info for each smoking category
        never_genes = {}
        ex_genes = {}
        cur_genes = {}
        
        # Load data for each smoking category
        for smoking in smoking_categories:
            data = load_celltype_data(cell_type, smoking)
            if data is None:
                continue
                
            # Create dictionary with gene name as key, and log2fc and p-value as values
            gene_dict = {}
            for _, row in data.iterrows():
                gene_name = row['names']
                log2fc = row['logfoldchanges']
                pval = row['pvals_adj']
                gene_dict[gene_name] = {'log2fc': log2fc, 'padj': pval}
            
            # Store in the appropriate dictionary
            if smoking == 'Never':
                never_genes = gene_dict
            elif smoking == 'Ex':
                ex_genes = gene_dict
            elif smoking == 'Cur':
                cur_genes = gene_dict
        
        # Find genes with opposite regulation
        celltype_switches = []
        
        # Compare Never vs Ex
        common_genes_never_ex = set(never_genes.keys()) & set(ex_genes.keys())
        for gene in common_genes_never_ex:
            never_fc = never_genes[gene]['log2fc']
            never_padj = never_genes[gene]['padj']
            ex_fc = ex_genes[gene]['log2fc']
            ex_padj = ex_genes[gene]['padj']
            
            # Check for opposite regulation
            if ((never_fc > 0 and ex_fc < 0) or (never_fc < 0 and ex_fc > 0)) and \
               never_padj < 0.05 and ex_padj < 0.05:
                celltype_switches.append({
                    'Gene': gene,
                    'Cell_Type': cell_type,
                    'Never_log2FC': never_fc,
                    'Ex_log2FC': ex_fc,
                    'Cur_log2FC': cur_genes.get(gene, {}).get('log2fc', 0),
                    'Never_padj': never_padj,
                    'Ex_padj': ex_padj,
                    'Cur_padj': cur_genes.get(gene, {}).get('padj', 1),
                    'Comparison': 'Never_vs_Ex',
                    'Magnitude': abs(never_fc - ex_fc)
                })
        
        # Compare Never vs Cur
        common_genes_never_cur = set(never_genes.keys()) & set(cur_genes.keys())
        for gene in common_genes_never_cur:
            never_fc = never_genes[gene]['log2fc']
            never_padj = never_genes[gene]['padj']
            cur_fc = cur_genes[gene]['log2fc']
            cur_padj = cur_genes[gene]['padj']
            
            # Check for opposite regulation
            if ((never_fc > 0 and cur_fc < 0) or (never_fc < 0 and cur_fc > 0)) and \
               never_padj < 0.05 and cur_padj < 0.05:
                celltype_switches.append({
                    'Gene': gene,
                    'Cell_Type': cell_type,
                    'Never_log2FC': never_fc,
                    'Ex_log2FC': ex_genes.get(gene, {}).get('log2fc', 0),
                    'Cur_log2FC': cur_fc,
                    'Never_padj': never_padj,
                    'Ex_padj': ex_genes.get(gene, {}).get('padj', 1),
                    'Cur_padj': cur_padj,
                    'Comparison': 'Never_vs_Cur',
                    'Magnitude': abs(never_fc - cur_fc)
                })
        
        # Compare Ex vs Cur
        common_genes_ex_cur = set(ex_genes.keys()) & set(cur_genes.keys())
        for gene in common_genes_ex_cur:
            ex_fc = ex_genes[gene]['log2fc']
            ex_padj = ex_genes[gene]['padj']
            cur_fc = cur_genes[gene]['log2fc']
            cur_padj = cur_genes[gene]['padj']
            
            # Check for opposite regulation
            if ((ex_fc > 0 and cur_fc < 0) or (ex_fc < 0 and cur_fc > 0)) and \
               ex_padj < 0.05 and cur_padj < 0.05:
                celltype_switches.append({
                    'Gene': gene,
                    'Cell_Type': cell_type,
                    'Never_log2FC': never_genes.get(gene, {}).get('log2fc', 0),
                    'Ex_log2FC': ex_fc,
                    'Cur_log2FC': cur_fc,
                    'Never_padj': never_genes.get(gene, {}).get('padj', 1),
                    'Ex_padj': ex_padj,
                    'Cur_padj': cur_padj,
                    'Comparison': 'Ex_vs_Cur',
                    'Magnitude': abs(ex_fc - cur_fc)
                })
        
        # Log summary for current cell type
        if celltype_switches:
            logger.info(f"  Found {len(celltype_switches)} genes with opposite regulation in {cell_type}")
            
            # Sort by magnitude of change and log top genes
            sorted_switches = sorted(celltype_switches, key=lambda x: x['Magnitude'], reverse=True)
            top_switches = sorted_switches[:5]
            
            logger.info(f"  Top genes with opposite regulation in {cell_type}:")
            for switch in top_switches:
                gene = switch['Gene']
                comp = switch['Comparison']
                never_fc = switch['Never_log2FC']
                ex_fc = switch['Ex_log2FC']
                cur_fc = switch['Cur_log2FC']
                logger.info(f"    {gene} ({comp}): Never log2FC = {never_fc:.2f}, Ex log2FC = {ex_fc:.2f}, Cur log2FC = {cur_fc:.2f}")
            
            # Add to master list
            all_celltype_switches.extend(celltype_switches)
        else:
            logger.info(f"  No genes with opposite regulation found in {cell_type}")
    
    # Save all cell type-specific direction switch genes to CSV
    if all_celltype_switches:
        df = pd.DataFrame(all_celltype_switches)
        output_file = os.path.join(output_dir, "celltype_smoking_direction_switch_genes.csv")
        df.to_csv(output_file, index=False)
        logger.info(f"Saved {len(all_celltype_switches)} cell type-specific direction switch genes to {output_file}")
        
        # Create a summary by cell type
        summary = df.groupby('Cell_Type').size().reset_index(name='Count')
        summary_file = os.path.join(output_dir, "celltype_direction_switch_summary.csv")
        summary.to_csv(summary_file, index=False)
        logger.info(f"Saved cell type direction switch summary to {summary_file}")
    
    return all_celltype_switches

def create_comprehensive_report():
    """
    Create a comprehensive report with detailed information for LLM analysis.
    This function compiles all the analysis results into a detailed report that
    can be used by an LLM to generate comprehensive explanations and insights.
    """
    logger.info("\n=== CREATING COMPREHENSIVE REPORT FOR LLM ANALYSIS ===\n")
    
    # Create the report dictionary
    report = {
        "metadata": {
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "analysis_type": "Smoking-Tumor Interaction Analysis",
            "data_source": smoking_dir
        },
        "summary_statistics": {},
        "opposite_regulation": {},
        "cell_type_analysis": {},
        "top_genes": {}
    }
    
    # Add summary statistics
    celltype_counts_df = analyze_cell_type_deg_counts()
    
    # Calculate overall stats
    total_unique_cell_types = celltype_counts_df['Cell_Type'].nunique()
    smoking_summary = celltype_counts_df.groupby('Smoking_Status').agg({
        'Up_Regulated': 'sum',
        'Down_Regulated': 'sum',
        'Total_DEGs': 'sum'
    }).reset_index()
    
    report["summary_statistics"] = {
        "total_cell_types_analyzed": total_unique_cell_types,
        "overall_degs_by_smoking": smoking_summary.to_dict(orient='records')
    }
    
    # Add opposite regulation analysis
    all_switches = identify_opposite_regulation_genes()
    if all_switches:
        switch_df = pd.DataFrame(all_switches)
        
        # Get top genes by magnitude
        switch_df['magnitude'] = abs(switch_df['Never_log2FC'] - switch_df['Cur_log2FC'])
        top_switches = switch_df.sort_values('magnitude', ascending=False).head(50)
        
        report["opposite_regulation"] = {
            "total_genes_with_opposite_regulation": len(switch_df),
            "never_vs_ex_count": len(switch_df[switch_df['Pattern_Type'] == 'Never_vs_Ex']),
            "never_vs_cur_count": len(switch_df[switch_df['Pattern_Type'] == 'Never_vs_Cur']),
            "ex_vs_cur_count": len(switch_df[switch_df['Pattern_Type'] == 'Ex_vs_Cur']),
            "top_genes_with_opposite_regulation": top_switches.to_dict(orient='records')
        }
    
    # Add cell type-specific analysis
    celltype_switches = analyze_celltype_opposite_regulation()
    if celltype_switches:
        celltype_switch_df = pd.DataFrame(celltype_switches)
        
        # Group by cell type
        celltype_summary = celltype_switch_df.groupby('Cell_Type').agg({
            'Gene': 'count',
            'Magnitude': 'mean'
        }).reset_index()
        celltype_summary.columns = ['Cell_Type', 'Gene_Count', 'Average_Magnitude']
        
        # For each cell type, get top 5 genes
        top_genes_by_celltype = {}
        for cell_type in celltype_switch_df['Cell_Type'].unique():
            cell_data = celltype_switch_df[celltype_switch_df['Cell_Type'] == cell_type]
            top_5 = cell_data.sort_values('Magnitude', ascending=False).head(5).to_dict(orient='records')
            top_genes_by_celltype[cell_type] = top_5
        
        report["cell_type_analysis"] = {
            "cell_type_summary": celltype_summary.to_dict(orient='records'),
            "top_genes_by_cell_type": top_genes_by_celltype
        }
    
    # Add analysis of immune and inflammatory genes
    immune_genes = ['GBP1', 'GBP4', 'GBP5', 'SERPING1', 'IRF1', 'AKR1C3', 'CD38', 'ENG']
    immune_gene_data = {}
    
    # Load data
    comp_data = load_smoking_comparison_data()
    
    # Extract data for immune genes
    for gene in immune_genes:
        if gene in comp_data.index:
            row = comp_data.loc[gene]
            immune_gene_data[gene] = {
                'Never_log2FC': row.get('Never_log2FC', 0),
                'Ex_log2FC': row.get('Ex_log2FC', 0),
                'Cur_log2FC': row.get('Cur_log2FC', 0),
                'Never_padj': row.get('Never_padj', 1),
                'Ex_padj': row.get('Ex_padj', 1),
                'Cur_padj': row.get('Cur_padj', 1)
            }
    
    report["immune_inflammatory_genes"] = immune_gene_data
    
    # Save the comprehensive report
    report_file = os.path.join(output_dir, "smoking_tumor_comprehensive_report.json")
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"Saved comprehensive report to {report_file}")
    
    # Also save a CSV version of key components for easy viewing
    if all_switches:
        top_switches_file = os.path.join(output_dir, "top_smoking_direction_switch_genes.csv")
        top_switches.to_csv(top_switches_file, index=False)
        logger.info(f"Saved top direction switch genes to {top_switches_file}")
    
    if immune_gene_data:
        immune_df = pd.DataFrame.from_dict(immune_gene_data, orient='index')
        immune_file = os.path.join(output_dir, "immune_inflammatory_gene_expression.csv")
        immune_df.to_csv(immune_file)
        logger.info(f"Saved immune/inflammatory gene expression data to {immune_file}")
    
    return report

def main():
    """Main execution function."""
    logger.info("=== STARTING SMOKING-TUMOR INTERACTION ANALYSIS ===")
    logger.info(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Data directory: {smoking_dir}")
    logger.info(f"Output directory: {output_dir}")
    
    # Analyze cell type DEG counts
    analyze_cell_type_deg_counts()
    
    # Identify genes with opposite regulation patterns
    identify_opposite_regulation_genes()
    
    # Analyze cell type-specific opposite regulation
    analyze_celltype_opposite_regulation()
    
    # Create comprehensive report
    create_comprehensive_report()
    
    logger.info("\n=== ANALYSIS COMPLETE ===")
    logger.info(f"Results saved to {output_dir}")

if __name__ == "__main__":
    main()