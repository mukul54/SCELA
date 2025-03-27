import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from collections import Counter

# File paths
base_dir = "/home/mukul.ranjan/Documents/cancer/SCELA/age_analysis/tumor_vs_normal"
files = {
    "global": os.path.join(base_dir, "global_tumor_vs_normal_significant.csv"),
    "40-54": os.path.join(base_dir, "40-54_(Mid-Adult)/tumor_vs_normal_significant.csv"),
    "55-64": os.path.join(base_dir, "55-64_(Pre-Senior)/tumor_vs_normal_significant.csv"),
    "65-74": os.path.join(base_dir, "65-74_(Young_Elderly)/tumor_vs_normal_significant.csv"),
    "75-84": os.path.join(base_dir, "75-84_(Senior)/tumor_vs_normal_significant.csv")
}

# Read the cell type percentages
celltype_pct = pd.read_csv("/home/mukul.ranjan/Documents/cancer/SCELA/age_analysis/additional_plots/celltype_percent_by_age.csv")
celltype_pct = celltype_pct.set_index('Age_Group')

# Function to load and analyze DEG files
def load_degs(file_path):
    try:
        df = pd.read_csv(file_path)
        # Classify genes as up or down regulated in tumor
        df['regulation'] = df['logfoldchanges'].apply(lambda x: 'Up in Tumor' if x > 0 else 'Up in Normal')
        return df
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return pd.DataFrame()

# Load all DEG datasets
deg_data = {age: load_degs(path) for age, path in files.items()}

# 1. Basic stats for each age group
def get_deg_stats(data_dict):
    stats = {}
    for age, df in data_dict.items():
        if not df.empty:
            up_in_tumor = sum(df['logfoldchanges'] > 0)
            up_in_normal = sum(df['logfoldchanges'] < 0)
            stats[age] = {
                'total_degs': len(df),
                'up_in_tumor': up_in_tumor,
                'up_in_normal': up_in_normal,
                'ratio': up_in_tumor / up_in_normal if up_in_normal > 0 else float('inf')
            }
    return pd.DataFrame(stats).T

# 2. Find age-specific DEGs (unique to each age group)
def find_age_specific_genes(data_dict):
    all_genes = {age: set(df['names'].unique()) for age, df in data_dict.items() if age != 'global' and not df.empty}
    
    age_specific = {}
    for age, genes in all_genes.items():
        other_ages_genes = set()
        for other_age, other_genes in all_genes.items():
            if other_age != age:
                other_ages_genes.update(other_genes)
        
        unique_genes = genes - other_ages_genes
        age_specific[age] = unique_genes
    
    return age_specific

# 3. Find genes with age-dependent expression patterns
def find_age_dependent_genes(data_dict):
    # Genes present in at least 3 age groups
    age_groups = [k for k in data_dict.keys() if k != 'global']
    gene_presence = Counter()
    
    for age in age_groups:
        if not data_dict[age].empty:
            gene_presence.update(data_dict[age]['names'].unique())
    
    common_genes = [gene for gene, count in gene_presence.items() if count >= 3]
    
    # Check for trends across age groups
    trends = {}
    for gene in common_genes:
        fold_changes = {}
        for age in age_groups:
            if not data_dict[age].empty and gene in data_dict[age]['names'].values:
                fold_changes[age] = data_dict[age][data_dict[age]['names'] == gene]['logfoldchanges'].values[0]
        
        if len(fold_changes) >= 3:
            trends[gene] = fold_changes
    
    return trends

# 4. Find top DEGs for each age group
def get_top_degs(data_dict, n=20):
    top_genes = {}
    for age, df in data_dict.items():
        if not df.empty:
            # Sort by absolute fold change
            sorted_df = df.sort_values(by='logfoldchanges', key=abs, ascending=False)
            top_genes[age] = sorted_df.head(n)
    return top_genes

# 5. Identify cell types with U-shaped and other distribution patterns
def identify_distribution_patterns(celltype_df):
    """Identify cell types with different distribution patterns across age groups"""
    # Define age group order for analysis
    age_order = ['<40 (Young)', '40-54 (Mid-Adult)', '55-64 (Pre-Senior)', '65-74 (Young Elderly)', '75-84 (Senior)']
    
    # Filter to only include rows in the specified order
    cell_df = celltype_df.loc[celltype_df.index.intersection(age_order)]
    if len(cell_df) < 3:
        print("Not enough age groups for pattern analysis")
        return {}
    
    # Reindex to ensure proper order
    cell_df = cell_df.reindex(age_order, axis=0)
    
    # Initialize pattern dictionaries
    patterns = {
        'u_shaped': [],
        'inverse_u_shaped': [],
        'increasing': [],
        'decreasing': []
    }
    
    # Check each cell type for patterns
    for col in cell_df.columns:
        values = cell_df[col].values
        
        # Skip if values contain NaN
        if np.isnan(values).any():
            continue
        
        # Need at least 3 points to determine pattern
        if len(values) < 3:
            continue
        
        # Check for U-shaped pattern (first and last values higher than middle)
        middle_values = values[1:-1]
        if values[0] > np.mean(middle_values) and values[-1] > np.mean(middle_values):
            patterns['u_shaped'].append(col)
        
        # Check for inverse U-shaped pattern
        elif values[0] < np.mean(middle_values) and values[-1] < np.mean(middle_values):
            patterns['inverse_u_shaped'].append(col)
        
        # Check for increasing trend
        elif np.corrcoef(range(len(values)), values)[0, 1] > 0.7:
            patterns['increasing'].append(col)
        
        # Check for decreasing trend
        elif np.corrcoef(range(len(values)), values)[0, 1] < -0.7:
            patterns['decreasing'].append(col)
    
    return patterns

# 6. Plot cell type distribution patterns
def plot_distribution_patterns(celltype_df, patterns, output_dir):
    """Create plots for different cell type distribution patterns"""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Define age group order for analysis
    age_order = ['<40 (Young)', '40-54 (Mid-Adult)', '55-64 (Pre-Senior)', '65-74 (Young Elderly)', '75-84 (Senior)']
    
    # Filter to only include rows in the specified order
    cell_df = celltype_df.loc[celltype_df.index.intersection(age_order)]
    
    # Reindex to ensure proper order
    cell_df = cell_df.reindex(age_order, axis=0)
    
    # Create plots for each pattern type
    pattern_types = {
        'u_shaped': 'U-Shaped Distribution',
        'inverse_u_shaped': 'Inverse U-Shaped Distribution',
        'increasing': 'Increasing with Age',
        'decreasing': 'Decreasing with Age'
    }
    
    for pattern_type, title in pattern_types.items():
        cell_types = patterns.get(pattern_type, [])
        if not cell_types:
            continue
        
        # Limit to 6 cell types per plot for readability
        for i in range(0, len(cell_types), 6):
            subset = cell_types[i:i+6]
            
            plt.figure(figsize=(12, 8))
            
            for cell_type in subset:
                plt.plot(range(len(cell_df)), cell_df[cell_type].values, marker='o', linewidth=2, label=cell_type)
            
            plt.title(f"{title} Cell Types", fontsize=16)
            plt.xlabel('Age Group', fontsize=14)
            plt.ylabel('Percentage of Cells', fontsize=14)
            plt.xticks(range(len(cell_df)), cell_df.index, rotation=45)
            plt.legend(loc='best')
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            
            # Save the plot
            idx = i // 6 if i > 0 else ''
            plt.savefig(os.path.join(output_dir, f"{pattern_type}_cell_types{idx}.pdf"))
            plt.close()

# 7. Correlation between cell type percentages and DEG counts
def correlate_celltypes_with_degs(celltype_df, deg_stats):
    # Ensure we have matching age groups
    matching_ages = []
    for age_deg in deg_stats.index:
        for age_cell in celltype_df.index:
            if age_deg.split('_')[0] in age_cell or age_deg in age_cell:
                matching_ages.append((age_deg, age_cell))
    
    correlations = {}
    for cell_type in celltype_df.columns:
        cell_pcts = []
        deg_counts = []
        
        for age_deg, age_cell in matching_ages:
            cell_pcts.append(celltype_df.loc[age_cell, cell_type])
            deg_counts.append(deg_stats.loc[age_deg, 'total_degs'])
        
        if len(cell_pcts) >= 3:  # Need at least 3 points for meaningful correlation
            corr = np.corrcoef(cell_pcts, deg_counts)[0, 1]
            correlations[cell_type] = corr
    
    return pd.Series(correlations).sort_values(ascending=False)

# Run analyses
print("Computing DEG statistics...")
deg_stats = get_deg_stats(deg_data)
print("\nDEG Statistics by Age Group:")
print(deg_stats)

print("\nFinding age-specific genes...")
age_specific = find_age_specific_genes(deg_data)
for age, genes in age_specific.items():
    print(f"\n{age} specific genes ({len(genes)}):")
    print(", ".join(list(genes)[:20]) + ("..." if len(genes) > 20 else ""))

print("\nIdentifying top DEGs for each age group...")
top_degs = get_top_degs(deg_data)
for age, genes_df in top_degs.items():
    if not genes_df.empty:
        print(f"\nTop 20 DEGs for {age}:")
        print(genes_df[['names', 'logfoldchanges', 'pvals_adj']].to_string(index=False))

print("\nFinding genes with age-dependent expression patterns...")
age_dependent = find_age_dependent_genes(deg_data)
print(f"Found {len(age_dependent)} genes with age-dependent expression patterns")
if age_dependent:
    # Show a few examples
    examples = list(age_dependent.items())[:5]
    for gene, fold_changes in examples:
        print(f"\n{gene}: {fold_changes}")

# Correlate cell type percentages with DEG counts
print("\nCorrelating cell type percentages with DEG counts...")
correlations = correlate_celltypes_with_degs(celltype_pct, deg_stats)
print("\nTop cell types positively correlated with DEG counts:")
print(correlations.head(10))

# Print cell types with strongest correlation to DEG counts
print("\nCell types most strongly associated with tumor-normal differences across ages:")
for cell_type, corr in correlations.head(5).items():
    print(f"{cell_type}: {corr:.4f}")

# Identify cell types with U-shaped and other distributions
print("\nAnalyzing cell type distribution patterns across age groups...")
patterns = identify_distribution_patterns(celltype_pct)

print("\nCell types with U-shaped distribution:")
print(', '.join(patterns['u_shaped']))

print("\nCell types with inverse U-shaped distribution:")
print(', '.join(patterns['inverse_u_shaped']))

print("\nCell types with increasing trend with age:")
print(', '.join(patterns['increasing']))

print("\nCell types with decreasing trend with age:")
print(', '.join(patterns['decreasing']))

# Create plots of distribution patterns
output_dir = os.path.join(os.path.dirname(base_dir), "additional_plots")
plot_distribution_patterns(celltype_pct, patterns, output_dir)
print(f"\nDistribution pattern plots saved to {output_dir}")

# Optional - create visualizations
plt.figure(figsize=(12, 8))
up_in_tumor = [stats['up_in_tumor'] for stats in deg_stats.to_dict('records')]
up_in_normal = [stats['up_in_normal'] for stats in deg_stats.to_dict('records')]
x = np.arange(len(deg_stats))
width = 0.35

plt.bar(x - width/2, up_in_tumor, width, label='Up in Tumor')
plt.bar(x + width/2, up_in_normal, width, label='Up in Normal')
plt.xlabel('Age Group')
plt.ylabel('Number of DEGs')
plt.title('Tumor vs Normal DEGs by Age Group')
plt.xticks(x, deg_stats.index)
plt.legend()
plt.tight_layout()
plt.savefig("age_tumor_normal_degs.pdf")