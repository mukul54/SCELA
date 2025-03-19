"""
Configuration settings for the cancer analysis pipeline
"""

import os

# Data paths
DATA_DIR = "/l/users/mukul.ranjan/cancer/data/GSE131907"
RESULTS_DIR = os.path.join(DATA_DIR, "analysis_results")
VIS_DIR = os.path.join(RESULTS_DIR, "visualizations")

# Analysis parameters
N_PCS = 50           # Number of principal components to use
N_NEIGHBORS = 15     # Number of neighbors for neighborhood graph
RESOLUTION = 0.5     # Default resolution for Leiden clustering
RESOLUTIONS = [0.25, 0.5, 1.0]  # Multiple resolutions for cluster exploration
MARKER_GENE_RESOLUTIONS = [0.25]  # Resolutions for which to identify marker genes
DEFAULT_MARKER_RESOLUTION = 0.25  # Default resolution for marker gene identification
# Recomputation flags for different analysis steps
FORCE_RECOMPUTE_PREPROCESSING = True  # Force recomputation of preprocessing steps
FORCE_RECOMPUTE_CLUSTERING = True    # Force recomputation of clustering results
FORCE_RECOMPUTE_MARKERS = True       # Force recomputation of marker genes
FORCE_RECOMPUTE_DEMOGRAPHICS = True  # Force recomputation of demographic analysis
FORCE_RECOMPUTE_VISUALIZATION = True # Force regeneration of visualizations
N_TOP_GENES = 2000   # Number of highly variable genes to select
TOP_N_MARKERS = 5    # Number of top marker genes to show per cell type

# Visualization parameters
DYNAMIC_FIGURE_SIZE = True  # Whether to adjust figure sizes based on data dimensions
