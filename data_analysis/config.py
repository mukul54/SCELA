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
RESOLUTION = 0.5     # Resolution for Leiden clustering
N_TOP_GENES = 2000   # Number of highly variable genes to select
TOP_N_MARKERS = 5    # Number of top marker genes to show per cell type
