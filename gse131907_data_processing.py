#!/usr/bin/env python
# Convert GSE131907 data to HDF5 and prepare for cell type annotation analysis

import os
import re
import pandas as pd
import numpy as np
import h5py
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse
import tables
from time import time
import gc
import warnings
warnings.filterwarnings('ignore')

# File Paths
BASE_DATA_PATH = "/l/users/mukul.ranjan/cancer/data/GSE131907"
RAW_UMI_FILE = f"{BASE_DATA_PATH}/GSE131907_Lung_Cancer_raw_UMI_matrix.txt"
CELL_ANNOTATION_FILE = f"{BASE_DATA_PATH}/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
METADATA_FILE = f"{BASE_DATA_PATH}/GSE131907_metadata.csv"
OUTPUT_H5_FILE = f"{BASE_DATA_PATH}/GSE131907_processed.h5"
OUTPUT_H5AD_FILE = f"{BASE_DATA_PATH}/GSE131907_anndata.h5ad"

def convert_to_hdf5():
    """
    Convert the large UMI matrix to HDF5 format using chunking to manage memory
    """
    print("Starting conversion to HDF5 format...")
    start_time = time()
    
    # Check if output file already exists to avoid redundant processing
    if os.path.exists(OUTPUT_H5_FILE):
        print(f"HDF5 file already exists at {OUTPUT_H5_FILE}")
        return
    
    # First, read the header to get all column names
    print("Reading header...")
    with open(RAW_UMI_FILE, 'r') as f:
        header = f.readline().strip().split('\t')
    
    gene_column = header[0]
    sample_columns = header[1:]
    
    # Create HDF5 file
    print("Creating HDF5 file...")
    with h5py.File(OUTPUT_H5_FILE, 'w') as h5f:
        # Store column names
        h5f.create_dataset('gene_column', data=np.array([gene_column], dtype='S100'))
        h5f.create_dataset('sample_columns', data=np.array(sample_columns, dtype='S100'))
        
        # Process the data in chunks
        chunk_size = 1000  # Number of genes per chunk
        chunk_idx = 0
        genes = []
        chunk_data = []
        
        # Use float32 to save memory
        # (No longer using tables package)
        
        print(f"Processing data in chunks of {chunk_size} genes...")
        # Process file in chunks
        with open(RAW_UMI_FILE, 'r') as f:
            f.readline()  # Skip header
            
            for i, line in enumerate(f):
                if i % 1000 == 0:
                    print(f"Processed {i} genes...")
                
                values = line.strip().split('\t')
                gene_name = values[0]
                genes.append(gene_name)
                expression_values = np.array(values[1:], dtype=np.float32)
                
                chunk_data.append(expression_values)
                
                # Write chunk to HDF5 when it reaches chunk_size
                if len(chunk_data) >= chunk_size:
                    chunk_array = np.array(chunk_data)
                    dataset_name = f'chunk_{chunk_idx}'
                    h5f.create_dataset(dataset_name, data=chunk_array, 
                                      chunks=True, compression='gzip', 
                                      compression_opts=4)
                    
                    # Store the gene names for this chunk
                    gene_names_dataset = f'genes_chunk_{chunk_idx}'
                    h5f.create_dataset(gene_names_dataset, 
                                      data=np.array(genes, dtype='S100'))
                    
                    # Reset for next chunk
                    chunk_data = []
                    genes = []
                    chunk_idx += 1
                    gc.collect()  # Force garbage collection
                    
            # Don't forget the last chunk if there's anything left
            if len(chunk_data) > 0:
                chunk_array = np.array(chunk_data)
                dataset_name = f'chunk_{chunk_idx}'
                h5f.create_dataset(dataset_name, data=chunk_array, 
                                  chunks=True, compression='gzip', 
                                  compression_opts=4)
                
                # Store the gene names for this chunk
                gene_names_dataset = f'genes_chunk_{chunk_idx}'
                h5f.create_dataset(gene_names_dataset, 
                                  data=np.array(genes, dtype='S100'))
    
    elapsed_time = time() - start_time
    print(f"Conversion completed in {elapsed_time:.2f} seconds")
    print(f"HDF5 file created at: {OUTPUT_H5_FILE}")

def create_anndata_object(force_recreate_ann=False):
    """
    Create an AnnData object combining UMI data with metadata
    """
    print("Creating AnnData object...")
    start_time = time()
    
    # Check if output file already exists
    if os.path.exists(OUTPUT_H5AD_FILE) and not force_recreate_ann:
        print(f"AnnData file already exists at {OUTPUT_H5AD_FILE}")
        return sc.read_h5ad(OUTPUT_H5AD_FILE)
    
    # Load cell annotations
    print("Loading cell annotations...")
    cell_annotations = pd.read_csv(CELL_ANNOTATION_FILE, sep='\t', 
                                 compression='gzip', index_col=0)
    
    # Load metadata
    print("Loading metadata...")
    metadata = pd.read_csv(METADATA_FILE, index_col=0)
    
    # Clean and validate metadata columns
    required_columns = ['Samples', 'Sex', 'Age', 'Smoking', 'Stages', 'EGFR', 'Histology']
    missing_cols = [col for col in required_columns if col not in metadata.columns]
    if missing_cols:
        raise ValueError(f"Metadata missing required columns: {missing_cols}")
    
    # Convert smoking status to categorical
    metadata['Smoking'] = metadata['Smoking'].astype('category')
    
    # Clean stage information
    metadata['Stages'] = metadata['Stages'].str.replace(r'[^IVABC]', '', regex=True)
    
    # Get sample information from column names using the HDF5 file
    print("Extracting sample information...")
    with h5py.File(OUTPUT_H5_FILE, 'r') as h5f:
        # Verify that the file has the expected structure
        if 'sample_columns' not in h5f:
            print("Error: HDF5 file doesn't have the expected structure.")
            print("Re-creating the HDF5 file from scratch...")
            os.remove(OUTPUT_H5_FILE)
            convert_to_hdf5()
            # Reopen the file
            h5f = h5py.File(OUTPUT_H5_FILE, 'r')
            
        sample_columns = h5f['sample_columns'][:]
        sample_columns = [col.decode('utf-8') for col in sample_columns]
        
        # Find the first available chunk
        chunk_found = False
        first_chunk = None
        num_features = 0
        
        # Count total number of features (genes)
        for key in h5f.keys():
            if key.startswith('chunk_'):
                if not chunk_found:
                    first_chunk = h5f[key][:]
                    chunk_found = True
                
            if key.startswith('genes_chunk_'):
                num_features += len(h5f[key])
        
        if not chunk_found:
            print("Error: No data chunks found in HDF5 file.")
            print("The file may be corrupted. Re-creating from scratch...")
            os.remove(OUTPUT_H5_FILE)
            convert_to_hdf5()
            return create_anndata_object(force_recreate_ann=True)
        
        # Create a sparse matrix to save memory
        print(f"Creating sparse matrix with {num_features} genes and {len(sample_columns)} cells...")
        X = sparse.lil_matrix((num_features, len(sample_columns)), dtype=np.float32)
        
        # Fill the matrix with data from chunks
        row_idx = 0
        all_genes = []
        
        for chunk_idx in range(100):  # Assuming less than 100 chunks
            chunk_key = f'chunk_{chunk_idx}'
            gene_key = f'genes_chunk_{chunk_idx}'
            
            if chunk_key not in h5f:
                break
                
            chunk_data = h5f[chunk_key][:]
            chunk_genes = h5f[gene_key][:]
            chunk_genes = [gene.decode('utf-8') for gene in chunk_genes]
            
            num_genes_in_chunk = len(chunk_genes)
            X[row_idx:row_idx+num_genes_in_chunk, :] = chunk_data
            all_genes.extend(chunk_genes)
            
            row_idx += num_genes_in_chunk
    
    # Create observation dataframe with proper sample IDs
    obs_df = pd.DataFrame({'Samples': sample_columns})
    
    # Clean sample IDs by extracting core identifier
    obs_df['Samples'] = obs_df['Samples'].str.extract(r'(LUNG_[NT]\d+|LN_\d+|NS_\d+|EBUS_\d+|BRONCHO_\d+|EFFUSION_\d+)')[0]
    
    # Merge metadata using pandas merge
    obs_df = pd.merge(obs_df, metadata, on='Samples', how='left')
    obs_df.set_index(pd.Index(sample_columns), inplace=True)

    # Create var DataFrame (gene annotations)
    var_df = pd.DataFrame(index=all_genes)
    
    # Convert to compressed sparse row format for efficiency
    X_csr = X.tocsr()
    
    # Create AnnData object
    print("Creating final AnnData object...")
    adata = sc.AnnData(X=X_csr.T, obs=obs_df, var=var_df)  # Transpose to get cells as rows
    
    # Add cell annotations if available
    if cell_annotations is not None:
        common_cells = adata.obs_names.intersection(cell_annotations.index)
        if len(common_cells) > 0:
            for col in cell_annotations.columns:
                adata.obs[col] = [
                    cell_annotations.loc[cell, col] if cell in cell_annotations.index else np.nan
                    for cell in adata.obs_names
                ]
            print(f"Added cell annotations for {len(common_cells)} cells")
    
    # Save the AnnData object
    print(f"Saving AnnData object to {OUTPUT_H5AD_FILE}...")
    adata.write(OUTPUT_H5AD_FILE)
    
    elapsed_time = time() - start_time
    print(f"AnnData creation completed in {elapsed_time:.2f} seconds")
    
    return adata

def perform_basic_analysis(adata=None):
    """
    Perform basic analysis on the data for cell type annotation
    """
    if adata is None:
        # Load the AnnData object if not provided
        if os.path.exists(OUTPUT_H5AD_FILE):
            print(f"Loading AnnData from {OUTPUT_H5AD_FILE}...")
            adata = sc.read_h5ad(OUTPUT_H5AD_FILE)
        else:
            print("AnnData file not found. Creating it first...")
            adata = create_anndata_object(force_recreate_ann=True)
    
    print("\n=== Basic Analysis of GSE131907 Lung Cancer Dataset ===")
    print(f"Dataset dimensions: {adata.shape[0]} cells × {adata.shape[1]} genes")
    
    # Basic preprocessing
    print("\nPerforming basic preprocessing...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Calculate QC metrics
    print("Calculating QC metrics...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Normalize the data
    print("Normalizing data...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Identify highly variable genes
    print("Identifying highly variable genes...")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    # PCA and UMAP
    print("Running dimensionality reduction...")
    sc.pp.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    
    # Clustering
    print("Performing clustering...")
    sc.tl.leiden(adata, resolution=0.5)
    
    # Find marker genes
    print("Finding marker genes...")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    
    # Save the processed data
    print("Saving processed data...")
    adata.write(f"{BASE_DATA_PATH}/GSE131907_processed_analysis.h5ad")
    
    return adata

def main():
    """
    Main function to run the entire pipeline
    """
    print("\n=== GSE131907 Lung Cancer Dataset Processing Pipeline ===\n")
    
    # Step 1: Convert data to HDF5 format
    if not os.path.exists(OUTPUT_H5_FILE):
        convert_to_hdf5()
    else:
        print(f"HDF5 file already exists at {OUTPUT_H5_FILE}")
    force_recreate_ann = True
    # Step 2: Create AnnData object
    if not os.path.exists(OUTPUT_H5AD_FILE):
        adata = create_anndata_object()
    elif force_recreate_ann:
        print(f"force_recreate_ann flag is enabled recreating it")
        adata = create_anndata_object(force_recreate_ann)
    else:
        print(f"AnnData file already exists at {OUTPUT_H5AD_FILE}")
        adata = sc.read_h5ad(OUTPUT_H5AD_FILE)
    
    # Step 3: Perform basic analysis
    perform_basic_analysis(adata)
    
    print("\n=== Pipeline completed successfully! ===")
    print("\nNext steps for analysis:")
    print("1. Identify marker genes for each cluster")
    print("2. Annotate cell types based on established markers")
    print("3. Compare cell type composition between demographic groups")
    print("4. Create visualizations of cell type distributions")

if __name__ == "__main__":
    main()
