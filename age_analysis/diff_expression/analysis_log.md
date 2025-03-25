python age_specific_degs.py
Loading data from ./age_analysis/age_celltype_annotation.h5ad
Loaded 166804 cells with 27578 genes
Tissue type distribution:
Tissue_Type
Tumor     102513
Normal     64291
Name: count, dtype: int64

Analyzing cell type: T cytotoxic cells
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
/home/mukul.ranjan/.conda/envs/cancer/lib/python3.11/site-packages/numba/np/ufunc/parallel.py:371: NumbaWarning: The TBB threading layer requires TBB version 2021 update 6 or later i.e., TBB_INTERFACE_VERSION >= 12060. Found TBB_INTERFACE_VERSION = 12050. The TBB threading layer is disabled.
  warnings.warn(problem)
WARNING: saving figure to file figures/violin_T cytotoxic cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_T cytotoxic cells_55-64_down.pdf
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_T cytotoxic cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_T cytotoxic cells_65-74_down.pdf
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_T cytotoxic cells_75-84_up.pdf
WARNING: saving figure to file figures/violin_T cytotoxic cells_75-84_down.pdf
  Comparing <40 (Young) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_T cytotoxic cells_<40_up.pdf
WARNING: saving figure to file figures/violin_T cytotoxic cells_<40_down.pdf

Analyzing cell type: Gamma delta T cells
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Gamma delta T cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_Gamma delta T cells_65-74_down.pdf
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Gamma delta T cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_Gamma delta T cells_55-64_down.pdf
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Gamma delta T cells_75-84_up.pdf
WARNING: saving figure to file figures/violin_Gamma delta T cells_75-84_down.pdf
  Skipping <40 (Young) - insufficient cells (ref: 3071, target: 11)

Analyzing cell type: Cholangiocytes
  Skipping 55-64 (Pre-Senior) - insufficient cells (ref: 0, target: 3662)

Analyzing cell type: T memory cells
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_T memory cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_T memory cells_65-74_down.pdf
/home/mukul.ranjan/.conda/envs/cancer/lib/python3.11/site-packages/pandas/core/arraylike.py:399: RuntimeWarning: divide by zero encountered in log10
  result = getattr(ufunc, method)(*inputs, **kwargs)
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_T memory cells_75-84_up.pdf
WARNING: saving figure to file figures/violin_T memory cells_75-84_down.pdf
/home/mukul.ranjan/.conda/envs/cancer/lib/python3.11/site-packages/pandas/core/arraylike.py:399: RuntimeWarning: divide by zero encountered in log10
  result = getattr(ufunc, method)(*inputs, **kwargs)
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_T memory cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_T memory cells_55-64_down.pdf
/home/mukul.ranjan/.conda/envs/cancer/lib/python3.11/site-packages/pandas/core/arraylike.py:399: RuntimeWarning: divide by zero encountered in log10
  result = getattr(ufunc, method)(*inputs, **kwargs)

Analyzing cell type: T cells
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_T cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_T cells_65-74_down.pdf
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_T cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_T cells_55-64_down.pdf
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_T cells_75-84_up.pdf
WARNING: saving figure to file figures/violin_T cells_75-84_down.pdf
  Comparing <40 (Young) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_T cells_<40_up.pdf
WARNING: saving figure to file figures/violin_T cells_<40_down.pdf

Analyzing cell type: Plasma cells
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Plasma cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_Plasma cells_55-64_down.pdf
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Plasma cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_Plasma cells_65-74_down.pdf
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Plasma cells_75-84_up.pdf
WARNING: saving figure to file figures/violin_Plasma cells_75-84_down.pdf
  Skipping <40 (Young) - insufficient cells (ref: 406, target: 1)

Analyzing cell type: Macrophages
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Macrophages_65-74_up.pdf
WARNING: saving figure to file figures/violin_Macrophages_65-74_down.pdf
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Macrophages_55-64_up.pdf
WARNING: saving figure to file figures/violin_Macrophages_55-64_down.pdf
/home/mukul.ranjan/.conda/envs/cancer/lib/python3.11/site-packages/pandas/core/arraylike.py:399: RuntimeWarning: divide by zero encountered in log10
  result = getattr(ufunc, method)(*inputs, **kwargs)
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Macrophages_75-84_up.pdf
WARNING: saving figure to file figures/violin_Macrophages_75-84_down.pdf
  Comparing <40 (Young) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Macrophages_<40_up.pdf
WARNING: saving figure to file figures/violin_Macrophages_<40_down.pdf

Analyzing cell type: B cells naive
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_B cells naive_55-64_up.pdf
WARNING: saving figure to file figures/violin_B cells naive_55-64_down.pdf
/home/mukul.ranjan/.conda/envs/cancer/lib/python3.11/site-packages/pandas/core/arraylike.py:399: RuntimeWarning: divide by zero encountered in log10
  result = getattr(ufunc, method)(*inputs, **kwargs)
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_B cells naive_65-74_up.pdf
WARNING: saving figure to file figures/violin_B cells naive_65-74_down.pdf
/home/mukul.ranjan/.conda/envs/cancer/lib/python3.11/site-packages/pandas/core/arraylike.py:399: RuntimeWarning: divide by zero encountered in log10
  result = getattr(ufunc, method)(*inputs, **kwargs)
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_B cells naive_75-84_up.pdf
WARNING: saving figure to file figures/violin_B cells naive_75-84_down.pdf
/home/mukul.ranjan/.conda/envs/cancer/lib/python3.11/site-packages/pandas/core/arraylike.py:399: RuntimeWarning: divide by zero encountered in log10
  result = getattr(ufunc, method)(*inputs, **kwargs)

Analyzing cell type: Endothelial cells
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Endothelial cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_Endothelial cells_55-64_down.pdf
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Endothelial cells_75-84_up.pdf
WARNING: saving figure to file figures/violin_Endothelial cells_75-84_down.pdf
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Endothelial cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_Endothelial cells_65-74_down.pdf
  Comparing <40 (Young) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Endothelial cells_<40_up.pdf
WARNING: saving figure to file figures/violin_Endothelial cells_<40_down.pdf

Analyzing cell type: Myoepithelial cells
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Myoepithelial cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_Myoepithelial cells_55-64_down.pdf
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Myoepithelial cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_Myoepithelial cells_65-74_down.pdf
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Myoepithelial cells_75-84_up.pdf
WARNING: saving figure to file figures/violin_Myoepithelial cells_75-84_down.pdf
  Skipping <40 (Young) - insufficient cells (ref: 1598, target: 7)

Analyzing cell type: Neutrophils
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Neutrophils_75-84_up.pdf
WARNING: saving figure to file figures/violin_Neutrophils_75-84_down.pdf
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Neutrophils_65-74_up.pdf
WARNING: saving figure to file figures/violin_Neutrophils_65-74_down.pdf
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Neutrophils_55-64_up.pdf
WARNING: saving figure to file figures/violin_Neutrophils_55-64_down.pdf
  Skipping <40 (Young) - insufficient cells (ref: 1073, target: 8)

Analyzing cell type: Mammary epithelial cells
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Mammary epithelial cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_Mammary epithelial cells_55-64_down.pdf
/home/mukul.ranjan/.conda/envs/cancer/lib/python3.11/site-packages/pandas/core/arraylike.py:399: RuntimeWarning: divide by zero encountered in log10
  result = getattr(ufunc, method)(*inputs, **kwargs)
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Mammary epithelial cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_Mammary epithelial cells_65-74_down.pdf
/home/mukul.ranjan/.conda/envs/cancer/lib/python3.11/site-packages/pandas/core/arraylike.py:399: RuntimeWarning: divide by zero encountered in log10
  result = getattr(ufunc, method)(*inputs, **kwargs)
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Mammary epithelial cells_75-84_up.pdf
WARNING: saving figure to file figures/violin_Mammary epithelial cells_75-84_down.pdf
/home/mukul.ranjan/.conda/envs/cancer/lib/python3.11/site-packages/pandas/core/arraylike.py:399: RuntimeWarning: divide by zero encountered in log10
  result = getattr(ufunc, method)(*inputs, **kwargs)
  Skipping <40 (Young) - insufficient cells (ref: 5514, target: 7)

Analyzing cell type: NK cells
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_NK cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_NK cells_65-74_down.pdf
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_NK cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_NK cells_55-64_down.pdf
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_NK cells_75-84_up.pdf
WARNING: saving figure to file figures/violin_NK cells_75-84_down.pdf
  Skipping <40 (Young) - insufficient cells (ref: 930, target: 14)

Analyzing cell type: Pulmonary alveolar type II cells
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Pulmonary alveolar type II cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_Pulmonary alveolar type II cells_65-74_down.pdf
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Pulmonary alveolar type II cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_Pulmonary alveolar type II cells_55-64_down.pdf
  Skipping 75-84 (Senior) - insufficient cells (ref: 1586, target: 8)

Analyzing cell type: Dendritic cells
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Dendritic cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_Dendritic cells_55-64_down.pdf
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Dendritic cells_75-84_up.pdf
WARNING: saving figure to file figures/violin_Dendritic cells_75-84_down.pdf
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Dendritic cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_Dendritic cells_65-74_down.pdf
  Comparing <40 (Young) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Dendritic cells_<40_up.pdf
WARNING: saving figure to file figures/violin_Dendritic cells_<40_down.pdf

Analyzing cell type: Epithelial cells
  Skipping <40 (Young) - insufficient cells (ref: 0, target: 990)
  Skipping 75-84 (Senior) - insufficient cells (ref: 0, target: 1)

Analyzing cell type: Plasmacytoid dendritic cells
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Plasmacytoid dendritic cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_Plasmacytoid dendritic cells_65-74_down.pdf
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Plasmacytoid dendritic cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_Plasmacytoid dendritic cells_55-64_down.pdf
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Plasmacytoid dendritic cells_75-84_up.pdf
WARNING: saving figure to file figures/violin_Plasmacytoid dendritic cells_75-84_down.pdf
  Skipping <40 (Young) - insufficient cells (ref: 58, target: 1)

Analyzing cell type: Hepatic stellate cells
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Hepatic stellate cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_Hepatic stellate cells_65-74_down.pdf
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Hepatic stellate cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_Hepatic stellate cells_55-64_down.pdf
  Comparing <40 (Young) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Hepatic stellate cells_<40_up.pdf
WARNING: saving figure to file figures/violin_Hepatic stellate cells_<40_down.pdf
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Hepatic stellate cells_75-84_up.pdf
WARNING: saving figure to file figures/violin_Hepatic stellate cells_75-84_down.pdf

Analyzing cell type: Mast cells
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Mast cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_Mast cells_55-64_down.pdf
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Mast cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_Mast cells_65-74_down.pdf
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Mast cells_75-84_up.pdf
WARNING: saving figure to file figures/violin_Mast cells_75-84_down.pdf

Analyzing cell type: Clara cells
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Clara cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_Clara cells_55-64_down.pdf
  Skipping 65-74 (Young Elderly) - insufficient cells (ref: 24, target: 10)
  Skipping 75-84 (Senior) - insufficient cells (ref: 24, target: 2)

Analyzing cell type: B cells
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_B cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_B cells_65-74_down.pdf
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_B cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_B cells_55-64_down.pdf
  Comparing 75-84 (Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_B cells_75-84_up.pdf
WARNING: saving figure to file figures/violin_B cells_75-84_down.pdf
  Skipping <40 (Young) - insufficient cells (ref: 257, target: 1)

Analyzing cell type: Goblet cells
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Goblet cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_Goblet cells_55-64_down.pdf
  Skipping 65-74 (Young Elderly) - insufficient cells (ref: 65, target: 1)

Analyzing cell type: Peritubular myoid cells
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Peritubular myoid cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_Peritubular myoid cells_65-74_down.pdf

Analyzing cell type: Ependymal cells
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Ependymal cells_65-74_up.pdf
WARNING: saving figure to file figures/violin_Ependymal cells_65-74_down.pdf
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Ependymal cells_55-64_up.pdf
WARNING: saving figure to file figures/violin_Ependymal cells_55-64_down.pdf
  Skipping 75-84 (Senior) - insufficient cells (ref: 157, target: 10)
  Skipping <40 (Young) - insufficient cells (ref: 157, target: 4)

Analyzing cell type: Oligodendrocytes
  Comparing 65-74 (Young Elderly) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Oligodendrocytes_65-74_up.pdf
WARNING: saving figure to file figures/violin_Oligodendrocytes_65-74_down.pdf
  Comparing 55-64 (Pre-Senior) vs 40-54 (Mid-Adult)
WARNING: saving figure to file figures/violin_Oligodendrocytes_55-64_up.pdf
WARNING: saving figure to file figures/violin_Oligodendrocytes_55-64_down.pdf
  Skipping <40 (Young) - insufficient cells (ref: 110, target: 5)

Analyzing cell type: Enteroendocrine cells
  Skipping <40 (Young) - insufficient cells (ref: 226, target: 1)
/home/mukul.ranjan/Documents/cancer/SCELA/age_specific_degs.py:162: FutureWarning: Downcasting object dtype arrays on .fillna, .ffill, .bfill is deprecated and will change in a future version. Call result.infer_objects(copy=False) instead. To opt-in to the future behavior, set `pd.set_option('future.no_silent_downcasting', True)`
  summary_df = summary_df.fillna(0).astype(int)

Analysis complete!

Summary of differential expression analysis:
                                  55-64 (Pre-Senior)  ...  <40 (Young)
T cytotoxic cells                                112  ...           21
Gamma delta T cells                              181  ...            0
T memory cells                                   615  ...            0
T cells                                           46  ...           51
Plasma cells                                     259  ...            0
Macrophages                                      457  ...          956
B cells naive                                    202  ...            0
Endothelial cells                                581  ...         1105
Myoepithelial cells                              721  ...            0
Neutrophils                                      306  ...            0
Mammary epithelial cells                        5345  ...            0
NK cells                                          89  ...            0
Pulmonary alveolar type II cells                2827  ...            0
Dendritic cells                                  165  ...           79
Plasmacytoid dendritic cells                       6  ...            0
Hepatic stellate cells                           567  ...         1127
Mast cells                                       468  ...            0
Clara cells                                      337  ...            0
B cells                                           84  ...            0
Goblet cells                                    3430  ...            0
Peritubular myoid cells                            0  ...            0
Ependymal cells                                  133  ...            0
Oligodendrocytes                                 123  ...            0

[23 rows x 4 columns]
