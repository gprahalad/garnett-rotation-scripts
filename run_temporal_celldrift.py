# --------
# Details:  Prahalad Giridhar (pg21@sanger.ac.uk) - 16 January 2025
# Function: Import pre-generated CellDrift coefficients and use these to plot temporal gene expression patterns using 3 smoothing methods. 
# NB:       This is meant to be run on a CellDrift environment which functions on Python 3.7. The _temporal.py script of CellDrift has been edited, the version in 
#           /lustre/scratch125/casm/team215mg/pg21_rotation/software/miniforge3/envs/celldrift_git/lib/python3.7/site-packages/CellDrift/ should be used.
# --------

import scanpy as sc
import pandas as pd
import numpy as np
import CellDrift as ct
import random
import anndata as ad
import os

# load data
df_zscore = pd.read_csv('/lustre/scratch125/casm/team215mg/pg21_rotation/Temporal_CellDrift/Contrast_Coefficients_combined_zscores_.txt', sep = '\t', header = 0, index_col = 0)
df_meta = pd.read_csv('/lustre/scratch125/casm/team215mg/pg21_rotation/Temporal_CellDrift/Contrast_Coefficients_combined_metadata_.txt', sep = '\t', header = 0, index_col = 0)
reference_adata = sc.read_h5ad('/lustre/scratch125/casm/team215mg/pg21_rotation/Coefficients_CellDrift/selection_rep_1/CellDrift_object_time_24.h5ad')

# redefine timepoints
time_origin = [24, 48, 72]
time_new = [1, 2, 3]
time_dict = dict(zip(time_origin, time_new))
df_meta['time'] = [time_dict[i] for i in df_meta['time']]

# run FDA across cell types and perturbations
fda = ct.FDA(df_zscore, df_meta)

# define list of cell states and treatments
cell_states = np.unique(reference_adata.obs['final_group_labeled'])
treatments = np.unique(reference_adata.obs['drug'])
treatments = np.delete(treatments, [0])

print(cell_states)

# check if output directory would exist, if not create, then output clusters and figures per state and perturbation
for state in cell_states:
    if not os.path.exists(f'/lustre/scratch125/casm/team215mg/pg21_rotation/Temporal_CellDrift/{state}'):
        os.makedirs(f'/lustre/scratch125/casm/team215mg/pg21_rotation/Temporal_CellDrift/{state}')
    for drug in treatments:
        print(state, drug)
        if not os.path.exists(f'/lustre/scratch125/casm/team215mg/pg21_rotation/Temporal_CellDrift/{state}/{drug}'):
            os.makedirs(f'/lustre/scratch125/casm/team215mg/pg21_rotation/Temporal_CellDrift/{state}/{drug}')
        fd, genes = fda.create_fd_genes(genes = df_zscore.index.values, cell_type = state, perturbation = drug)
        df_cluster = ct.fda_cluster(fd, genes, seed = 42, n_clusters = 8,
                                    output_folder = f'/lustre/scratch125/casm/team215mg/pg21_rotation/Temporal_CellDrift/{state}/{drug}/')
        ct.draw_smoothing_clusters(fd, df_cluster, n_neighbors = 2, bandwidth = 1, cluster_key = 'clusters_fuzzy', 
                                    output_folder = f'/lustre/scratch125/casm/team215mg/pg21_rotation/Temporal_CellDrift/{state}/{drug}/')
        #ct.run_anova_test(fd, genes, state, treatments)
