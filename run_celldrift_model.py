# --------
# Details:  Prahalad Giridhar (pg21@sanger.ac.uk) - 15 January 2025
# Function: Format and run the CellDrift model.
# NB:       This script must be run on a CellDrift environment which functions on Python 3.7. The _model.py script of CellDrift has been edited, the version in
#           /lustre/scratch125/casm/team215mg/pg21_rotation/software/miniforge3/envs/celldrift_git/lib/python3.7/site-packages/CellDrift/ should be used.
# --------

import scanpy as sc
import pandas as pd
import numpy as np
import CellDrift as ct
import random
import anndata as ad

adata = sc.read('/lustre/scratch125/casm/team215mg/pg21_rotation/new_per_cellstate_subset_COLO_005_130225.h5ad')

adata.obs['time'] = adata.obs['time'].str[:2]
adata.obs['time'] = adata.obs['time'].astype(int)
adata.obs['final_group_labeled'] = adata.obs['final_group_labeled'].replace({'Stress, EMT': 'Stress'})
adata.obs['final_group_labeled'] = adata.obs['final_group_labeled'].replace({'Hybrid/Fetal-like': 'Fetal'})
adata.obs['final_group_labeled'] = adata.obs['final_group_labeled'].replace({'Hybrid/Stem-like': 'Stem'})
adata.obs['final_group_labeled'] = adata.obs['final_group_labeled'].replace({'Secretory-like': 'Secretory'})

adata = ct.setup_celldrift(
    adata,
    cell_type_key = 'final_group_labeled',
    perturb_key = 'drug',
    time_key = 'time',
    control_name = 'DMSO',
    perturb_name = None,
    size_factor_key = 'total_counts',
    batch_key = None,
    n_reps = 10,
    n_cells_perBlock = 100,
    use_pseudotime = False,
    min_cells_perGene = 150
)

# check the counts of each cell type
print(adata.obs['final_group_labeled'].value_counts())

adata = ct.model_timescale(
    adata,
    n_processes = 16,
    chunksize = 100,
    pairwise_contrast_only = False,
    adjust_batch = False
)