# --------
# Details:  Prahalad Giridhar (pg21@sanger.ac.uk) - 9 January 2025
# Function: Create a semi-random subset of adata for one organoid for testing models. This script only includes HVGs, and then force-stratifies by drug treatment and timepoint,
#           within which 3000 random cells are selected.
# NB:       This is meant to be run on a CellDrift environment which functions on Python 3.7, and as such Scanpy version <3.11. If upgrading you should change the subsampling function
#           to the non-deprecated function.
# --------

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

# Import adata object (one organoid)
adata = sc.read("/lustre/scratch125/casm/team215mg/pg21_rotation/final_sample_filtered_parse_pp_COLO_005.h5ad")

# Import annotation file to annotate adata.obs
annotated_meta_path = '/lustre/scratch125/casm/team215mg/pg21_rotation/final_meta_parse_processed_annotated_states.csv'
annotated_meta_df = pd.read_csv(annotated_meta_path, sep=',', index_col=0)
print('Imported cell states succesfully')

# Add annotated cell states to adata.obs 
adata.obs = adata.obs.merge(how='left', right=annotated_meta_df[['module_sample', 'sample_grouped_mod', 'grouped_mod', 'final_group_labeled']], left_index=True, right_index=True)
print('Merged cell states into adata.obs succesfully')

# Define drug treatments for subsetting whole data by treatment
treatments = ["DMSO", "MRTX1133", "MRTX1133_afatinib", "MRTX1133_SHP099", "afatinib", "SHP099", "none"]

# Define timepoints for subsetting whole adata by timepoint
timepoints = ["24hr", "48hr", "72hr"]

# Define cell states for stratifying subsetted data by
cell_states = np.unique(adata.obs['final_group_labeled'])
print(f'Cell states identified: {cell_states}')

# Subset the whole adata object such that adata.var only includes HVGs
adata = adata[:, adata.var['highly_variable'] == True]
print('Filtered adata.var by HVGs succesfully')
print(f'Number of cells pre-filtering = {adata.shape[0]}')

# Using the adata object that only includes HVGs, loop through adata by drug treatment and timepoints per treatment. For each subset, pull out 3000 cells if possible, then append together.
adatas_cells = []
num_cells = 1000
for drug in treatments:
    adata_drug = adata[adata.obs['drug'] == drug, :].copy()
    print(f'Processing treatment {drug}')
    for time in timepoints:
        adata_time = adata_drug[adata_drug.obs['time'] == time, :].copy()
        print(f'Processing timepoint {time}')
        for state in cell_states:
            adata_state = adata_time[adata_time.obs['final_group_labeled'] == state]
            subset_cells = adata_state.shape[0]
            if subset_cells > 0:
                if subset_cells > num_cells:
                    subsample_adata = sc.pp.subsample(adata_state, n_obs = num_cells, copy = True)
                    adatas_cells.append(subsample_adata)
                    print(f'Succesfully subsampled treatment {drug} at {time} in {state} with {num_cells} cells')
                else: 
                    subsample_adata = sc.pp.subsample(adata_state, n_obs = subset_cells, copy = True)
                    adatas_cells.append(subsample_adata)
                    print(f'Succesfully subsampled treatment {drug} at {time} in {state} with {subset_cells} cells')
            else:
                print(f'Skipped treatment {drug} at {time} in {state} due to cells = 0')

# Write out subsetted adata
final_adata_subset = ad.concat(adatas_cells) 
final_adata_subset.write('new_per_cellstate_subset_COLO_005_130125.h5ad')