# --------
# Details:  Prahalad Giridhar (pg21@sanger.ac.uk) taken from code snippet from Szen Toh (st25@sanger.ac.uk), 9 Jan 2025
# Function: Import cell states and subset by taking small number of genes and cells for understanding data structures and playing with CellDrift. This step is included in other 
#           adata splitting scripts, so this script should be ignored if not being used for the purposes of very restricted testing.
# --------

import pandas as pd
import anndata as ad
import scanpy as sc

# Import adata object (one organoid)
adata = sc.read("/lustre/scratch125/casm/team215mg/pg21_rotation/final_sample_filtered_parse_pp_COLO_005.h5ad")

# Import annotation file to annotate adata.obs
annotated_meta_path = '/lustre/scratch125/casm/team215mg/pg21_rotation/final_meta_parse_processed_annotated_states.csv'
annotated_meta_df = pd.read_csv(annotated_meta_path, sep=',', index_col=0)
print('Imported cell states succesfully')

# Add annotated cell states to adata.obs 
adata.obs = adata.obs.merge(how='left', right=annotated_meta_df[['module_sample', 'sample_grouped_mod', 'grouped_mod', 'final_group_labeled']], left_index=True, right_index=True)
print('Merged cell states into adata.obs succesfully')

# Subset adata in a rudimentary way for data handling testing
adata_subset = adata[0:1000, 0:300]
adata_subset.write_h5ad('test_subset_adata.h5ad')