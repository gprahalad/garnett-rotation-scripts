import scanpy as sc
import anndata as ad
import pandas as pd

# Import adata object (one organoid)
adata = sc.read("/lustre/scratch125/casm/team215mg/pg21_rotation/final_sample_filtered_parse_pp_COLO_005.h5ad")

# Import annotation file to annotate adata.obs
annotated_meta_path = '/lustre/scratch125/casm/team215mg/pg21_rotation/final_meta_parse_processed_annotated_states.csv'
annotated_meta_df = pd.read_csv(annotated_meta_path, sep=',', index_col=0)
print('Imported cell states succesfully')

# Add annotated cell states to adata.obs 
adata.obs = adata.obs.merge(how='left', right=annotated_meta_df[['module_sample', 'sample_grouped_mod', 'grouped_mod', 'final_group_labeled']], left_index=True, right_index=True)
print('Merged cell states into adata.obs succesfully')

# Subset the whole adata object such that adata.var only includes HVGs
adata = adata[:, adata.var['highly_variable'] == True]
print('Filtered adata.var by HVGs succesfully')
print(f'Number of cells pre-filtering = {adata.shape[0]}')

# Randomly subset by a fraction of the total
adata_subset = sc.pp.subsample(adata, fraction=0.2)
adata_subset.write_h5ad('random_fraction_subset_COLO_005.h5ad')
print(f'Number of cells post-filtering = {adata_subset.shape[0]}')