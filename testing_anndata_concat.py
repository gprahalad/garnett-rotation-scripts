import anndata as ad
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix
import random
import string
import os
import argparse

def parse_args():
    """Parses command-line options for main()."""
    summary = 'Splitting and re-concatenation of AnnData files for testing purposes'
    parser = argparse.ArgumentParser(description=summary)
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-working_dir',
                               required = True,
                               type = str,
                               help = 'Path to working directory for retrieving test split adatas.')
    try:    
        opts = parser.parse_args()
    except:
        sys.exit(0)
    
    return opts

args = parse_args()

counts = csr_matrix(np.random.poisson(1, size=(100, 100)), dtype=np.float32)
adata = ad.AnnData(counts)

region_name_list = string.ascii_lowercase[:10]
adata.obs['regionName'] = random.choices(region_name_list, k = 100)

for region_name in np.unique(adata.obs['regionName']):
    adata_subset = adata[adata.obs['regionName']==region_name, :]
    adata_subset.write_h5ad(f'split_{region_name}_region.h5ad')

split_files = os.listdir(args.working_dir)

adatas = []
for split_adatas in split_files:
    small_adata = sc.read_h5ad(split_adatas)
    adatas.append(small_adata)
    print(split_adatas, small_adata.shape[0], small_adata.obsm[0])

adata_concat = ad.concat(adatas)
adata_concat.obs_names_make_unique()
adata_concat.write_h5ad(os.path.join(args.working_dir, f'recombined_adatas.h5ad'))

print("concat", adata_concat.shape[0])