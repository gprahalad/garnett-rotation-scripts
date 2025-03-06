import pandas as pd
import scanpy as sc
import numpy as np

degs = pd.read_csv('/lustre/scratch125/casm/team215mg/pg21_rotation/pbulk_results/G12D_DEG_list.csv')
reference_adata = sc.read_h5ad('/lustre/scratch125/casm/team215mg/pg21_rotation/Coefficients_CellDrift/selection_rep_1/CellDrift_object_time_24.h5ad')


# define list of cell states, treatments and genes
cell_states = np.unique(reference_adata.obs['final_group_labeled'])
treatments = np.unique(reference_adata.obs['drug'])
treatments = np.delete(treatments, [0])
gene_list = degs['gene']

for state in cell_states:
    for drug in treatments:
        genes_and_clusters = pd.read_csv(f'/lustre/scratch125/casm/team215mg/pg21_rotation/Temporal_CellDrift/{state}/{drug}/FDA_clusters.txt', sep = '\t')
        genes_and_clusters.loc[gene_list.isin(genes_and_clusters), ['genes', 'clusters_fuzzy']]