# --------
# Details:  Prahalad Giridhar (pg21@sanger.ac.uk), code provided by Szen Toh (st25@sanger.ac.uk) - 10 January 2025
# Function: Pre-processing and normalization to use on AnnData files after subsetting - in this case used for normalizing after subsetting to a 54000 cells from a single organoid.
# --------

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import csv
import gzip
import os
import scipy.io
import seaborn as sns
import math
from sklearn import preprocessing
import matplotlib.pyplot as plt
import random
import sys
from scipy import sparse
from sklearn.metrics import adjusted_rand_score


plt.rcParams["axes.grid"] = False

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

random.seed(1)

#sample_name = sys.argv[1]

# insert path here as required
# Load, normalize and log-transform count data
adata = sc.read_h5ad('/lustre/scratch125/casm/team215mg/pg21_rotation/random_subset_COLO_005.h5ad')
print("Read successful!")

# use raw counts 
adata.X = adata.layers["counts"].copy()

print('Done copy counts!')

del adata.layers['norm']
del adata.layers['log_norm']
print(len(adata.obs.index))

cell_type_genes = ['ANPEP','PPARG','AGR2','GCG','CHGA','SCGN','S100B','SLC26A3','IAPP', 'TFF3', 'SST', 'KRT20', 'GFRA3','FABP1','CA1','FABP2','VIL1', 
'CDH17','KLF5','MUC13','PHGR1','SLC51B', 'ACSL5','AQP8','SLC39A5','MAOA','DGAT1','CHGB','CPE','MUC2','ITLN1','SPINK4']

# include genes for adsorptive - 5, secretory - 4, intestine - 4, stemness - 6, dev - 5, EMT - 5, neuroendo - 5, squamous - 5
impt_genes = ['ALDOB', 'FABP2', 'ADH4', 'DMBT1', 'CCL25', 'TFF1', 'MUC5AC', 'CLCA1', 'CDX1', 'URAD', 'GDF15', 'NOX1', 'LGR5', 'ASCL2', 'IGF2', 'SMOC2', 'RGMB', 'SOX9', 'NKD1', 'BMP4', 
              'NR2F1', 'PROX1', 'WNT5B', 'CDH2', 'VIM', 'CLU', 'FN1', 'CDH1', 'NEUROD1', 'INSM1', 'RFX6', 'SCG3', 'KRT31', 'KRT31', 'FGF3', 'HEPHL1', 'PRR9', 'KRT5']

# Perform PCA on highly variable genes

# normalize 
sc.pp.normalize_total(adata, target_sum=1e4)
adata.layers["norm"] = adata.X.copy()

sc.pp.log1p(adata)
adata.layers["log_norm"] = adata.X.copy()

# Required to identify HVGs but this pre-processing script is being used for normalizing data that already has labelled HVGs

'''
# find HVGs 
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# select relevant genes as highly variable 
adata.var.loc[cell_type_genes, 'highly_variable'] = True
adata.var.loc[impt_genes, 'highly_variable'] = True
'''

# calculate PCA and find cumulative variance explained
def pca_cml_var(adata_sample):
    adata_pca = sc.tl.pca(adata_sample, n_comps = 800, copy = True)
    
    # Cumulative variance explained:
    cml_var_explained = np.cumsum(adata_pca.uns['pca']['variance_ratio'])
    x = range(len(adata_pca.uns['pca']['variance_ratio']))
    y = cml_var_explained
    plt.scatter(x, y, s=4)
    plt.xlabel('PC')
    plt.ylabel('Cumulative variance explained')
    plt.title('Cumulative variance explained by PCs')
    
    # set the minimum cumulative fraction of variance explained
    min_cml_frac = 0.7
    # find the number of PCs that together explain the minimal cum. fraction of variance
    n_pcs = next(idx for idx,cml_frac in enumerate(cml_var_explained) if cml_frac > min_cml_frac)
    print('Number of PCs that together explain a fraction of ' + str(min_cml_frac) + ' of the variance: ' + str(n_pcs))

    return adata_pca, n_pcs

adata, n_pcs = pca_cml_var(adata)

# subset to PCs that explain 70% of variation
sc.tl.pca(adata, n_comps=n_pcs)

# run neighbors and UMAP
sc.pp.neighbors(adata, n_neighbors=15)

sc.tl.umap(adata)