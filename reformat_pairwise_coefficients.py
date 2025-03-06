import numpy as np
import pandas as pd

df_pairwiseComp_all = pd.read_csv('/lustre/scratch125/casm/team215mg/pg21_rotation/scripts/Coefficients_CellDrift/genes_771_ncellblock_10/Contrast_coefficients_time_25.txt',
                                   sep = '\t')

for gene in df_pairwiseComp_all
