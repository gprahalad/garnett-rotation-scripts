import pandas as pd
import numpy as np

coefficients = pd.DataFrame()

# define timepoints as necessary
timepoints = ['24', '48', '72']

# loop through all reps to include and timepoints, modify contrast value so it is unique per rep, concatenate
for rep in range(1,11):
    for time in timepoints:
        temp_coeff = pd.read_csv(f'/lustre/scratch125/casm/team215mg/pg21_rotation/Coefficients_CellDrift/selection_rep_{rep}/Contrast_Coefficients_time_{time}.txt', sep = '\t')
        temp_coeff['contrast'] = temp_coeff['contrast'] + f'_{rep}'
        temp_coeff['time'] = time
        coefficients = pd.concat([coefficients, temp_coeff], ignore_index = True)

# fix cell type labels for model and simplify perturbation column
coefficients['cell_type'] = coefficients['cell_type'].str.split(',', n = 1).str.get(-1).str.split("'").str.get(-2)
coefficients['perturbation'] = coefficients['perturbation'].str.split(',', n = 1).str.get(0).str.split("'").str.get(-2)
coefficients['cell_type'] = coefficients['cell_type'].map({'Fetal':0, 'Stem':1, 'Cycling':2, 'Secretory':3, 'Stress':4, 'Respiration':5})

# filter to desired conditions, ie. drug combinations
coefficients = coefficients[(coefficients['perturbation'] == 'MRTX1133_SHP099') | (coefficients['perturbation'] == 'MRTX1133_afatinib')]

# pivot to have genes as columns instead of per row, row is now contrast per rep
pivot = coefficients.pivot_table('z', ['contrast', 'cell_type', 'time'], 'gene').fillna(0)
pivot.reset_index(drop = False, inplace = True)

# calc deltas per gene per rep, mega inefficient unpythonic stackoverflow heart attack inducer
for contrast in np.unique(coefficients['contrast']):
    coeff_filtered = pivot[pivot['contrast'] == contrast]
    for gene in np.unique(coefficients['gene']):
        z_72 = coeff_filtered.loc[coeff_filtered['time'] == '72', [gene]]
        z_48 = coeff_filtered.loc[coeff_filtered['time'] == '48', [gene]]
        z_24 = coeff_filtered.loc[coeff_filtered['time'] == '24', [gene]]
        pivot.loc[((pivot['contrast'] == contrast) & (pivot['time'] == '24')), [f'{gene}_delta']] = z_24[gene].values[0] - z_48[gene].values[0]
        pivot.loc[((pivot['contrast'] == contrast) & (pivot['time'] == '48')), [f'{gene}_delta']] = z_48[gene].values[0] - z_24[gene].values[0]
        pivot.loc[((pivot['contrast'] == contrast) & (pivot['time'] == '72')), [f'{gene}_delta']] = z_72[gene].values[0] - z_48[gene].values[0]

pivot.to_csv('/lustre/scratch125/casm/team215mg/pg21_rotation/new_all_pivot.tsv', sep = '\t')
