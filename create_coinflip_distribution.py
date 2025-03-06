# --------
# Details:  Prahalad Giridhar (pg21@sanger.ac.uk) - 24 February 2025
# Function: Identify collinearity between genes and other DE genes, then coinflip in a randomly distributed order to evaluate resulting model performance.
#           The purpose of doing so is to see if collinear genes are interchangeable, as doing so would allow for informative biomarkers to be selected.
# --------

import numpy as np
import random
import pandas as pd
import xgboost as xgb
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

coefficients = pd.DataFrame()
timepoints = ['24', '48', '72']

for rep in range(1,11):
    for time in timepoints:
        temp_coeff = pd.read_csv(f'/lustre/scratch125/casm/team215mg/pg21_rotation/Coefficients_CellDrift/selection_rep_{rep}/Contrast_Coefficients_time_{time}.txt', sep = '\t')
        temp_coeff['contrast'] = temp_coeff['contrast'] + f'_{rep}'
        temp_coeff['time'] = time
        coefficients = pd.concat([coefficients, temp_coeff], ignore_index = True)

coefficients['cell_type'] = coefficients['cell_type'].str.split(',', n = 1).str.get(-1).str.split("'").str.get(-2)
coefficients['perturbation'] = coefficients['perturbation'].str.split(',', n = 1).str.get(0).str.split("'").str.get(-2)
coefficients['cell_type'] = coefficients['cell_type'].map({'Fetal':0, 'Stem':1, 'Cycling':2, 'Secretory':3, 'Stress':4, 'Respiration':5})

coefficients = coefficients[(coefficients['perturbation'] == 'MRTX1133_SHP099') | (coefficients['perturbation'] == 'MRTX1133_afatinib')]

accuracies = []

pivot = pd.read_csv('/lustre/scratch125/casm/team215mg/pg21_rotation/new_all_pivot.tsv', sep = '\t')

pivot_24 = pivot[pivot['time'] == 24]
pivot_48 = pivot[pivot['time'] == 48]
pivot_72 = pivot[pivot['time'] == 72]
y_24 = pivot_24['cell_type'].to_numpy()
y_48 = pivot_48['cell_type'].to_numpy()
y_72 = pivot_72['cell_type'].to_numpy()

for gene in np.unique(coefficients['gene']):
    pivot_24[f'{gene}_delta'] = pivot_48[f'{gene}_delta'].values
    pivot_48[f'{gene}_delta'] = pivot_72[f'{gene}_delta'].values

pivot_24.drop(['Unnamed: 0', 'contrast', 'cell_type', 'time'], axis = 1, inplace = True)
X_24 = pivot_24.to_numpy()
pivot_48.drop(['Unnamed: 0', 'contrast', 'cell_type', 'time'], axis = 1, inplace = True)
X_48 = pivot_48.to_numpy()
pivot_72.drop(['Unnamed: 0', 'contrast', 'cell_type', 'time'], axis = 1, inplace = True)
X_72 = pivot_72.to_numpy()

X_24_train, X_24_test, y_24_train, y_24_test = train_test_split(pivot_24, y_24, test_size = 0.2)
X_48_train, X_48_test, y_48_train, y_48_test = train_test_split(pivot_48, y_48, test_size = 0.2)
X_72_train, X_72_test, y_72_train, y_72_test = train_test_split(pivot_72, y_72, test_size = 0.2)

for threshold in [0.7]:
    for i in range(1, 201):
        corr_matrix = pivot_24.corr(method = 'pearson')
        corr_matrix = corr_matrix.stack()
        corr_matrix = corr_matrix[((corr_matrix > threshold)) & (corr_matrix != 1)]
        corr_matrix = corr_matrix.to_frame()
        corr_gene_list = np.unique(corr_matrix.index.get_level_values(0))
        random.shuffle(corr_gene_list)

        for gene in corr_gene_list:
            for corr in corr_matrix[corr_matrix.index.get_level_values(0) == gene].index.get_level_values(1):
                if (random.randint(0, 1) == 0) & (corr in np.unique(corr_matrix.index.get_level_values(0))) & (gene in np.unique(corr_matrix.index.get_level_values(0))):
                    corr_matrix.drop(labels = corr, axis = 0, inplace = True)
                if (random.randint(0, 1) == 1) & (gene in np.unique(corr_matrix.index.get_level_values(0))):
                    corr_matrix.drop(labels = gene, axis = 0, inplace = True)
        
        pivot_24_dropped = pivot_24.drop(labels = np.unique(corr_matrix.index.get_level_values(1)), axis = 1)
        pivot_48_dropped = pivot_48.drop(labels = np.unique(corr_matrix.index.get_level_values(1)), axis = 1)
        pivot_72_dropped = pivot_72.drop(labels = np.unique(corr_matrix.index.get_level_values(1)), axis = 1)

        xgb_tuned = XGBClassifier(learning_rate = 0.2,
                            max_depth = 16,
                            n_estimators = 460,
                            subsample = 0.3)
        xgb_tuned.fit(pivot_24_dropped, y_24)

        predictions = xgb_tuned.predict(pivot_72_dropped)
        accuracy = accuracy_score(y_72, predictions)

        accuracies.append(accuracy * 100)

fig, axs = plt.subplots(1, sharey = True, tight_layout = True)
axs.hist(accuracies, bins = 10)
plt.savefig('coinflip_distribution.svg')