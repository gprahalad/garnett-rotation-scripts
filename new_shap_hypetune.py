import pandas as pd
import numpy as np
from shaphypetune import BoostSearch, BoostRFE, BoostRFA, BoostBoruta
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from xgboost import XGBClassifier
from scipy import stats

coefficients = pd.DataFrame()
timepoints = ['24', '48', '72']

for rep in range(1,11):
    for time in timepoints:
        temp_coeff = pd.read_csv(f'/lustre/scratch125/casm/team215mg/pg21_rotation/Coefficients_CellDrift/original/selection_rep_{rep}/Contrast_Coefficients_time_{time}.txt', sep = '\t')
        temp_coeff['contrast'] = temp_coeff['contrast'] + f'_{rep}'
        temp_coeff['time'] = time
        coefficients = pd.concat([coefficients, temp_coeff], ignore_index = True)

coefficients['cell_type'] = coefficients['cell_type'].str.split(',', n = 1).str.get(-1).str.split("'").str.get(-2)
coefficients['perturbation'] = coefficients['perturbation'].str.split(',', n = 1).str.get(0).str.split("'").str.get(-2)
coefficients['cell_type'] = coefficients['cell_type'].map({'Fetal':0, 'Stem':1, 'Cycling':2, 'Secretory':3, 'Stress':4, 'Respiration':5})
coefficients = coefficients[(coefficients['perturbation'] == 'MRTX1133_SHP099') | (coefficients['perturbation'] == 'MRTX1133_afatinib')]

pivot = pd.read_csv('/lustre/scratch125/casm/team215mg/pg21_rotation/all_pivot.tsv', sep = '\t')

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

xgb_tuned = XGBClassifier(learning_rate = 0.2,
                          max_depth = 18,
                          n_estimators = 421,
                          subsample = 0.2,
                          early_stopping_rounds = 20)

xgb_tuned.fit(pivot_24, y_24, eval_set = [(pivot_72, y_72)])

param_dist = {
    'learning_rate': stats.uniform(0.01, 0.2),
    'max_depth': [10, 15],
    'subsample': stats.uniform(0, 1),
    'reg_alpha': stats.uniform(0, 5),
    'reg_lambda': stats.uniform(0, 5)
}

model = BoostRFE(xgb_tuned, param_grid = param_dist, n_iter = 200,
                 min_features_to_select = 1, step = 1, n_jobs = 5000,
                 sampling_seed = 0)
model.fit(pivot_24, y_24, eval_set = [(pivot_72, y_72)], verbose = 0)

predictions = model.predict(pivot_72)
accuracy = accuracy_score(y_72, predictions)

print(f'Accuracy of tuned model = {accuracy * 100}')
print(model.best_params_)
print(model.best_score_)
print(pivot_24.columns[model.support_])