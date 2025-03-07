# garnett-rotation-scripts
Scripts for my rotation in the Garnett group during which I created models to predict cell state distribution of organoids after drug treatment.

### Order of operations
1. Subset data from single organoid (draft_split_adata.py)
2. Run CellDrift model (run_celldrift_model.py) and if you want trajectories use the temporal functions (run_temporal_celldrift.py)
3. Use either random_forest.ipynb or more_cells_xgb.ipynb for modelling

The rest of the scripts in this directory are either auxiliary functions based off what I have in the notebooks (such as make_pivot.py) or attempts with feature selection or optimization I did not think were very impactful in the end (Hyperopt or SHAP score evaluation). I have still included the latter in case they are of interest in the future and to document some of the things I tried that were not as successful.
