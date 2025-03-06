import os
import pandas as pd
import re

DGE_result_file_path = '/lustre/scratch125/casm/team215mg/pg21_rotation/pbulk_results/'

all_DGE_results_kras = []

for filename in os.listdir(DGE_result_file_path):
    if filename.startswith('G12') and filename.endswith("_all_genes_annot.csv"):
        if not filename.endswith("org_all_genes_annot.csv"):
            DGE_result_path = os.path.join(DGE_result_file_path, filename)
            print(f'reading file {DGE_result_path}')
        
            DGE_result_annot = pd.read_csv(DGE_result_path, sep=',')
        
            sample_condition = re.sub('(.*)_(.*)_(.*)_(.*)_all_genes_annot.csv','\\1', os.path.basename(filename))
            print(f'Reading {sample_condition}')
        
            DGE_result_annot['sample_condition'] = sample_condition
        
            DGE_result_annot[['mut_status','condition']] = DGE_result_annot['sample_condition'].str.split('_', n=1, expand=True)
            DGE_result_annot['time'] = DGE_result_annot['condition'].str.extract(r'(\d+(\s*(?:hr)))')[0]
            DGE_result_annot['drug'] = DGE_result_annot['condition'].str.extract(r'(.*?)_(\d+)')[0]
            DGE_result_annot = DGE_result_annot[DGE_result_annot['diffexpressed'] != 'NO']
            DGE_result_annot = DGE_result_annot[DGE_result_annot['mut_status'] == 'G12D']

            all_DGE_results_kras.append(DGE_result_annot)

all_DGE_kras = pd.concat(all_DGE_results_kras)
all_DGE_kras.to_csv('G12D_DEG_list.csv')