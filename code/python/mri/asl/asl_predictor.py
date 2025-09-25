import pandas as pd
import numpy as np

# Paths and filenames
data_partition_path = '/data/'
study_data_path = f'{data_partition_path}/UNSAM/CovidProject2/'
results_path = f'{study_data_path}/DataAnalysis/'
response_excel_filename = f'{study_data_path}/Respuestas.xlsx'
summary_excel_filename = f'{study_data_path}/ResumenRespuestas.xlsx'
volunteers_excel_filename = f'{study_data_path}/VoluntariosProyectoCovidProlongado.xlsx'

# Load ASL CSVs
asl_csv_paths = '/home/martin/data/UNSAM/CovidProject2/DataAnalysis/ASL/CleanData/'
scov_gm_filename = f'{asl_csv_paths}/scov_total_gm_clean.csv'

scov_table = pd.read_csv(scov_gm_filename)

# Match groups
summary_table = pd.read_excel(summary_excel_filename, sheet_name='resumenTotal')
indices_scov = []
indices_non_cov = []

for i, pid in enumerate(scov_table['participant_id']):
    ind = summary_table.index[summary_table['ID'] == pid].tolist()
    if ind:
        indices_scov.append(ind[0])
    else:
        print(f'Warning: Subject {pid} not found')
        indices_non_cov.append(i)

summary_table_matched = summary_table.loc[indices_scov].reset_index(drop=True)
group = summary_table_matched['Grupo']
scov_table['group'] = group
scov_table['sex'] = summary_table_matched['Genero']
scov_table['age'] = summary_table_matched['Edad']

# Simple oversampling
indices_control = scov_table['group'] == 'CONTROL'
mat_to_concate = pd.concat([scov_table[indices_control]] * 3, ignore_index=True)
scov_table = pd.concat([scov_table, mat_to_concate], ignore_index=True)
indices_shuffle = np.random.permutation(scov_table.shape[0])
scov_table_shuffled = scov_table.iloc[indices_shuffle].reset_index(drop=True)
group_shuffled = scov_table_shuffled['group']