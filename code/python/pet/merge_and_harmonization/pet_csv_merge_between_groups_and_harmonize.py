import os
import site
import subprocess
import sys
import shlex
import nibabel as nib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Configuration
pet_only = False  # Set to True if you want to preprocess only PET data without MRI
pet_mri = True  # Set to True if you want to preprocess PET data with MRI

# Output path for merged summary csv files
merged_summary_csv_path = './data/processed_pet/merged_summary_csv/'
os.makedirs(merged_summary_csv_path, exist_ok=True)

# Path UNSAM COVID data
pet_covid_unsam_path = '/home/martin/data_imaging/CovidProject/FullStudy/PET//'
preprocessed_path = os.path.join(pet_covid_unsam_path, 'preprocessed')
preprocessed_with_mri_path = os.path.join(pet_covid_unsam_path, 'ProcesadosCOVID_PET_MRI_09_2025')
# Output path for summary csv file
covid_unsam_summary_csv_path = os.path.join(preprocessed_path, 'rois_uptake')
covid_unsam_demographics_filename = './data/study_summary.csv'

# Path UNSAM Controls data
pet_controls_unsam_path = '/home/martin/data_imaging/PET/CEUNIM_BrainControls/'
preprocessed_controls_unsam_path = os.path.join(pet_controls_unsam_path, 'preprocessed')
# Output path for summary csv file
controls_unsam_summary_csv_path = os.path.join(preprocessed_controls_unsam_path, 'rois_uptake')
controls_unsam_demographics_filename = os.path.join(pet_controls_unsam_path, 'lista_estudios_pet_v01_anon.csv')

# Paths CERMEP reference data
cermep_path = '/home/martin/data_imaging/PET/CERMEP/ALL/derivatives/'
preprocessed_path = os.path.join(cermep_path, 'preprocessed')
preprocessed_with_mri_path = os.path.join(cermep_path, 'preprocessed_with_mri')
# Output path for summary csv file
cermep_summary_csv_path = os.path.join(cermep_path, 'rois_uptake')
cermep_demographics_filename = os.path.join(cermep_path, '../participants.tsv')

# Filenames for summary csv files
summary_uptake_csv_filename = 'summary_roi_uptake_values.csv'
summary_suvr_mean_csv_filename = 'summary_suvr_mean_values.csv'
summary_suvr_cerebellum_csv_filename = 'summary_suvr_cerebellum_values.csv'

# Create an empty dataframe to store the summary
df_summary = pd.DataFrame()
df_summary_suvr_mean = pd.DataFrame()
df_summary_suvr_cerebellum = pd.DataFrame()

# Read the summary csv files for each groups and merge them into a single dataframe

################# Read the summary csv files for UNSAM COVID data ##########################
#################################################################################################
df_summary_unsam = pd.read_csv(os.path.join(covid_unsam_summary_csv_path, summary_uptake_csv_filename))
df_summary_suvr_mean_unsam = pd.read_csv(os.path.join(covid_unsam_summary_csv_path, summary_suvr_mean_csv_filename))
df_summary_suvr_cerebellum_unsam = pd.read_csv(os.path.join(covid_unsam_summary_csv_path, summary_suvr_cerebellum_csv_filename))
df_demographics_unsam = pd.read_csv(covid_unsam_demographics_filename)
# Normalize ID formatting before filtering/merging to avoid dtype and whitespace mismatches.
df_summary_unsam['ID'] = df_summary_unsam['ID'].astype(str).str.strip()
df_summary_suvr_mean_unsam['ID'] = df_summary_suvr_mean_unsam['ID'].astype(str).str.strip()
df_summary_suvr_cerebellum_unsam['ID'] = df_summary_suvr_cerebellum_unsam['ID'].astype(str).str.strip()
df_demographics_unsam['ID'] = df_demographics_unsam['ID'].astype(str).str.strip()
unsam_ids = set(df_summary_unsam['ID'])
df_demographics_unsam = df_demographics_unsam[df_demographics_unsam['ID'].isin(unsam_ids)]
# Get fields: Edad,Genero,Altura ,Peso
df_demographics_unsam = df_demographics_unsam[['ID', 'Edad', 'Genero', 'Peso']]
# Rename to age', 'sex', 'weight
df_demographics_unsam.rename(columns={'Edad': 'age', 'Genero': 'sex', 'Peso': 'weight'}, inplace=True)
df_demographics_unsam['age'] = pd.to_numeric(df_demographics_unsam['age'], errors='coerce')
df_demographics_unsam['weight'] = pd.to_numeric(df_demographics_unsam['weight'], errors='coerce')
df_demographics_unsam['sex'] = df_demographics_unsam['sex'].astype(str).str.strip()
# Leave only the IDS with PET data, i.e., those that are in the summary csv files
# Add these fields to the summary dataframes after the ID column
df_summary_unsam = pd.merge(df_summary_unsam, df_demographics_unsam, on='ID', how='left')
df_summary_suvr_mean_unsam = pd.merge(df_summary_suvr_mean_unsam, df_demographics_unsam, on='ID', how='left')
df_summary_suvr_cerebellum_unsam = pd.merge(df_summary_suvr_cerebellum_unsam, df_demographics_unsam, on='ID', how='left')



################# Read the summary csv files for UNSAM Controls data ##########################
#################################################################################################
df_summary_controls_unsam = pd.read_csv(os.path.join(controls_unsam_summary_csv_path, summary_uptake_csv_filename))
df_summary_suvr_mean_controls_unsam = pd.read_csv(os.path.join(controls_unsam_summary_csv_path, summary_suvr_mean_csv_filename))
df_summary_suvr_cerebellum_controls_unsam = pd.read_csv(os.path.join(controls_unsam_summary_csv_path, summary_suvr_cerebellum_csv_filename))
df_demographics_controls_unsam = pd.read_csv(controls_unsam_demographics_filename)
# Normalize ID formatting before filtering/merging to avoid dtype and whitespace mismatches.
df_summary_controls_unsam['ID'] = df_summary_controls_unsam['ID'].astype(str).str.strip()
df_summary_suvr_mean_controls_unsam['ID'] = df_summary_suvr_mean_controls_unsam['ID'].astype(str).str.strip()
df_summary_suvr_cerebellum_controls_unsam['ID'] = df_summary_suvr_cerebellum_controls_unsam['ID'].astype(str).str.strip()
df_demographics_controls_unsam['anonimizado'] = df_demographics_controls_unsam['anonimizado'].astype(str).str.strip()
unsam_ids = set(df_summary_controls_unsam['ID'])
df_demographics_controls_unsam = df_demographics_controls_unsam[df_demographics_controls_unsam['anonimizado'].isin(unsam_ids)]
# Get fields: Edad,Genero,Altura ,Peso
df_demographics_controls_unsam = df_demographics_controls_unsam[['anonimizado', 'edad', 'sexo', 'peso']]
# Rename anonimizado to ID
#df_demographics_controls_unsam.rename(columns={'anonimizado': 'ID'}, inplace=True)
# Rename to age', 'sex', 'weight
df_demographics_controls_unsam.rename(columns={'edad': 'age', 'sexo': 'sex', 'peso': 'weight'}, inplace=True)
df_demographics_controls_unsam['age'] = pd.to_numeric(df_demographics_controls_unsam['age'], errors='coerce')
df_demographics_controls_unsam['weight'] = pd.to_numeric(df_demographics_controls_unsam['weight'], errors='coerce')
df_demographics_controls_unsam['sex'] = df_demographics_controls_unsam['sex'].astype(str).str.strip()

# Leave only the IDS with PET data, i.e., those that are in the summary csv files
# Add these fields to the summary dataframes after the ID column
df_summary_controls_unsam = pd.merge(df_summary_controls_unsam, df_demographics_controls_unsam, left_on='ID', right_on='anonimizado', how='left')
df_summary_suvr_mean_controls_unsam = pd.merge(df_summary_suvr_mean_controls_unsam, df_demographics_controls_unsam, left_on='ID', right_on='anonimizado', how='left')
df_summary_suvr_cerebellum_controls_unsam = pd.merge(df_summary_suvr_cerebellum_controls_unsam, df_demographics_controls_unsam, left_on='ID', right_on='anonimizado', how='left')
# Remove anonimizado column after merge
df_summary_controls_unsam.drop(columns=['anonimizado'], inplace=True)
df_summary_suvr_mean_controls_unsam.drop(columns=['anonimizado'], inplace=True)
df_summary_suvr_cerebellum_controls_unsam.drop(columns=['anonimizado'], inplace=True)



################# Read the summary csv files for CERMEP reference data ##########################
#################################################################################################
df_summary_cermep = pd.read_csv(os.path.join(cermep_summary_csv_path, summary_uptake_csv_filename))
df_summary_suvr_mean_cermep = pd.read_csv(os.path.join(cermep_summary_csv_path, summary_suvr_mean_csv_filename))
df_summary_suvr_cerebellum_cermep = pd.read_csv(os.path.join(cermep_summary_csv_path, summary_suvr_cerebellum_csv_filename))

# CERMEP TSV has a trailing delimiter in some rows; read with explicit columns to prevent field shift.
cermep_columns = [
    'participant_id', 'birth_date', 'age', 'sex', 'weight',
    'size', 'injected_dose', 'activity_in_ICV_mask', 'handedness', 'comments', '_extra'
]
df_demographics_cermep = pd.read_csv(
    cermep_demographics_filename,
    sep='\t',
    header=0,
    names=cermep_columns,
    usecols=['participant_id', 'age', 'sex', 'weight']
)
df_demographics_cermep.rename(columns={'participant_id': 'ID'}, inplace=True)
df_demographics_cermep['ID'] = df_demographics_cermep['ID'].astype(str).str.strip()
df_summary_cermep['ID'] = df_summary_cermep['ID'].astype(str).str.strip()
df_summary_suvr_mean_cermep['ID'] = df_summary_suvr_mean_cermep['ID'].astype(str).str.strip()
df_summary_suvr_cerebellum_cermep['ID'] = df_summary_suvr_cerebellum_cermep['ID'].astype(str).str.strip()
df_demographics_cermep['age'] = pd.to_numeric(df_demographics_cermep['age'], errors='coerce')
df_demographics_cermep['weight'] = pd.to_numeric(df_demographics_cermep['weight'], errors='coerce')
df_demographics_cermep['sex'] = df_demographics_cermep['sex'].astype(str).str.strip()
# Add CERMEP demographics to the summary dataframes
df_summary_cermep = pd.merge(df_summary_cermep, df_demographics_cermep, on='ID', how='left')
df_summary_suvr_mean_cermep = pd.merge(df_summary_suvr_mean_cermep, df_demographics_cermep, on='ID', how='left')
df_summary_suvr_cerebellum_cermep = pd.merge(df_summary_suvr_cerebellum_cermep, df_demographics_cermep, on='ID', how='left')

# Merge the dataframes
df_summary = pd.concat([df_summary_unsam, df_summary_controls_unsam, df_summary_cermep], ignore_index=True)
df_summary_suvr_mean = pd.concat([df_summary_suvr_mean_unsam, df_summary_suvr_mean_controls_unsam, df_summary_suvr_mean_cermep], ignore_index=True)
df_summary_suvr_cerebellum = pd.concat([df_summary_suvr_cerebellum_unsam, df_summary_suvr_cerebellum_controls_unsam, df_summary_suvr_cerebellum_cermep], ignore_index=True)

# Add the group column to the summary dataframes, after the ID column, with values 'unsam' or 'cermep' depending on the source of the data
df_summary.insert(1, 'group', ['LC' if id in df_summary_unsam['ID'].values else ('Controls' if id in df_summary_controls_unsam['ID'].values else 'Controls') for id in df_summary['ID']])
df_summary_suvr_mean.insert(1, 'group', ['LC' if id in df_summary_suvr_mean_unsam['ID'].values else ('Controls' if id in df_summary_suvr_mean_controls_unsam['ID'].values else 'Controls') for id in df_summary_suvr_mean['ID']])
df_summary_suvr_cerebellum.insert(1, 'group', ['LC' if id in df_summary_suvr_cerebellum_unsam['ID'].values else ('Controls' if id in df_summary_suvr_cerebellum_controls_unsam['ID'].values else 'Controls') for id in df_summary_suvr_cerebellum['ID']])

# Add the site column to the summary dataframes, after the group column, with values 'unsam' or 'cermep' depending on the source of the data
df_summary.insert(2, 'site', ['CEUNIM' if id in df_summary_unsam['ID'].values or id in df_summary_controls_unsam['ID'].values else 'CERMEP' for id in df_summary['ID']])
df_summary_suvr_mean.insert(2, 'site', ['CEUNIM' if id in df_summary_suvr_mean_unsam['ID'].values or id in df_summary_suvr_mean_controls_unsam['ID'].values else 'CERMEP' for id in df_summary_suvr_mean['ID']])
df_summary_suvr_cerebellum.insert(2, 'site', ['CEUNIM' if id in df_summary_suvr_cerebellum_unsam['ID'].values or id in df_summary_suvr_cerebellum_controls_unsam['ID'].values else 'CERMEP' for id in df_summary_suvr_cerebellum['ID']])

# Now move the age, sex and weight columns to be after the site column
age_col = df_summary.pop('age')
sex_col = df_summary.pop('sex')
weight_col = df_summary.pop('weight')
df_summary.insert(3, 'age', age_col)
df_summary.insert(4, 'sex', sex_col)
df_summary.insert(5, 'weight', weight_col)

# Write the merged summary dataframes to csv files  
df_summary.to_csv(os.path.join(merged_summary_csv_path, summary_uptake_csv_filename), index=False)
df_summary_suvr_mean.to_csv(os.path.join(merged_summary_csv_path, summary_suvr_mean_csv_filename), index=False)
df_summary_suvr_cerebellum.to_csv(os.path.join(merged_summary_csv_path, summary_suvr_cerebellum_csv_filename), index=False)

# Generate a table with the demographics for each group and site, with the median and standard deviation for age and weight, and the counts for sex, and save it to a csv file in the merged_summary_csv_path
demographics_summary = df_summary.groupby(['group', 'site']).agg(
    age_median=('age', 'median'),
    age_std=('age', 'std'),
    weight_median=('weight', 'median'),
    weight_iqr=('weight', lambda x: x.quantile(0.75) - x.quantile(0.25)),
    sex_counts=('sex', lambda x: x.value_counts().to_dict())
).reset_index()
demographics_summary.to_csv(os.path.join(merged_summary_csv_path, 'demographics_summary.csv'), index=False)

# Now harmonize the SUVR values for each structure in the summary_suvr_mean and summary_suvr_cerebellum dataframes, using neuroharmonize
from neuroHarmonize import harmonizationLearn
def harmonize_suvr(df_input, reference_group='CEUNIM'):
    # The harmonization model requires numeric feature values and matching row counts between the
    # SUVR matrix and the demographic covariates. Drop incomplete demographics and any non-feature
    # metadata columns that survive a previous merge (e.g. 'anonimizado').
    df = df_input.copy()
    metadata_cols = {'ID', 'group', 'site', 'age', 'sex', 'weight', 'anonimizado'}

    df['sex'] = df['sex'].astype(str).str.strip().str.upper()
    df['sex'] = df['sex'].map({'F': 0, 'M': 1, 'FEMALE': 0, 'MALE': 1})

    # Keep only the numeric regional uptake measures; metadata columns are not valid predictors.
    suvr_columns = [col for col in df.columns if col not in metadata_cols]
    feature_df = df[suvr_columns].apply(pd.to_numeric, errors='coerce')

    suvr_values = feature_df.to_numpy(dtype=float)
    covars = df[['age', 'sex', 'site']].copy() # do not include weight
    covars['site'] = covars['site'].astype(str).str.strip()
    covars = covars.rename(columns={'site': 'SITE'})
    covars['SITE'] = covars['SITE'].astype('category')
    covars['sex'] = covars['sex'].astype(float)
    covars['age'] = pd.to_numeric(covars['age'], errors='coerce')

    # Write covars for debugging
    covars.to_csv(os.path.join(merged_summary_csv_path, 'covars_for_harmonization.csv'), index=False)
    # run harmonization and store the adjusted data
    my_model, suvr_values_harmonized = harmonizationLearn(suvr_values, covars, ref_batch=reference_group)
    # Create a new dataframe with the harmonized SUVR values
    df_harmonized = df_input.copy()
    df_harmonized[suvr_columns] = suvr_values_harmonized
    return df_harmonized

# Harmonize the SUVR values for each structure in the summary_suvr_mean and summary_suvr_cerebellum dataframes
df_summary_suvr_mean_harmonized = harmonize_suvr(df_summary_suvr_mean, reference_group='CEUNIM')
df_summary_suvr_cerebellum_harmonized = harmonize_suvr(df_summary_suvr_cerebellum, reference_group='CEUNIM')
# Write the harmonized summary dataframes to csv files
df_summary_suvr_mean_harmonized.to_csv(os.path.join(merged_summary_csv_path, 'summary_suvr_mean_harmonized.csv'), index=False)
df_summary_suvr_cerebellum_harmonized.to_csv(os.path.join(merged_summary_csv_path, 'summary_suvr_cerebellum_harmonized.csv'), index=False)

# Plot the distribution of age, sex and weight for both groups, in a single figure with three subplots, one for each variable, and save the figure to the merged_summary_csv_path
import matplotlib.pyplot as plt
import seaborn as sns

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Age distribution
sns.histplot(data=df_summary, x='age', hue='group', multiple='stack', ax=axes[0])
axes[0].set_title('Age Distribution')

# Sex distribution
sns.countplot(data=df_summary, x='sex', hue='group', ax=axes[1])
axes[1].set_title('Sex Distribution')

# Weight distribution
sns.histplot(data=df_summary, x='weight', hue='group', multiple='stack', ax=axes[2])
axes[2].set_title('Weight Distribution')

plt.tight_layout()
plt.savefig(os.path.join(merged_summary_csv_path, 'demographics_distribution.png'))
plt.close()


# Plot each structure in an individual file split by group and site.
site_palette = {
    'CEUNIM': '#d73027',
    'CERMEP': '#1f78b4'
}


def sanitize_filename(name):
    allowed = set('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_')
    normalized = ''.join(ch if ch in allowed else '_' for ch in name)
    return normalized.strip('_')


def plot_structures_per_file(df, dataset_name):
    # Create a subfolder to sabe the plots for each dataset
    plots_path = os.path.join(merged_summary_csv_path, dataset_name)
    os.makedirs(plots_path, exist_ok=True)
    structure_columns = [col for col in df.columns if col not in ['ID', 'group', 'site', 'age', 'sex', 'weight']]

    for structure in structure_columns:
        plot_df = df[['group', 'site', structure]].copy()
        plot_df[structure] = pd.to_numeric(plot_df[structure], errors='coerce')
        plot_df = plot_df.dropna(subset=[structure])
        if plot_df.empty:
            continue

        fig, ax = plt.subplots(figsize=(7, 6))
        sns.boxplot(
            data=plot_df,
            x='group',
            y=structure,
            hue='site',
            palette=site_palette,
            dodge=True,
            width=0.6,
            fliersize=0,
            ax=ax
        )
        sns.stripplot(
            data=plot_df,
            x='group',
            y=structure,
            hue='site',
            palette=site_palette,
            dodge=True,
            alpha=0.65,
            size=3,
            jitter=0.18,
            linewidth=0,
            ax=ax
        )

        handles, labels = ax.get_legend_handles_labels()
        if handles and labels:
            unique = dict(zip(labels, handles))
            ax.legend(unique.values(), unique.keys(), title='Site', loc='best', frameon=False)

        if (plot_df[structure] > 0).all():
            ax.set_yscale('log')

        ax.set_title(f'{dataset_name}: {structure} by group and site')
        ax.set_xlabel('Group')
        ax.set_ylabel('SUVR')
        plt.tight_layout()

        filename = f'{dataset_name}_{sanitize_filename(structure)}.png'
        plt.savefig(os.path.join(plots_path, filename), dpi=150)
        plt.close(fig)


plot_structures_per_file(df_summary_suvr_mean, 'suvr_mean')
plot_structures_per_file(df_summary_suvr_cerebellum, 'suvr_cerebellum')
