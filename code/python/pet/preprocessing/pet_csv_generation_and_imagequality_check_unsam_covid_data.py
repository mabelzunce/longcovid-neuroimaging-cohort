import os
import subprocess
import sys
import shlex
import nibabel as nib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Script that reads the preprocessed data and generates multiple snapshots for image quality check. It uses the UNSAMPETQuantificationTools to preprocess the data and generate snapshots.

# Configuration
pet_only = False  # Set to True if you want to preprocess only PET data without MRI
pet_mri = True  # Set to True if you want to preprocess PET data with MRI

# Paths
pet_covid_unsam_path = '/home/martin/data_imaging/CovidProject/FullStudy/PET//'
preprocessed_path = os.path.join(pet_covid_unsam_path, 'preprocessed')
preprocessed_with_mri_path = os.path.join(pet_covid_unsam_path, 'ProcesadosCOVID_PET_MRI_09_2025')

subjects_to_exclude = ['CP0106']  # List of subjects to exclude from processing

# Output path for checked data
checked_path = os.path.join(preprocessed_path, 'checked')
os.makedirs(checked_path, exist_ok=True)
# Output path for summary csv file
summary_csv_path = os.path.join(preprocessed_path, 'rois_uptake')
os.makedirs(summary_csv_path, exist_ok=True)

# UNSAMPETQuantificationTools output subfolders:
final_results_subfolder = 'Final'
final_intensity_normalized_results_subfolder = 'Final/Normalized'
final_csv_results_subfolder = 'Final/CSV'
spatial_normalized_filename = 'spatial_normalized_pet_image.nii.gz'
spatial_intensity_normalized_mean_filename = 'intensity_normalization_mean.nii.gz'
spatial_intensity_normalized_cerebellum_filename = 'intensity_normalization_cerebellum.nii.gz'
roi_uptake_csv_filename = 'subject_normalization_values_Hammers.csv'


summary_uptake_csv_filename = 'summary_roi_uptake_values.csv'
summary_suvr_mean_csv_filename = 'summary_suvr_mean_values.csv'
summary_suvr_cerebellum_csv_filename = 'summary_suvr_cerebellum_values.csv'

# Create an empty dataframe to store the summary
df_summary = pd.DataFrame()
df_summary_suvr_mean = pd.DataFrame()
df_summary_suvr_cerebellum = pd.DataFrame()

# Get all subject folders
subject_dirs = sorted([d for d in os.listdir(preprocessed_path) if d.startswith('CP') and d not in subjects_to_exclude])
for subject in subject_dirs:
    subject_nifti_path = os.path.join(preprocessed_path, subject)
    
    # Check if the FDG NIfTI file exists
    fdg_preprocessed_nifti_file = os.path.join(subject_nifti_path, final_results_subfolder, spatial_normalized_filename)
    if not os.path.exists(fdg_preprocessed_nifti_file):
        print(f'FDG NIfTI file not found for {subject}, skipping.')
        continue

    # Read the nifti image and generate snapshots for image quality check:
    img = nib.load(fdg_preprocessed_nifti_file)
    data = img.get_fdata()

    # Generate snapshots
    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(3, 5, figure=fig)
    for i, slice_idx in enumerate(np.linspace(data.shape[2] // 3, data.shape[2] * 2 //3, 5, dtype=int)):
        ax = fig.add_subplot(gs[i // 5, i % 5])
        ax.imshow(np.flipud(np.transpose(data[:, :, slice_idx])), cmap='gray_r')
        ax.set_title(f'Slice {slice_idx}')
        ax.axis('off')
    # Now coronal slices
    for i, slice_idx in enumerate(np.linspace(data.shape[1] // 3, data.shape[1] * 2 //3, 5, dtype=int)):
        ax = fig.add_subplot(gs[(i +5) // 5, (i + 5) % 5])
        ax.imshow(np.flipud(np.transpose(data[:, slice_idx, :])), cmap='gray_r')
        ax.set_title(f'Coronal Slice {slice_idx}')
        ax.axis('off')
    # Now sagittal slices
    for i, slice_idx in enumerate(np.linspace(data.shape[0] // 3, data.shape[0] * 2 //3, 5, dtype=int)):
        ax = fig.add_subplot(gs[(i +10) // 5, (i + 10) % 5])
        ax.imshow(np.flipud(np.transpose(data[slice_idx, :, :])), cmap='gray_r')
        ax.set_title(f'Sagittal Slice {slice_idx}')
        ax.axis('off')
    plt.suptitle(f'Image Quality Check: {subject}')
    # Save the figure instead of showing it
    output_figure_path = os.path.join(checked_path, f'{subject}_image_quality_check.png')
    os.makedirs(checked_path, exist_ok=True)
    plt.savefig(output_figure_path)
    plt.close(fig)

    # Now read csv:
    csv_file_path = os.path.join(subject_nifti_path, final_csv_results_subfolder, roi_uptake_csv_filename)
    # Csv headers:
    csv_headers = ['Structure', 'n_label', 'Regional uptake mean values', 'Normalization to total brain mean value', 'Normalization to cerebellum uptake values']
    if os.path.exists(csv_file_path):
        # Read the CSV file and print the first few lines for verification
        df_uptake = pd.read_csv(csv_file_path, names=csv_headers, header=0)
        # Create a df with the ID and the mean values for each structure all in one row, use as headers
        # ID and all the names taken from the csv Structure column, and the values from the 'Regional uptake mean values' column
        df_uptake_row = pd.DataFrame(df_uptake['Regional uptake mean values'].values.reshape(1, -1), columns=df_uptake['Structure'].values)
        df_uptake_row.insert(0, 'ID', subject)
        # Append the row to the summary dataframe
        df_summary = pd.concat([df_summary, df_uptake_row], ignore_index=True)
        # Create a df with the ID and the mean values for each structure all in one row, use as headers
        # ID and all the names taken from the csv Structure column, and the values from the 'Normalization to total brain mean value' column
        df_suvr_mean_row = pd.DataFrame(df_uptake['Normalization to total brain mean value'].values.reshape(1, -1), columns=df_uptake['Structure'].values)
        df_suvr_mean_row.insert(0, 'ID', subject)
        # Append the row to the summary dataframe
        df_summary_suvr_mean = pd.concat([df_summary_suvr_mean, df_suvr_mean_row], ignore_index=True)
        # Create a df with the ID and the mean values for each structure all in one row, use as headers
        # ID and all the names taken from the csv Structure column, and the values from the 'Normalization to cerebellum uptake values' column
        df_suvr_cerebellum_row = pd.DataFrame(df_uptake['Normalization to cerebellum uptake values'].values.reshape(1, -1), columns=df_uptake['Structure'].values)
        df_suvr_cerebellum_row.insert(0, 'ID', subject)
        # Append the row to the summary dataframe
        df_summary_suvr_cerebellum = pd.concat([df_summary_suvr_cerebellum, df_suvr_cerebellum_row], ignore_index=True)
# Write csv summary file:
df_summary.to_csv(os.path.join(summary_csv_path, summary_uptake_csv_filename), index=False)
df_summary_suvr_mean.to_csv(os.path.join(summary_csv_path, summary_suvr_mean_csv_filename), index=False)
df_summary_suvr_cerebellum.to_csv(os.path.join(summary_csv_path, summary_suvr_cerebellum_csv_filename), index=False)


 

    