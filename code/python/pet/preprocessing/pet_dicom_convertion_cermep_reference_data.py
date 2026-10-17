import os
import subprocess
from pathlib import Path

dicom_path = '/home/martin/data_imaging/PET/CERMEP/ALL/sourcedata/'
nifti_path = '/home/martin/data_imaging/PET/CERMEP/ALL/derivatives/nifti/'

# Get all subject folders
subject_dirs = sorted([d for d in os.listdir(dicom_path) if d.startswith('sub-')])

for subject in subject_dirs:
    subject_dicom_path = os.path.join(dicom_path, subject)
    subject_nifti_path = os.path.join(nifti_path, subject)
    
    # Create subject output directory
    os.makedirs(subject_nifti_path, exist_ok=True)
    
    # Process each modality folder
    modalities = ['CT', 'FDG', 'FLAIR', 'T1']
    
    for modality in modalities:
        modality_dicom_path = os.path.join(subject_dicom_path, modality)
        
        if not os.path.exists(modality_dicom_path):
            continue
        
        # Output file path
        output_nifti = os.path.join(subject_nifti_path, f'{subject}_{modality}.nii.gz')
        
        # Convert DICOM to NIfTI using dcm2niix
        try:
            subprocess.run(
                ['dcm2niix', '-o', subject_nifti_path, '-f', f'{subject}_{modality}', modality_dicom_path],
                check=True,
                capture_output=True
            )
            print(f'Successfully converted {subject}/{modality} to NIfTI')
        except subprocess.CalledProcessError as e:
            print(f'Error converting {subject}/{modality}: {e}')

