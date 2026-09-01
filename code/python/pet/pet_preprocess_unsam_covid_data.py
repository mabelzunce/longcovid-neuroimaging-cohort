import os
import subprocess
import sys
import shlex

# Configuration
pet_only = True  # Set to True if you want to preprocess only PET data without MRI
pet_mri = False  # Set to True if you want to preprocess PET data with MRI

# Paths
UNSAMPETQuantificationTools_path = '/home/martin/data/UNSAM/PET/UNSAMPETQuantificationTools/'
petquant_script = os.path.join(UNSAMPETQuantificationTools_path, 'main.py')
pet_template = os.path.join(UNSAMPETQuantificationTools_path, 'data/atlas/MNI152_PET_1mm.nii')

if not os.path.exists(petquant_script):
    raise FileNotFoundError(f'PET quantification script not found: {petquant_script}')

if not os.path.exists(pet_template):
    raise FileNotFoundError(f'PET template not found: {pet_template}')

pet_path = '/home/martin/data_imaging/CovidProject/FullStudy/PET/'
nifti_path = os.path.join(pet_path, 'Nifti_statics')
preprocessed_path = os.path.join(pet_path, 'preprocessed')
preprocessed_with_mri_path = os.path.join(pet_path, 'preprocessed_with_mri')

# Get all subject files
subject_files = sorted([f for f in os.listdir(nifti_path) if f.startswith('CP') and f.endswith('.nii.gz')])

for pet_file in subject_files:
    # Filenames have format CP0023_s_fdgpet.nii.gz
    # Get subject ID from the filename, splitting strs by '_' and taking the first part
    pet_file_split = pet_file.split('_')
    subject = pet_file_split[0]
    
    
    # Check if the FDG NIfTI file exists
    fdg_nifti_file = os.path.join(nifti_path, pet_file)
    if not os.path.exists(fdg_nifti_file):
        print(f'FDG NIfTI file not found for {subject}, skipping.')
        continue

    if pet_only:
        # output path: (The path is created by the pipeline with the -s name)
        #subject_preprocessed_path = os.path.join(preprocessed_path, subject)
        # Create subject output directory
        #os.makedirs(subject_preprocessed_path, exist_ok=True)
        # Call the preprocessing script for the subject using UNSAMPETQuantificationTools: "PETQuantification [-h] [-m MRI] [-p PET] [-o OUTPUT] [-s SUBJECT]
        #              [-t PET_TEMPLATE] [-f FREESURFER_DIR] [--no-flirt]
        #              [--no-ant] [--no-freesurfer]""
        # First without MRI
        try:
            cmd = [
                sys.executable,
                petquant_script,
                '-p', fdg_nifti_file,
                '-o', preprocessed_path,
                '-s', subject,
                '-t', pet_template,
                '--no-flirt'
            ]
            print(f'Running command for {subject}: {shlex.join(cmd)}')

            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True
            )
            # print the output of the command
            print(f'Output for {subject}:\n{result.stdout.decode("utf-8", errors="ignore")}')
            print(f'Successfully preprocessed {subject}')
        except subprocess.CalledProcessError as e:
            print(f'Error preprocessing {subject}: {e}')
            if e.stderr:
                print(e.stderr.decode('utf-8', errors='ignore'))

    # Then with MRI
    if pet_mri:
        # Check if the T1 NIfTI file exists
        t1_nifti_file = os.path.join(subject_nifti_path, 'anat', f'{subject}_T1w.nii.gz')
        if not os.path.exists(t1_nifti_file):
            print(f'T1 NIfTI file not found for {subject}, skipping.')
            continue
        subject_preprocessed_with_mri_path = os.path.join(preprocessed_with_mri_path, subject)
        os.makedirs(subject_preprocessed_with_mri_path, exist_ok=True)
        try:
            cmd = [
                sys.executable,
                petquant_script,
                '-m', t1_nifti_file,
                '-p', fdg_nifti_file,
                '-o', subject_preprocessed_with_mri_path,
                '-s', subject,
                '-t', pet_template,
                '--no-flirt'
            ]
            print(f'Running command for {subject}: {shlex.join(cmd)}')

            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True
            )
            # print the output of the command
            print(f'Output for {subject}:\n{result.stdout.decode("utf-8", errors="ignore")}')
            print(f'Successfully preprocessed {subject}')
        except subprocess.CalledProcessError as e:
            print(f'Error preprocessing {subject}: {e}')
            if e.stderr:
                print(e.stderr.decode('utf-8', errors='ignore'))

    