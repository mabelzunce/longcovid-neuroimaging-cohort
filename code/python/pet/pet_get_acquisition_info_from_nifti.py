import os
import nibabel as nib

# Script that reads each nifti image from the PET folder and extracts the acquisition information from the nifti header. The extracted information is then saved in a text file for each subject.
pet_path = '/home/martin/data_imaging/CovidProject/FullStudy/PET/'
nifti_path = os.path.join(pet_path, 'Nifti_statics')


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

    # Extract acquisition information from the NIfTI header
    img = nib.load(fdg_nifti_file)
    header = img.header

    # print the full header information
    print(f'Acquisition information for {subject}:')
    print(header)