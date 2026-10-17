import os
import shutil

input_dir = r"D:/LongCOVID/PET/Nifti/"
quantification_dir = "E:/ProcesadasPETSol/ProcesadosCOVID_PET_MRI_09_2025/"
output_dir = r"D:/LongCOVID/PET/Nifti_statics/"

os.makedirs(output_dir, exist_ok=True)

for subject_id in os.listdir(quantification_dir):
    subject_path = os.path.join(quantification_dir, subject_id)
    if os.path.isdir(subject_path):
        nifti_path = os.path.join(subject_path, "nifti")
        if os.path.exists(nifti_path) and os.path.isdir(nifti_path):
            # Look for PET file
            pet_file = f"{subject_id}_rec-acstat-PET.nii.gz"
            pet_src = os.path.join(nifti_path, pet_file)
            
            # Look for T1 file
            t1_file = f"{subject_id}_t1mprage.nii.gz"
            t1_src = os.path.join(nifti_path, t1_file)
            
            # Copy PET file if exists
            if os.path.exists(pet_src):
                pet_dst = os.path.join(output_dir, f"{subject_id}_FDG_pet.nii.gz")
                shutil.copy2(pet_src, pet_dst)
                print(f"Copied PET: {pet_src} to {pet_dst}")
            else:
                print(f"PET file not found for subject {subject_id}: {pet_file}")
            
            # Copy T1 file if exists
            if os.path.exists(t1_src):
                t1_dst = os.path.join(output_dir, f"{subject_id}_t1mprage.nii.gz")
                shutil.copy2(t1_src, t1_dst)
                print(f"Copied T1: {t1_src} to {t1_dst}")
            else:
                print(f"T1 file not found for subject {subject_id}: {t1_file}")
        else:
            print(f"No nifti folder found for subject {subject_id}")