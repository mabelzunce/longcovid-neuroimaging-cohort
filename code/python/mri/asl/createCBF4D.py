import SimpleITK as sitk
import numpy as np
import os
import pandas as pd
from collections import Counter
import shutil

processed_images = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/ASLProcessedCOVID/derivatives/ExploreASL/Population'
output_path_image_cbf_4d =  'cbf_total_covid_final.nii.gz'
output_path_image_mask_vascular_4d =  'mask_vascular_total_covid_final.nii.gz'
output_path_csv =  'subjects_total_covid_final.csv'
output_path_cbf_images =  '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/FinalCBFImages'

if not os.path.exists(output_path_cbf_images):
    os.mkdir(output_path_cbf_images)

all_files = os.listdir(processed_images)

images = []
subjects = []
groups = []

for file in all_files:
    if file.startswith('qCBF_sub'):
        images.append(file)

images.sort()

# First image
first_image = sitk.ReadImage(os.path.join(processed_images, images[0]))
first_image_size = first_image.GetSize()

# # Crear un array 4D vacío
img_cbf_array_4d = np.zeros((*first_image_size, len(images)), dtype=np.float32)
mask_array_4d = np.zeros((*first_image_size, len(images)), dtype=np.uint8) # uint porque es binaria
#
for i,image in enumerate(images):
    subject = image[9:15]

    CBF_image= os.path.join(processed_images, image)

    if not os.path.exists(CBF_image):
        continue

    img = sitk.ReadImage(CBF_image)
    img.CopyInformation(first_image)

    if img.GetSize() != first_image.GetSize():
        raise ValueError(f"Las dimensiones de {subject} no coinciden con la primera imagen.")

    img_cbf_array_4d[:, :, :, i] = sitk.GetArrayFromImage(img) # Agrego a la imagen 4D
    shutil.copy(CBF_image, os.path.join(output_path_cbf_images, image))

    # -------------------------------
    # Leer máscara correspondiente
    # -------------------------------
    # Buscar por patrón MaskVascular_sub-{subject}_1_ASL_1.nii.gz
    mask_filename = f"MaskVascular_sub-{subject}_1_ASL_1.nii.gz"
    mask_path = os.path.join(processed_images, mask_filename)

    if not os.path.exists(mask_path):
        raise FileNotFoundError(f"Máscara no encontrada para {subject}")

    mask = sitk.ReadImage(mask_path)
    mask.CopyInformation(first_image)

    if mask.GetSize() != first_image_size:
        raise ValueError(f"Dimensión incompatible en máscara de {subject}")

    mask_array_4d[:, :, :, i] = sitk.GetArrayFromImage(mask)
    shutil.copy(mask_path, os.path.join(output_path_cbf_images, mask_filename))

    subjects.append(subject)


# # Crea el DataFrame con sujetos y grupos
df = pd.DataFrame({'Subject': subjects})
df.to_csv(os.path.join(output_path_cbf_images, output_path_csv), index=False)

# Crear Imagen CBF 4D
image_cbf_4d = sitk.GetImageFromArray(img_cbf_array_4d)
image_cbf_4d.CopyInformation(first_image)

sitk.WriteImage(image_cbf_4d, os.path.join(output_path_cbf_images, output_path_image_cbf_4d))

# Crear Imagen Mask Vascular 4D

image_mask_4d = sitk.GetImageFromArray(mask_array_4d)
image_mask_4d.CopyInformation(first_image)

sitk.WriteImage(image_mask_4d, os.path.join(output_path_cbf_images, output_path_image_mask_vascular_4d))
