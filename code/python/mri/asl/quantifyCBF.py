import SimpleITK as sitk
import os
import pandas as pd
import shutil

from statsmodels.sandbox.distributions.genpareto import quant

import image_utils
import quantification_utils as quant
import visualization_utils as viz

processed_images = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/ASLProcessedCOVID/derivatives/ExploreASL/Population'
path_atlas_files = '/home/sol/BrainTools/ExploreASL/External/Atlases'
output_path_cbf_images =  '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/ASL_BIDS/FinalCBFImages'

results_dir = 'ResultsQuantification'

name_Hammers_atlas = 'Hammers.nii.gz'
name_lobes_atlas = 'MNI_Structural.nii.gz'

name_Hammers_atlas_labels = 'Hammers.tsv'
name_lobes_atlas_labels = 'MNI_Structural.tsv'

# Create dir for results
results_paths = os.path.join(output_path_cbf_images, results_dir)
if not os.path.exists(results_paths):
    os.mkdir(results_paths)

# Read Atlas images and labels
atlas_hammers_img = image_utils.read_image(os.path.join(path_atlas_files, name_Hammers_atlas))
atlas_lobes_img = image_utils.read_image(os.path.join(path_atlas_files, name_lobes_atlas))

atlas_array_hammers = sitk.GetArrayFromImage(atlas_hammers_img)
atlas_array_lobes  = sitk.GetArrayFromImage(atlas_hammers_img)

df_atlas_hammers = pd.read_csv(os.path.join(path_atlas_files, name_Hammers_atlas_labels),  sep='\t', header=None)
df_atlas_lobes = pd.read_csv(os.path.join(path_atlas_files, name_lobes_atlas_labels), sep='\t', header=None)

# Initialize results dataframes
cbf_df_hammers = pd.DataFrame()
scov_df_hammers = pd.DataFrame()
mad_df_hammers = pd.DataFrame()
madnorm_df_hammers = pd.DataFrame()

# Initialize dataframes for Lobes
cbf_df_lobes = pd.DataFrame()
scov_df_lobes = pd.DataFrame()
mad_df_lobes = pd.DataFrame()
madnorm_df_lobes = pd.DataFrame()

# Initialize results dataframes for UNMASKED (whole brain)
cbf_df_hammers_unmasked = pd.DataFrame()
scov_df_hammers_unmasked = pd.DataFrame()
mad_df_hammers_unmasked = pd.DataFrame()
madnorm_df_hammers_unmasked = pd.DataFrame()

cbf_df_lobes_unmasked = pd.DataFrame()
scov_df_lobes_unmasked = pd.DataFrame()
mad_df_lobes_unmasked = pd.DataFrame()
madnorm_df_lobes_unmasked = pd.DataFrame()

# Flatten si es una fila única con todos los labels
labels_list_hammers  = df_atlas_hammers.values.flatten().tolist()
labels_list_lobes  = df_atlas_lobes.values.flatten().tolist()

labels_dict_hammers = {i+1: label for i, label in enumerate(labels_list_hammers)}
labels_dict_lobes = {i+1: label for i, label in enumerate(labels_list_lobes)}

all_files = os.listdir(processed_images)
cbf_images = []



# List all files
for file in all_files:
    if file.startswith('qCBF_sub'):
        cbf_images.append(file)

# First image
first_image_cbf = image_utils.read_image(os.path.join(processed_images, cbf_images[0]))
first_image_size = first_image_cbf.GetSize()

# Iterate across all subjects
for i,cbf_image in enumerate(cbf_images):
    subject = cbf_image[9:15] # Extract CP
    mask_filename = f"MaskVascular_sub-{subject}_1_ASL_1.nii.gz"
    image_masked_path = os.path.join(results_paths, f"Image_Masked-{subject}.nii.gz")
    image_masked_visualization_path = os.path.join(results_paths, f"Image_Masked-{subject}.png")

    image_atlas_Hammers_visualization_path = os.path.join(results_paths, f"Image_Atlas_Hammers-{subject}.png")
    image_atlas_lobes_visualization_path = os.path.join(results_paths, f"Image_Atlas_lobes-{subject}.png")



    # CBF and Mask File
    mask_image_path = os.path.join(processed_images, mask_filename)
    CBF_image_path = os.path.join(processed_images, cbf_image)

    # Read Images
    cbf_img_cbf = image_utils.read_image(CBF_image_path)
    mask_img = image_utils.read_image(mask_image_path)

    # Check dimension
    image_utils.check_same_size(cbf_img_cbf, first_image_cbf, subject, "First image CBF")
    image_utils.check_same_size(mask_img, first_image_cbf, subject, "First image CBF")

    # Apply mask
    image_masked = image_utils.multiply_mask_image(cbf_img_cbf, mask_img)

    image_utils.write_image(image_masked, image_masked_path)
    viz.show_brain_views(image_masked, cmap='hot', save_path=image_masked_visualization_path, show=False)
    viz.show_brain_views_with_mask(image_masked, atlas_hammers_img, cmap_image='gray', cmap_mask='tab20', alpha=0.4,
                                   title='CBF with Atlas overlay', show=False, save_path=image_atlas_Hammers_visualization_path)

    viz.show_brain_views_with_mask(image_masked, atlas_lobes_img, cmap_image='gray', cmap_mask='tab20', alpha=0.4,
                                   title='CBF with Atlas overlay', show=False, save_path=image_atlas_lobes_visualization_path)

    # Quantify CBF and add to a dataframe
    image_masked_array = sitk.GetArrayFromImage(image_masked)

    # Quantify Hammers
    results_hammers = quant.quantify_all_regions(image_masked_array, atlas_array_hammers, labels_dict_hammers, subject,
                                           "Hammers")

    cbf_df_hammers = pd.concat([cbf_df_hammers, pd.DataFrame(results_hammers['cbf'])], ignore_index=True)
    scov_df_hammers = pd.concat([scov_df_hammers, pd.DataFrame(results_hammers['scov'])], ignore_index=True)
    mad_df_hammers = pd.concat([mad_df_hammers, pd.DataFrame(results_hammers['mad'])], ignore_index=True)
    madnorm_df_hammers = pd.concat([madnorm_df_hammers, pd.DataFrame(results_hammers['madnorm'])], ignore_index=True)

    # Quantify Lobes
    results_lobes = quant.quantify_all_regions(image_masked_array, atlas_array_lobes, labels_dict_lobes, subject, "Lobes")

    cbf_df_lobes = pd.concat([cbf_df_lobes, pd.DataFrame(results_lobes['cbf'])], ignore_index=True)
    scov_df_lobes = pd.concat([scov_df_lobes, pd.DataFrame(results_lobes['scov'])], ignore_index=True)
    mad_df_lobes = pd.concat([mad_df_lobes, pd.DataFrame(results_lobes['mad'])], ignore_index=True)
    madnorm_df_lobes = pd.concat([madnorm_df_lobes, pd.DataFrame(results_lobes['madnorm'])], ignore_index=True)

    # Quantify unmasked image
    cbf_img_array = sitk.GetArrayFromImage(cbf_img_cbf)

    # Hammers
    results_hammers_unmasked = quant.quantify_all_regions(cbf_img_array, atlas_array_hammers, labels_dict_hammers,
                                                          subject, "Hammers")
    cbf_df_hammers_unmasked = pd.concat([cbf_df_hammers_unmasked, pd.DataFrame(results_hammers_unmasked['cbf'])],
                                        ignore_index=True)
    scov_df_hammers_unmasked = pd.concat([scov_df_hammers_unmasked, pd.DataFrame(results_hammers_unmasked['scov'])],
                                         ignore_index=True)
    mad_df_hammers_unmasked = pd.concat([mad_df_hammers_unmasked, pd.DataFrame(results_hammers_unmasked['mad'])],
                                        ignore_index=True)
    madnorm_df_hammers_unmasked = pd.concat(
        [madnorm_df_hammers_unmasked, pd.DataFrame(results_hammers_unmasked['madnorm'])], ignore_index=True)

    # Lobes
    results_lobes_unmasked = quant.quantify_all_regions(cbf_img_array, atlas_array_lobes, labels_dict_lobes, subject,
                                                        "Lobes")
    cbf_df_lobes_unmasked = pd.concat([cbf_df_lobes_unmasked, pd.DataFrame(results_lobes_unmasked['cbf'])],
                                      ignore_index=True)
    scov_df_lobes_unmasked = pd.concat([scov_df_lobes_unmasked, pd.DataFrame(results_lobes_unmasked['scov'])],
                                       ignore_index=True)
    mad_df_lobes_unmasked = pd.concat([mad_df_lobes_unmasked, pd.DataFrame(results_lobes_unmasked['mad'])],
                                      ignore_index=True)
    madnorm_df_lobes_unmasked = pd.concat([madnorm_df_lobes_unmasked, pd.DataFrame(results_lobes_unmasked['madnorm'])],
                                          ignore_index=True)

# Save Hammers
cbf_df_hammers.to_csv(os.path.join(results_paths, 'CBF_Hammers.csv'), index=False)
scov_df_hammers.to_csv(os.path.join(results_paths, 'SCOV_Hammers.csv'), index=False)
mad_df_hammers.to_csv(os.path.join(results_paths, 'MAD_Hammers.csv'), index=False)
madnorm_df_hammers.to_csv(os.path.join(results_paths, 'MADNorm_Hammers.csv'), index=False)

# Save Lobes
cbf_df_lobes.to_csv(os.path.join(results_paths, 'CBF_Lobes.csv'), index=False)
scov_df_lobes.to_csv(os.path.join(results_paths, 'SCOV_Lobes.csv'), index=False)
mad_df_lobes.to_csv(os.path.join(results_paths, 'MAD_Lobes.csv'), index=False)
madnorm_df_lobes.to_csv(os.path.join(results_paths, 'MADNorm_Lobes.csv'), index=False)

# Save UNMASKED results
cbf_df_hammers_unmasked.to_csv(os.path.join(results_paths, 'CBF_Hammers_UNMASKED.csv'), index=False)
scov_df_hammers_unmasked.to_csv(os.path.join(results_paths, 'SCOV_Hammers_UNMASKED.csv'), index=False)
mad_df_hammers_unmasked.to_csv(os.path.join(results_paths, 'MAD_Hammers_UNMASKED.csv'), index=False)
madnorm_df_hammers_unmasked.to_csv(os.path.join(results_paths, 'MADNorm_Hammers_UNMASKED.csv'), index=False)

cbf_df_lobes_unmasked.to_csv(os.path.join(results_paths, 'CBF_Lobes_UNMASKED.csv'), index=False)
scov_df_lobes_unmasked.to_csv(os.path.join(results_paths, 'SCOV_Lobes_UNMASKED.csv'), index=False)
mad_df_lobes_unmasked.to_csv(os.path.join(results_paths, 'MAD_Lobes_UNMASKED.csv'), index=False)
madnorm_df_lobes_unmasked.to_csv(os.path.join(results_paths, 'MADNorm_Lobes_UNMASKED.csv'), index=False)
