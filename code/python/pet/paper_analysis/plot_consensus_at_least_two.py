

import SimpleITK as sitk
import numpy as np
import pandas as pd
from nilearn import plotting, datasets, surface
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LinearSegmentedColormap
import os

consensus_rois_path = os.path.join('./results/pet/paper_analysis/', 'consensus/')
full_consensus_csv = os.path.join(consensus_rois_path, "consensus_rois.csv")
two_consensus_csv = os.path.join(consensus_rois_path, "consensus_at_least_two_datasets_rois.csv")

path_mni_152 = os.path.join("./data/pet/atlases",
    "MNI152_T1_1mm_brain.nii.gz"
)

path_hammers = os.path.join(
    "./data/pet/atlases",
    "Hammers_mith-n30r95-MaxProbMap-gm-MNI152-SPM12.nii.gz"
)
labels_hammers_csv = os.path.join(
    "./data/pet/atlases",
    "labels_Hammers.csv"
)

plots_dir = os.path.join(consensus_rois_path, "plots")
if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)


def resample_sitk(input_image, reference_image):
    """
    Resamples a SimpleITK image to match the spacing,
    size, and orientation of a reference image.

    Parameters
    ----------
    input_image  : SimpleITK.Image
        The input image to be resampled.
    reference_image : SimpleITK.Image
        The reference image whose spacing, size,
        and orientation will be used for resampling.

    Returns
    ---------
    resampled_image: simpleITK image object

    """
    # Create an identity transform
    identity_transform = sitk.Transform()

    # Define the interpolator to use
    interpolator = sitk.sitkLinear

    # Define the default pixel value for points outside the input image
    defaultPixelValue = 0

    # Define the output pixel type
    outputPixelType = sitk.sitkFloat32

    # Resample the input image
    resampled_image = sitk.Resample(input_image, reference_image, identity_transform, interpolator, defaultPixelValue,
                                    outputPixelType)
    return resampled_image











# Main
# read atlases
hammers_image = sitk.ReadImage(path_hammers)
mni_image = sitk.ReadImage(path_mni_152)
# resample hammers atlas to MNI space
resampled_hammers = resample_sitk(hammers_image, mni_image)
array_hammers = sitk.GetArrayFromImage(resampled_hammers)

# read hammer label
labels_hammers_df = pd.read_csv(labels_hammers_csv)
label_name_to_id = dict(zip(labels_hammers_df['structure'], labels_hammers_df['n_label']))


# Read the consensus ROIs CSV file
consensus_df = pd.read_csv(two_consensus_csv)
labels = consensus_df['ROI'].values

# Plot neurocloud effect sizes
array_zmap = np.zeros_like(array_hammers, dtype=np.float32)
effect_sizes = consensus_df['neurocloud_effect'].values
for label, effect_size in zip(labels, effect_sizes):
    # get the label id from the label name:
    label_id = label_name_to_id.get(label, None)
    if label_id is not None:
        array_zmap[array_hammers == label_id] = effect_size
zmap_image = sitk.GetImageFromArray(array_zmap)
zmap_image.CopyInformation(resampled_hammers)
# Changed output format from .png to .nii.gz to support float data
path_zmap = os.path.join(plots_dir, "consensus_at_least_two_map_neurocloud_effect_size.nii.gz")
sitk.WriteImage(zmap_image, path_zmap)
cmap_blue_white_red = LinearSegmentedColormap.from_list(
    "blue_white_red",
    [
        (0, 0, 1),   # Azul (negativo)
        (1, 1, 1),   # Blanco (0)
        (1, 0, 0)    # Rojo (positivo)
    ]
)
path_plot_zscores = os.path.join(
    plots_dir,
    "consensus_at_least_two_with_neurocloud_effect.png"
)
plotting.plot_stat_map(
    path_zmap,
    bg_img=path_mni_152,
    cmap=cmap_blue_white_red,
    display_mode="ortho",
    cut_coords=(18, 4, 19),
    threshold=0.0,
    vmax=1.5,
    vmin=-1.5,
    black_bg=False,
    annotate=False,
    dim=0,
    alpha=1.0
)

plt.savefig(
    path_plot_zscores,
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# Plot ceunim effect sizes
array_zmap = np.zeros_like(array_hammers, dtype=np.float32)
effect_sizes = consensus_df['ceunim_beta'].values # we use the beta values for ceunim effect sizes
for label, effect_size in zip(labels, effect_sizes):
    # get the label id from the label name:
    label_id = label_name_to_id.get(label, None)
    if label_id is not None:
        array_zmap[array_hammers == label_id] = effect_size
zmap_image = sitk.GetImageFromArray(array_zmap)
zmap_image.CopyInformation(resampled_hammers)
# Changed output format from .png to .nii.gz to support float data
path_zmap = os.path.join(plots_dir, "consensus_at_least_two_map_ceunim_effect_size.nii.gz")
sitk.WriteImage(zmap_image, path_zmap)
path_plot_zscores = os.path.join(
    plots_dir,
    "consensus_at_least_two_with_ceunim_effect.png"
)
plotting.plot_stat_map(
    path_zmap,
    bg_img=path_mni_152,
    cmap=cmap_blue_white_red,
    display_mode="ortho",
    cut_coords=(18, 4, 19),
    threshold=0.0,
    vmax=0.4,
    vmin=-0.4,
    black_bg=False,
    annotate=False,
    dim=0,
    alpha=1.0
)

plt.savefig(
    path_plot_zscores,
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# Plot cermep effect sizes
array_zmap = np.zeros_like(array_hammers, dtype=np.float32)
effect_sizes = consensus_df['cermep_beta'].values # we use the beta values for cermep effect sizes
for label, effect_size in zip(labels, effect_sizes):
    # get the label id from the label name:
    label_id = label_name_to_id.get(label, None)
    if label_id is not None:
        array_zmap[array_hammers == label_id] = effect_size
zmap_image = sitk.GetImageFromArray(array_zmap)
zmap_image.CopyInformation(resampled_hammers)
# Changed output format from .png to .nii.gz to support float data
path_zmap = os.path.join(plots_dir, "consensus_at_least_two_map_cermep_effect_size.nii.gz")
sitk.WriteImage(zmap_image, path_zmap)
path_plot_zscores = os.path.join(
    plots_dir,
    "consensus_at_least_two_with_cermep_effect.png"
)
plotting.plot_stat_map(
    path_zmap,
    bg_img=path_mni_152,
    cmap=cmap_blue_white_red,
    display_mode="ortho",
    cut_coords=(18, 4, 19),
    threshold=0.0,
    vmax=0.4,
    vmin=-0.4,
    black_bg=False,
    annotate=False,
    dim=0,
    alpha=1.0
)

plt.savefig(
    path_plot_zscores,
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# Now plot a plot_glass_brain for ceunim effect sizes
path_plot_glass_brain = os.path.join(
    plots_dir,
    "consensus_at_least_two_with_ceunim_effect_glass_brain.png"
)
plotting.plot_glass_brain(
    path_zmap,
    display_mode="lyrz",
    colorbar=True,
    cmap=cmap_blue_white_red,
    threshold=0.0,
    vmax=0.4,
    vmin=-0.4,
    black_bg=False,
    plot_abs=False
)
plt.savefig(
    path_plot_glass_brain,
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# Now plot a plot_glass_brain for cermep effect sizes
path_plot_glass_brain = os.path.join(
    plots_dir,
    "consensus_at_least_two_with_cermep_effect_glass_brain.png"
)
plotting.plot_glass_brain(
    path_zmap,
    display_mode="lyrz",
    colorbar=True,
    cmap=cmap_blue_white_red,
    threshold=0.0,
    vmax=0.4,
    vmin=-0.4,
    black_bg=False,
    plot_abs=False
)
plt.savefig(
    path_plot_glass_brain,
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# Plot a surface representation of the CEUNIM effect sizes
from neuromaps.datasets import fetch_fslr
from neuromaps.transforms import mni152_to_fslr
from surfplot import Plot

# Project Hammers ROI labels to fsLR 32k without averaging integer labels
lh_gii, rh_gii = mni152_to_fslr(
    path_hammers, 
    fslr_density="32k", 
    method="nearest"
)
lh_parc = lh_gii.agg_data()
rh_parc = rh_gii.agg_data()
surfaces = fetch_fslr(density="32k")
lh_surf, rh_surf = surfaces["inflated"]

# Map each projected ROI to its CEUNIM effect size
ceunim_effects_by_label = {}
for label, effect_size in zip(labels, consensus_df["ceunim_beta"].values):
    label_id = label_name_to_id.get(label)
    if label_id is not None:
        ceunim_effects_by_label[label_id] = effect_size

lh_effects = np.full(lh_parc.shape, np.nan, dtype=np.float32)
rh_effects = np.full(rh_parc.shape, np.nan, dtype=np.float32)
for label_id, effect_size in ceunim_effects_by_label.items():
    lh_effects[lh_parc == label_id] = effect_size
    rh_effects[rh_parc == label_id] = effect_size

# Initialize surfplot and render effect sizes with the volume-plot color scale
p = Plot(lh_surf, rh_surf, views=["lateral", "medial"])
p.add_layer(
    {"left": lh_effects, "right": rh_effects},
    cmap=cmap_blue_white_red,
    color_range=(-0.4, 0.4),
    cbar=True,
    cbar_label="CEUNIM effect size",
    zero_transparent=False
)

# Add outlines to retain the Hammers ROI boundaries
#p.add_layer(
#    {"left": lh_parc, "right": rh_parc},
#    cmap="gray",
#    as_outline=True,
#    cbar=False
#)

fig = p.build()
fig.show()
fig.savefig(
    os.path.join(plots_dir, "consensus_at_least_two_with_ceunim_effect_surface.png"),
    dpi=300,
    bbox_inches="tight"
)