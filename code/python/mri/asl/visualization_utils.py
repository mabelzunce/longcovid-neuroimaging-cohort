import matplotlib.pyplot as plt
import SimpleITK as sitk
import numpy as np

def show_brain_views(image, cmap='hot', alpha=1.0, title='Image visualization', save_path=None, show=True):
    """
    Displays and optionally saves three orthogonal views (axial, coronal, sagittal) of a brain image.

    Args:
        image (sitk.Image): input image
        cmap (str): colormap to use
        alpha (float): transparency
        title (str): plot title
        save_path (str, optional): if provided, saves the figure
        show (bool, optional): if True, displays the plot
    """
    img_array = sitk.GetArrayFromImage(image)

    axial_slice = img_array.shape[0] // 2
    coronal_slice = img_array.shape[1] // 2
    sagittal_slice = img_array.shape[2] // 2

    plt.figure(figsize=(12, 4))

    # Axial
    plt.subplot(1, 3, 1)
    plt.imshow(np.flipud(img_array[axial_slice, :, :]), cmap=cmap, alpha=alpha)
    plt.title('Axial view')
    plt.axis('off')

    # Coronal
    plt.subplot(1, 3, 2)
    plt.imshow(np.flipud(img_array[:, coronal_slice, :]), cmap=cmap, alpha=alpha)
    plt.title('Coronal view')
    plt.axis('off')

    # Sagittal
    plt.subplot(1, 3, 3)
    plt.imshow(np.flipud(img_array[:, :, sagittal_slice]), cmap=cmap, alpha=alpha)
    plt.title('Sagittal view')
    plt.axis('off')

    plt.suptitle(title)

    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        print(f"Figure saved to {save_path}")

    if show:
        plt.show()
    else:
        plt.close()

def show_brain_views_with_mask(image, mask, cmap_image='gray', cmap_mask='tab20', alpha=0.3,
                               title='Brain with Mask overlay', save_path=None, show=True):
    """
    Displays and optionally saves three orthogonal views (axial, coronal, sagittal) of a brain image
    with a mask overlay.

    Args:
        image (sitk.Image): background image (e.g. anatomy or CBF)
        mask (sitk.Image): mask or atlas to overlay
        cmap_image (str): colormap for background image
        cmap_mask (str): colormap for mask overlay
        alpha (float): transparency of mask overlay
        title (str): plot title
        save_path (str, optional): if provided, saves the figure
        show (bool, optional): if True, displays the plot
    """
    img_array = sitk.GetArrayFromImage(image)
    mask_array = sitk.GetArrayFromImage(mask)

    axial_slice = img_array.shape[0] // 2
    coronal_slice = img_array.shape[1] // 2
    sagittal_slice = img_array.shape[2] // 2

    plt.figure(figsize=(12, 4))

    # Axial
    plt.subplot(1, 3, 1)
    plt.imshow(np.flipud(img_array[axial_slice, :, :]), cmap=cmap_image)
    plt.imshow(np.flipud(mask_array[axial_slice, :, :]), cmap=cmap_mask, alpha=alpha)
    plt.title('Axial view')
    plt.axis('off')

    # Coronal
    plt.subplot(1, 3, 2)
    plt.imshow(np.flipud(img_array[:, coronal_slice, :]), cmap=cmap_image)
    plt.imshow(np.flipud(mask_array[:, coronal_slice, :]), cmap=cmap_mask, alpha=alpha)
    plt.title('Coronal view')
    plt.axis('off')

    # Sagittal
    plt.subplot(1, 3, 3)
    plt.imshow(np.flipud(img_array[:, :, sagittal_slice]), cmap=cmap_image)
    plt.imshow(np.flipud(mask_array[:, :, sagittal_slice]), cmap=cmap_mask, alpha=alpha)
    plt.title('Sagittal view')
    plt.axis('off')

    plt.suptitle(title)

    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        print(f"Figure saved to {save_path}")

    if show:
        plt.show()
    else:
        plt.close()
