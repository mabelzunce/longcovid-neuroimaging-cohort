import SimpleITK as sitk
import os

def read_image(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"No se encontró el archivo: {path}")
    return sitk.ReadImage(path)

def write_image(image, path):
    sitk.WriteImage(image, path)

def multiply_mask_image(image, mask):
    mask_bin = sitk.Cast(mask > 0, image.GetPixelID())
    return sitk.Multiply(image, mask_bin)

def copy_image_information(target_img, reference_img):
    """
    Copies spatial information (origin, spacing, direction) from the reference image to the target image.

    Args:
        target_img (sitk.Image): image to which the information will be copied
        reference_img (sitk.Image): reference image providing the spatial information

    Returns:
        sitk.Image: target image with updated spatial information
    """
    target_img.CopyInformation(reference_img)
    return target_img

def check_same_size(img1, img2, name1="Image 1", name2="Image 2"):
    """
    Checks that two images have the same size. Raises an error if they do not.

    Args:
        img1 (sitk.Image): first image to compare
        img2 (sitk.Image): second image to compare
        name1 (str): name of the first image (for error message)
        name2 (str): name of the second image (for error message)

    Raises:
        ValueError: if the images have different sizes
    """
    if img1.GetSize() != img2.GetSize():
        raise ValueError(f"The dimensions of {name1} do not match {name2}.")
