import numpy as np

def quantify_region(image_array, region_mask):
    """
    Calculates mean, SCOV, MAD, and MAD normalized for a region.

    Args:
        image_array (np.ndarray): masked image as numpy array
        region_mask (np.ndarray): boolean mask for the region

    Returns:
        tuple: (mean_cbf, scov, mad, madnorm)
    """
    values = image_array[region_mask]
    values = values[~np.isnan(values)]

    if len(values) > 0:
        mean_cbf = np.mean(values)
        scov = np.std(values) / mean_cbf if mean_cbf != 0 else np.nan
        mad = np.median(np.abs(values - np.median(values)))
        madnorm = mad / np.median(values) if np.median(values) != 0 else np.nan
    else:
        mean_cbf, scov, mad, madnorm = [np.nan] * 4

    return mean_cbf, scov, mad, madnorm


def quantify_all_regions(image_array, atlas_array, labels_dict, subject, atlas_name=None):
    """
    Quantifies all regions in an atlas for a given subject's image.

    Args:
        image_array (np.ndarray): masked image array
        atlas_array (np.ndarray): atlas image array
        labels_dict (dict): {label_number: region_name}
        subject (str): subject ID
        atlas_name (str, optional): name of the atlas (e.g. 'Hammers', 'Lobes')

    Returns:
        dict of lists: {'cbf':[], 'scov':[], 'mad':[], 'madnorm':[]}
    """
    results = {'cbf': [], 'scov': [], 'mad': [], 'madnorm': []}

    for label_value, region_name in labels_dict.items():
        region_mask = (atlas_array == label_value)
        mean_cbf, scov, mad, madnorm = quantify_region(image_array, region_mask)

        base_dict = {'subject': subject, 'region': region_name}
        if atlas_name:
            base_dict['atlas'] = atlas_name

        results['cbf'].append({**base_dict, 'cbf': mean_cbf})
        results['scov'].append({**base_dict, 'scov': scov})
        results['mad'].append({**base_dict, 'mad': mad})
        results['madnorm'].append({**base_dict, 'madnorm': madnorm})

    return results
