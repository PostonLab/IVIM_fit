import os
import sys

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn.image import load_img, new_img_like

# Load the label names from the TSV file
label_names_df = pd.read_csv(snakemake.input.label_lookup, sep="\t")
label_names = label_names_df["name"].values

# Load DWI in T1 space
dwi_path = snakemake.input.dwi_t1_space

# Load the mask image (assuming the mask image path remains unchanged)
mask_path = snakemake.input.mask
mask_img = load_img(mask_path)
mask_data = mask_img.get_fdata()

# Load the segmented image
segmented_img = load_img(snakemake.input.t1_seg)

out_dir = snakemake.params.out_dir

# load and init b-values
bval = np.genfromtxt(snakemake.input.bval)


def select_data_in_one_direction(data, bval, lower_idx, upper_idx):
    # Select a subset of the data based on the same bvec
    data_1dir = data[:, lower_idx:upper_idx]

    # Insert S0 into the subset data
    data_1dir = np.insert(data_1dir, 0, data[:, 0], axis=1)

    return data_1dir


def calculate_mean_values(masked_volume, segmented_img):
    # Function to calculate mean values within each ROI in the segmented_img

    masked_data = masked_volume.get_fdata()

    # Get unique labels (excluding background)
    unique_labels = np.unique(segmented_img.get_fdata())
    unique_labels = unique_labels[
        unique_labels != 0
    ]  # Assuming 0 is background

    mean_values = []

    # Calculate mean for each label
    for label in unique_labels:
        # Create a mask for the current label
        label_mask = segmented_img.get_fdata() == label
        label_values = masked_data[label_mask]

        # Calculate mean for this label
        if label_values.size > 0:
            mean_val = np.mean(label_values)
        else:
            mean_val = np.nan  # Assign NaN if no values

        mean_values.append(mean_val)

    return mean_values


# Load the specific DWI image
dwi_img = load_img(dwi_path)
dwi_data = dwi_img.get_fdata()

num_volumes = dwi_data.shape[-1]  # Get the number of volumes (t)

total_mean = []

for i in range(num_volumes):
    # Get the ith volume
    volume = dwi_data[..., i]

    vol_img = nib.Nifti1Image(volume, dwi_img.affine)

    # masked_volume = volume * mask_data
    masked_volume = new_img_like(
        vol_img, np.where(mask_data > 0, vol_img.get_fdata(), 0)
    )

    # Calculate mean values for the current image
    mean_values = calculate_mean_values(masked_volume, segmented_img)

    total_mean.append(mean_values)

# ROI_mean_data is a 2D array (volumes, average for each ROI)
ROI_mean_data = (np.array(total_mean)).T

# Convert the array to a pandas DataFrame
avg_dwi = pd.DataFrame(ROI_mean_data)

avg_dwi.to_csv(snakemake.output.avg_dwi, index=False, header=False)
