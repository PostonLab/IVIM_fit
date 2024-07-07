import os
import subprocess
import sys
from multiprocessing import Pool

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

repo_url = "https://github.com/merelvdthiel/TF2.4_IVIM-MRI_CodeCollection.git"  # IVIM fit repo
destination_dir = "."  # desired directory o fthe repo

# subprocess.run(["git", "clone", repo_url])

sys.path.insert(0, "TF2.4_IVIM-MRI_CodeCollection")

# load and init b-values
bval = np.genfromtxt(snakemake.input.bval)

from src.standardized.IAR_LU_biexp import IAR_LU_biexp
from src.standardized.OGC_AmsterdamUMC_Bayesian_biexp import (
    OGC_AmsterdamUMC_Bayesian_biexp,
)
from src.wrappers.ivim_fit import ivim_fit
from src.wrappers.OsipiBase import OsipiBase

# We can import all algorithms and select the ones we want to use using the config file.
# Only 2 is hard coded for now.
algorithm1 = IAR_LU_biexp()
algorithm2 = OGC_AmsterdamUMC_Bayesian_biexp()

algo_names = ["IAR_LU_biexp", "OGC_AmsterdamUMC_biexp_segmented"]
algorithms = [algorithm1, algorithm2]


def apply_mask(dwi_path, mask_path):
    # Load Diffusion Weighted Image (DWI) and Brain Mask
    dwi_img = nib.load(dwi_path)
    mask_img = nib.load(mask_path)

    # Convert mask into a binary format
    mask_data_bin = np.where(mask_img.get_fdata() > 0, 1, 0)

    # Apply the mask to the DWI
    masked_data = np.multiply(
        dwi_img.dataobj[..., :], mask_data_bin[..., np.newaxis]
    )

    # Create a new Nifti image
    # masked_img = nib.Nifti1Image(masked_data, dwi_img.affine, dwi_img.header)

    return masked_data, dwi_img


masked_data, dwi_img = apply_mask(snakemake.input.dwi, snakemake.input.mask)


def make_dir(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)


def preprocess_and_normalize(data, bval, lower_idx, upper_idx):
    # Select a subset of the 4D image data based on the same bvec
    data_1dir = data[:, :, :, lower_idx:upper_idx]

    # Insert S0 into the subset data
    data_1dir = np.insert(data_1dir, 0, data[:, :, :, 0], axis=3)

    # Reshape data for fitting
    sx, sy, sz, n_bval = data_1dir.shape
    X_dw = np.reshape(data_1dir, (sx * sy * sz, n_bval))

    # Select only relevant values, delete background and noise, and normalise data
    # Identify the b=0 signal
    selsb = np.array(np.unique(bval)) == 0

    # Calculate mean S0 signal, replace nan values with 0
    S0 = np.nanmean(X_dw[:, selsb], axis=1)
    S0[S0 != S0] = 0
    S0 = np.squeeze(S0)

    # Identify valid ids based on S0 signal (greater than half of its median value)
    valid_id = S0 > (0.5 * np.median(S0[S0 > 0]))

    # Normalize data based on valid id
    data_norm = X_dw[valid_id, :]

    return data_norm, valid_id, [sx, sy, sz]


# apply algorithm 1 to the diffusion data to obtain f, D* and D for each voxel in the slice


def fit_and_map_params(
    data, bval, lower_idx, upper_idx, algorithm, algorithm_name
):
    """data_norm, valid_id, dimensions = preprocess_and_normalize(
        data, bval, lower_idx, upper_idx
    )
    sx, sy, sz = dimensions[0], dimensions[1], dimensions[2]

    del data

    # Fit your model
    maps = OsipiBase.osipi_fit(algorithm, data_norm, np.unique(bval))

    # Extract each parameter's array
    f_array = maps[:, 0]
    Dstar_array = maps[:, 1]
    D_array = maps[:, 2]

    # Create parameter maps and reshape each map to the original 3D form
    f_map = np.zeros([sx * sy * sz])
    f_map[valid_id] = f_array[0 : sum(valid_id)]
    f_map = np.reshape(f_map, [sx, sy, sz])

    Dstar_map = np.zeros([sx * sy * sz])
    Dstar_map[valid_id] = Dstar_array[0 : sum(valid_id)]
    Dstar_map = np.reshape(Dstar_map, [sx, sy, sz])

    D_map = np.zeros([sx * sy * sz])
    D_map[valid_id] = D_array[0 : sum(valid_id)]
    D_map = np.reshape(D_map, [sx, sy, sz])"""

    # Generate and return tuple including unique string identifier and result
    unique_index = "slice-{}to{}".format(lower_idx, upper_idx)

    out_dir = snakemake.params.out_dir

    savedir = f"{out_dir}_{algorithm_name}/{unique_index}"

    # make_dir(savedir)

    # save these volumes as nii.gz files

    """nib.save(
        nib.Nifti1Image(f_map, dwi_img.affine, dwi_img.header),
        f"{savedir}/F.nii.gz",
    )
    nib.save(
        nib.Nifti1Image(Dstar_map, dwi_img.affine, dwi_img.header),
        f"{savedir}/Dstar.nii.gz",s
    )
    nib.save(
        nib.Nifti1Image(D_map, dwi_img.affine, dwi_img.header),
        f"{savedir}/D.nii.gz",
    )"""

    return savedir


# Create a worker function
def worker(idx_range):
    lower_idx, upper_idx = idx_range

    savedir = fit_and_map_params(
        masked_data, bval, lower_idx, upper_idx, algorithm, algorithm_name
    )

    return savedir


def process_data_in_chunks(masked_data, bval, algorithm, algorithm_name):
    max_bval = len(bval)
    bval_count = len(np.unique(bval))
    num_cores = 4  # Specify the number of cores you want to use

    # Initialize Pool with specific number of cores
    with Pool(num_cores) as p:

        # Prepare the arguments for each iteration
        indexes = [
            (i, min(i + bval_count - 1, max_bval))
            for i in range(1, max_bval, bval_count)
        ]

        # Run the process in parallel
        param_maps_list = p.map(worker, indexes)

    return param_maps_list


df = pd.DataFrame()

for algorithm, algorithm_name in zip(algorithms, algo_names):
    param_maps_list = process_data_in_chunks(
        masked_data, bval, algorithm, algorithm_name
    )
    df[algorithm_name] = param_maps_list

df.to_csv(snakemake.output.fitted, index=False)
