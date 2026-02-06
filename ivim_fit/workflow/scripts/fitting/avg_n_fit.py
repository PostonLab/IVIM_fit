import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "TF2.4_IVIM-MRI_CodeCollection")

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

# Load the label names from the TSV file
label_names_df = pd.read_csv(snakemake.input.label_lookup, sep="\t")
label_names = label_names_df["name"].values

out_dir = snakemake.params.out_dir

# load and init b-values
bval = np.genfromtxt(snakemake.input.bval)


def select_data_in_one_direction(data, bval, lower_idx, upper_idx):
    # Select a subset of the data based on the same bvec
    data_1dir = data[:, lower_idx:upper_idx]

    # Insert S0 into the subset data
    data_1dir = np.insert(data_1dir, 0, data[:, 0], axis=1)

    return data_1dir


# Reading average DWI data by ROI (In T1 space)
ROI_mean_data_csv = pd.read_csv(snakemake.input.avg_dwi_t1_space, header=None)
ROI_mean_data = ROI_mean_data_csv.values

max_bval = len(bval)
bval_count = len(np.unique(bval))

# Prepare the arguments for each iteration
indexes = [
    (i, min(i + bval_count - 1, max_bval))
    for i in range(1, max_bval, bval_count)
]

# Defining a data frame to hold paths of averaged csvs
df_paths = pd.DataFrame()

# Looping through ivim fitting algorithms
for algorithm, algorithm_name in zip(algorithms, algo_names):

    f_map = []
    Dstar_map = []
    D_map = []

    # Looping through unique directions and fit ivim model
    for ind in indexes:
        lower_idx, upper_idx = ind

        data_selected = select_data_in_one_direction(
            ROI_mean_data, bval, lower_idx, upper_idx
        )

        # Fitting ivim models for each direction
        maps = OsipiBase.osipi_fit(algorithm, data_selected, np.unique(bval))

        f_array = maps.get("f")
        f_map.append(f_array)

        Dstar_array = maps.get("D*")
        Dstar_map.append(Dstar_array)

        D_array = maps.get("D")
        D_map.append(D_array)

    f_df = pd.DataFrame(np.array(f_map), columns=label_names)
    f_df["metric"] = "f"
    Dstar_df = pd.DataFrame(np.array(Dstar_map), columns=label_names)
    Dstar_df["metric"] = "Dstar"
    D_df = pd.DataFrame(np.array(D_map), columns=label_names)
    D_df["metric"] = "D"

    # Concatenate the DataFrames vertically (row-wise)
    concatenated_df = pd.concat([f_df, Dstar_df, D_df], ignore_index=True)

    savedir = f"{out_dir}_{algorithm_name}"

    concatenated_df.to_csv(f"{savedir}/avg_n_fit.csv", index=False)

    df_paths[algorithm_name] = [f"{savedir}/avg_n_fit.csv"]


df_paths.to_csv(snakemake.output.avg_n_fitted, index=False)
