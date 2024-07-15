import nibabel as nib
import numpy as np
import pandas as pd

out_dir = snakemake.params.out_dir
df_concat = pd.DataFrame()


def concat_nii(nifti_images, algorithm, metric):
    # Load all Nifti files into a Python list and concatenate along the 4th dimension
    data = [np.expand_dims(img.get_fdata(), axis=3) for img in nifti_images]
    data_concat = np.concatenate(data, axis=3)

    # Save new 4D Image
    img_concat = nib.Nifti1Image(data_concat, nifti_images[0].affine)

    file_name = f"{out_dir}_{algorithm}/concatted_{metric}.nii.gz"
    nib.save(img_concat, file_name)

    return file_name


def concat_and_save_4D_nii(csv_path):
    # Load in csv file
    df = pd.read_csv(csv_path)

    # Loop over our algorithms
    for algorithm in df.columns:

        D_nifti_images = [
            nib.load(fp + "/D.nii.gz") for fp in df[algorithm].dropna()
        ]
        Dstar_nifti_images = [
            nib.load(fp + "/Dstar.nii.gz") for fp in df[algorithm].dropna()
        ]
        F_nifti_images = [
            nib.load(fp + "/F.nii.gz") for fp in df[algorithm].dropna()
        ]

        D_file = concat_nii(D_nifti_images, algorithm, "D")
        Dstar_file = concat_nii(Dstar_nifti_images, algorithm, "Dstar")
        F_file = concat_nii(F_nifti_images, algorithm, "F")

        file_paths = [D_file, Dstar_file, F_file]
        print(algorithm)
        print(file_paths)

        df_concat[algorithm] = file_paths


concat_and_save_4D_nii(snakemake.input.fitted_dwi)
df_concat.to_csv(snakemake.output.concatted, index=False)
