import nibabel as nib
import numpy as np
import pandas as pd


def process_row(row):
    metric_file = row["Metric"]
    print(f"Calculating the mean for {metric_file}")
    image_nifti = nib.load(metric_file)
    image = image_nifti.get_fdata()

    # Calculate the mean image across the fourth dimension
    mean_image = np.mean(image, axis=-1)

    # Save the 3D mean image
    mean_nifti = nib.Nifti1Image(
        mean_image, affine=image_nifti.affine
    )  # Use affine if available
    output_filename = metric_file.replace(".nii.gz", "_mean.nii.gz")
    nib.save(mean_nifti, output_filename)

    return output_filename


# Load the CSV file
df = pd.read_csv(snakemake.input.resampled)

# Create a loop to get mean for each row
output_filename_list = []
for index, row in df.iterrows():

    file_names = process_row(row)
    output_filename_list.append(file_names)

df_mean = pd.DataFrame()

df_mean["Algorithm"] = df["Algorithm"].values
df_mean["Mean_metric"] = output_filename_list

df_mean.to_csv(snakemake.output.mean, index=False)
