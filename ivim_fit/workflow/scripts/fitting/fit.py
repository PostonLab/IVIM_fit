import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
import subprocess

repo_url = "https://github.com/merelvdthiel/TF2.4_IVIM-MRI_CodeCollection.git"  # IVIM fit repo 
destination_dir = "../../lib"  # desired directory o fthe repo

#subprocess.run(["git", "clone", repo_url])

sys.path.insert(0, 'TF2.4_IVIM-MRI_CodeCollection')
from utilities.data_simulation.Download_data import download_data
import nibabel
import os

# Load the data.

# load and init b-values
bval = np.genfromtxt(snakemake.input.bval)

#load nifti
data = nibabel.load(snakemake.input.dwi)
datas = data.get_fdata()

from src.wrappers.ivim_fit import ivim_fit
from src.wrappers.OsipiBase import OsipiBase

currentdir = os.path.dirname(os.getcwd())
dirs = os.listdir('TF2.4_IVIM-MRI_CodeCollection/src/standardized')
dirs.remove('__init__.py')
print(dirs)

from src.standardized.IAR_LU_biexp import IAR_LU_biexp
from src.standardized.OGC_AmsterdamUMC_Bayesian_biexp import OGC_AmsterdamUMC_Bayesian_biexp




import nibabel as nib
import numpy as np

def apply_mask(dwi_path, mask_path):
    # Load Diffusion Weighted Image (DWI) and Brain Mask
    dwi_img = nib.load(dwi_path)
    mask_img = nib.load(mask_path)

    # Ensure that both images are in the same space
    #print(mask_img.shape)
    #print(dwi_img.shape)
    #assert dwi_img.shape == mask_img.shape, "Images are not in the same space"

    # Get the image data as numpy arrays
    dwi_data = dwi_img.get_fdata()
    mask_data = mask_img.get_fdata()

    # Convert mask into a binary format
    mask_data_bin = np.where(mask_img.get_fdata()>0, 1, 0)

    # Apply the mask to the DWI
    masked_data = np.multiply(dwi_img.dataobj[..., :],
                              mask_data_bin[..., np.newaxis])

    # Create a new Nifti image
    masked_img = nib.Nifti1Image(masked_data, dwi_img.affine, dwi_img.header)

    return masked_img, masked_data


#masked_img, masked_data = apply_mask(snakemake.input.dwi, snakemake.input.mask)
#nib.save(masked_img, "masked_dwi.nii")



"""#choose a voxel
x=60
y=60
z=30
sx, sy, sz, n_bval = datas.shape

print(datas.shape)

data_vox=np.squeeze(datas[x,y,z,:])

# normalise data
selsb = np.array(bval) == 0
print(selsb)
S0 = np.nanmean(data_vox[selsb], axis=0).astype('<f')
data_vox = data_vox / S0

print(data_vox)

direction = 6 #choose: 1, 2, 3, 4, 5, or 6
signal_1dir=data_vox[1:9]
print(signal_1dir)
signal_1dir=np.insert(signal_1dir, 0, 1)

print(signal_1dir)"""


"""

algorithm1=IAR_LU_biexp()
algorithm2=OGC_AmsterdamUMC_Bayesian_biexp()


#fit the IVIM model for 1 direction

datas = masked_data

data_1dir=datas[:,:,:,1:9] #pick the signal with the same bvec
data_1dir=np.insert(data_1dir, 0 , datas[:,:,:,0], axis=3) #add S0
#reshape data for fitting
sx, sy, sz, n_bval = data_1dir.shape
X_dw = np.reshape(data_1dir, (sx * sy * sz, n_bval))
#select only relevant values, delete background and noise, and normalise data
selsb = np.array(np.unique(bval)) == 0
S0 = np.nanmean(X_dw[:, selsb], axis=1)
S0[S0 != S0] = 0
S0=np.squeeze(S0)
valid_id = (S0 > (0.5 * np.median(S0[S0 > 0])))
data_norm = X_dw[valid_id, :]
data_norm.shape

#apply algorithm 1 to the diffusion data to obtain f, D* and D for each voxel in the slice
maps = OsipiBase.osipi_fit(algorithm1,data_norm,np.unique(bval))

f_array=maps[:,0]
Dstar_array=maps[:,1]
D_array=maps[:,2]

f_map = np.zeros([sx * sy * sz])
f_map[valid_id] = f_array[0:sum(valid_id)]
f_map = np.reshape(f_map, [sx, sy, sz])

Dstar_map = np.zeros([sx * sy * sz])
Dstar_map[valid_id] = Dstar_array[0:sum(valid_id)]
Dstar_map = np.reshape(Dstar_map, [sx, sy, sz])

D_map = np.zeros([sx * sy * sz])
D_map[valid_id] = D_array[0:sum(valid_id)]
D_map = np.reshape(D_map, [sx, sy, sz])

# save these volumes as nii.gz files
savedir=('.')
nibabel.save(nibabel.Nifti1Image(f_map, data.affine, data.header),'{folder}/f.nii.gz'.format(folder = savedir))
nibabel.save(nibabel.Nifti1Image(Dstar_map, data.affine, data.header),'{folder}/Dstar.nii.gz'.format(folder = savedir))
nibabel.save(nibabel.Nifti1Image(D_map, data.affine, data.header),'{folder}/D.nii.gz'.format(folder = savedir))

"""
with open(snakemake.output.test, "w") as file:  # "w" for write mode
    # Write your data here
    file.write("This is some text to save in the file.")