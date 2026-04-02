<div align="center">

# 🧠 IVIM_fit

**Preprocessing and biophysical modeling of diffusion MRI data**

[![Python](https://img.shields.io/badge/Python-3.9--3.11-blue?logo=python&logoColor=white)](https://www.python.org/)
[![Snakemake](https://img.shields.io/badge/Snakemake-≥7.20-green)](https://snakemake.readthedocs.io/)
[![Docker](https://img.shields.io/badge/Docker-dimu2h%2Fivim__fit-blue?logo=docker)](https://hub.docker.com/r/dimu2h/ivim_fit)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/Docs-postonlab.github.io-purple)](https://postonlab.github.io/IVIM_docs/)

</div>

---

## What is IVIM?

**Intravoxel Incoherent Motion (IVIM)** is a biophysical model applied to diffusion MRI data that separates two distinct signal components within each voxel:

- 🩸 **Perfusion** — signal from water molecules moving pseudo-randomly through the microvascular network (capillary blood flow)
- 🔬 **Diffusion** — signal from true Brownian motion of water in tissue

IVIM analysis requires **multi-shell diffusion MRI** acquired with multiple b-values, including **low b-values (b < 200 s/mm²)** that are sensitive to the fast-moving perfusion component, and **higher b-values** that capture the slower tissue diffusion component. This is in contrast to standard DTI, which typically uses a single non-zero b-value shell.

By fitting the IVIM model, this pipeline estimates three key parameters per voxel:

| Parameter | Symbol | Description |
|-----------|--------|-------------|
| Perfusion fraction | *f* | Fraction of signal from capillary blood flow |
| Pseudo-diffusion coefficient | *D\** | Rate of water movement in microvasculature |
| True diffusion coefficient | *D* | Rate of Brownian motion in tissue |

This pipeline handles full preprocessing of raw BIDS diffusion data (denoising, eddy current correction, motion correction) before fitting the IVIM model, and outputs parameter maps in standard space.

---

## 🚀 Quick Start

### Requirements
- Python 3.9–3.11
- [Singularity / Apptainer](https://apptainer.org/) (for tool containers)
- BIDS-formatted diffusion MRI data with multi-shell acquisition including low b-values

---

## 💻 Local Installation

```bash
# 1. Clone the repository
git clone https://github.com/PostonLab/IVIM_fit.git
cd IVIM_fit

# 2. Run the setup script (installs all dependencies automatically)
bash setup_env.sh

# 3. Activate the virtual environment
source .venv/bin/activate

# 4. Dry run to verify everything is working
ivim_fit <path to bids data> <path to output> participant -np

# 5. Run the full workflow
ivim_fit <path to bids data> <path to output> participant -c all --use-singularity
```

> **Note:** The `-c` flag sets the number of cores. Use `-c all` to use all available cores.
> 
> **Note:** On first run, all required Singularity containers will be downloaded automatically.
>
> **Note:** If you see a warning about `datrie` failing to install during setup, this is harmless — the pipeline works correctly without it.

---

## 🐳 Docker

For users who prefer a fully containerized environment:

```bash
# Pull the image
docker pull dimu2h/ivim_fit:latest

# Run the workflow
docker run --rm \
  -v <path to bids data>:/bids_input:ro \
  -v <path to output>:/output \
  --privileged \
  dimu2h/ivim_fit:latest \
  /bids_input /output participant -c all --use-singularity
```

---

## 🖥️ Sherlock HPC Cluster (Stanford)

### Option A — Install and run with SLURM job submission

This approach sets up the Python environment once and uses Snakemake's SLURM executor to submit each subject as a separate job, enabling parallel processing across subjects.

```bash
# 1. Start an interactive session
salloc

# 2. Load Python
ml python/3.11.0

# 3. Clone the repository
git clone https://github.com/PostonLab/IVIM_fit.git
cd IVIM_fit

# 4. Run the setup script
bash setup_env.sh

# 5. Activate the virtual environment
source .venv/bin/activate

# 6. Dry run
ivim_fit <path to bids data> <path to output> participant -np

# 7. Run with SLURM (submits each subject as a separate job)
ivim_fit <path to bids data> <path to output> participant \
    --executor slurm \
    --use-singularity \
    --jobs 50 \
    --default-resources slurm_partition=normal mem_mb=8000
```

> **Note:** The virtual environment is stored inside the project directory (`.venv`) to avoid home directory quota issues on Sherlock.
>
> **Note:** When returning to run again, just reload with `ml python/3.11.0` and activate with `source .venv/bin/activate`.

---

### Option B — Run via Apptainer (no setup required)

If you prefer not to install any dependencies, pull the container directly and submit as a SLURM batch job.

```bash
# 1. Set cache to scratch to avoid home quota issues
export APPTAINER_CACHEDIR=/scratch/$USER/apptainer_cache

# 2. Pull the container
apptainer pull ivim_fit.sif docker://dimu2h/ivim_fit:latest

# 3. Dry run
apptainer run \
  --bind <path to bids data>:/bids_input \
  --bind <path to output>:/output \
  ivim_fit.sif \
  /bids_input /output participant -np
```

Create a SLURM batch script `run_ivim.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=ivim_fit
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal

apptainer run \
  --bind <path to bids data>:/bids_input \
  --bind <path to output>:/output \
  ivim_fit.sif \
  /bids_input /output participant -c 8 --use-singularity
```

Then submit:
```bash
sbatch run_ivim.sh
```

---

## 🤝 How to Contribute

```bash
# 1. Clone the repository
git clone https://github.com/PostonLab/IVIM_fit.git

# 2. Checkout a new branch and make your changes
git checkout -b my-feature

# 3. Before committing, run the quality checks
poe quality

# 4. Push your changes and open a pull request on GitHub
```

---

## 📄 Acknowledgement

This project utilizes code adapted from the [snakedwi](https://github.com/akhanf/snakedwi) pipeline, an open-source project for diffusion MRI preprocessing. Their work provided a valuable foundation for building the preprocessing steps of the IVIM_fit pipeline.

---

<div align="center">
Made with ❤️ by the <a href="https://postonlab.github.io">Poston Lab</a> at Stanford University
</div>