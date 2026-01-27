# IVIM_fit

https://postonlab.github.io/IVIM_docs/

Preprocessing diffusion MRI data and fitting IVIM models

## How to use:

  1. Clone the repository and cd into it.
  2. Install poetry using `pip install poetry`
  2. Install dependencies with `poetry install`, then use `poetry shell` to activate the virtual environment.
  3. Navigate to ivim_fit directory `cd ivim_fit` where you will see the `run.py` file.
  4. First run the workflow dry-run with `python ./run.py <path to bids data> <path to output> participant -np`.
  5. If the dry run is successful, we can run the workflow with `python ./run.py <path to bids data> <path to output> participant -c all --use-singularity`. Where `-c` is for number of cores you want to use. `all` will use everything available.
  6. When running for the first time, it will install all the singularity containers. 

## How to run in Sherlock cluster at Stanford:

  1. Activate an interractive shell using `salloc`
  1. Load python3 using `ml python/3.9.0`
  2. Create a python environment `python3.9 -m venv <name>`
  3. Activate the environment `source py3/bin/activate`
  4. Install dependencies with `poetry install`.(If poetry is not installed, install it using `python3 -m pip install poetry`. If needed install `six` package using `python3 -m pip install six`)
  5. Deactivate the existing python environment using `deactivate`. Then use `poetry shell` to activate the new ivim virtual environment. 
  6. Navigate to `ivim_fit` directory using `cd ivim_fit` where you will see the `run.py` file.
  7. First run the workflow dry-run with `python ./run.py <path to bids data> <path to output> participant -np`.
  8. If the dry run is successful, we can run the workflow with `python ./run.py <path to bids data> <path to output> participant -c all --use-singularity`. Where `-c` is for number of cores you want to use. `all` will use everything available.
  9. If running again, first load python 3.9 using `ml python/3.9.0`, and then run `poetry shell` to activae the python environment.

## How to contribute:

  1. Clone the repository with `git clone https://github.com/PostonLab/IVIM_fit.git`
  2. Checkout a branch and make changes.
  3. Before doing `git add .`, run `poe quality`. This will automatically run formatting and make necessary changes. 
  4. If everything looks pretty, proceed with pushing changes. 
  5. Go to github and make a pull request for your changes. 

## Acknowledgement 

This project utilizes code adapted from the snakedwi pipeline ([https://github.com/akhanf/snakedwi]), an open-source project for diffusion MRI preprocessing. Their work provided a valuable foundation for building the preprocessing steps of IVIM_fit pipeline.
