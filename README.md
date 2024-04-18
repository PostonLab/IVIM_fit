# IVIM_fit

Preprocessing diffusion MRI data and fitting IVIM models

## How to use:

  1. Install poetry
  2. Install dependencies with `poetry install`, then use `poetry shell` to activate the virtual environment.
  3. Navigate to ivim_fit directory `cd ivim_fit` where you will see the `run.py` file.
  4. First run the workflow dry-run with `./run.py <path to bids data> <path to output> participant -np`.
  5. If the dry run is successful, we can run the workflow with `./run.py <path to bids data> <path to output> participant -c all`. Where `-c` is for number of cores you want to use. `all` will use everything available.

## How to contribute:

  1. Clone the repository with `git clone https://github.com/PostonLab/IVIM_fit.git`
  2. Checkout a branch and make changes.
  3. Before doing `git add .`, run `poe quality`. This will automatically run formatting and make necessary changes. 
  4. If everything looks pretty, proceed with pushing changes. 
  5. Go to github and make a pull request for your changes. 

## Acknowledgement 

This project utilizes code adapted from the snakedwi pipeline ([https://github.com/akhanf/snakedwi]), an open-source project for diffusion MRI preprocessing. Their work provided a valuable foundation for building the preprocessing steps of IVIM_fit pipeline.