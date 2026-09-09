# IVIM_fit
Preprocessing diffusion MRI data and fitting IVIM models

## How to use locally:
1. Clone the repository in Sherlock using `git clone https://github.com/PostonLab/IVIM_fit.git`. `cd` into the cloned repository, then `cd` into the ivim_fit directory within the repository. You should see a `run.py` file there.
```
git clone https://github.com/PostonLab/IVIM_fit.git
```
```
cd ivim_fit
```
2. Install poetry version 1.8.5 using `python3 -m pip install poetry==1.8.5`
```
python3 -m pip install poetry==1.8.5
```
3. Install dependencies with `poetry install`, then use `poetry shell` to activate the virtual environment.
```
poetry install
```
```
poetry shell
```
4. First run the workflow dry-run with python `./run.py <path to bids data> <path to output> participant -np`.
```
./run.py <insert path to your bids data> <insert path to your output> participant -np
```
5. If the dry run is successful, we can run the workflow with `python ./run.py <path to bids data> <path to output> participant -c all --use-singularity`. Where `-c` is for number of cores you want to use. `all` will use everything available.
```
python ./run.py <path to bids data> <path to output> participant -c all --use-singularity
```
6. When running for the first time, it will install all the singularity containers.

## How to run in Sherlock cluster at Stanford:
1. Clone the repository in Sherlock using `git clone https://github.com/PostonLab/IVIM_fit.git`. `cd` into the cloned repository, then `cd` into the ivim_fit directory within the repository. You should see a `run.py` file there.
```
git clone https://github.com/PostonLab/IVIM_fit.git
```
```
cd ivim_fit
```
(If you see the error: `fatal: unable to access 'https://github.com/PostonLab/IVIM_fit.git/': error setting certificate verify locations:  CAfile: /etc/ssl/certs/ca-certificates.crt CApath: none`, first manually point git to the CA bundle with `git config --global http.sslCAInfo /etc/ssl/certs/ca-bundle.crt`, then try `git clone` again.) 
```
git config --global http.sslCAInfo /etc/ssl/certs/ca-bundle.crt
```
1. Activate an interractive shell using `salloc`
```
salloc
```
2. Load python3 using `ml python/3.9.0`
```
ml python/3.9.0
```
3. Create a python environment `python3.9 -m venv <name>`
```
python3.9 -m venv <name>
```
4. Activate the environment `source <venv_name>/bin/activate`
```
source <venv_name>/bin/activate
```
6. Install dependencies with `poetry install`. [If poetry version 1.8.5 is not already installed, install it using `python3 -m pip install poetry==1.8.5`, then rerun `poetry install` (recommended), or `wget -O - https://install.python-poetry.org | python3 - --version 1.8.5`]. If needed install six package using `python3 -m pip install six`)
```
poetry install
```
``` title="If poetry not already installed:"
python3 -m pip install poetry==1.8.5
```
7. Deactivate the existing python environment using `deactivate`. Then use `poetry shell` to activate the new ivim virtual environment.
```
deactivate
```
```
poetry shell
```
8. Navigate to `ivim_fit` directory using `cd ivim_fit` where you will see the `run.py` file.
```
cd ivim_fit
```
9. First run the workflow dry-run with `python ./run.py <path to bids data> <path to output> participant -np`.
```
python ./run.py <path to bids data> <path to output> participant -np
```
10. If the dry run is successful, we can run the workflow with `python ./run.py <path to bids data> <path to output> participant -c all --use-singularity`. Where `-c` is for number of cores you want to use. `all` will use everything available.
```
python ./run.py <path to bids data> <path to output> participant -c all --use-singularity
```
11. If running again, first load python 3.9 using `ml python/3.9.0`, and then run `poetry shell` to activae the python environment.
```
ml python/3.9.0
```
```
poetry shell
```

## How to contribute:
1. Clone the repository with `git clone https://github.com/PostonLab/IVIM_fit.git`
```
git clone https://github.com/PostonLab/IVIM_fit.git
```
2. Checkout a branch and make changes.
3. Before doing `git add .`, run `poe quality`. This will automatically run formatting and make necessary changes.
4. If everything looks pretty, proceed with pushing changes.
5. Go to github and make a pull request for your changes.

## Acknowledgement
This project utilizes code adapted from the snakedwi pipeline ([https://github.com/akhanf/snakedwi]), an open-source project for diffusion MRI preprocessing. Their work provided a valuable foundation for building the preprocessing steps of IVIM_fit pipeline.









