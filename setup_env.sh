#!/bin/bash

# Load required modules on Sherlock
module load python/3.11.0  # use a specific version, not just python/3.9.0

# Install poetry if not already installed
if ! command -v poetry &> /dev/null; then
    pip install poetry
fi

# Tell poetry to create the venv inside the project directory
# This avoids home directory quota issues
poetry config virtualenvs.in-project true

# Install dependencies
poetry install --only main

echo ""
echo "Setup complete. To activate the environment run:"
echo "  source .venv/bin/activate"
echo ""
echo "Then run the workflow with:"
echo "  cd ivim_fit"
echo "  python ./run.py <bids_dir> <output_dir> participant -np"
