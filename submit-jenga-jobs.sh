#!/bin/bash

## submit-jenga-jobs.sh by Rohan Maddamsetti.
## Usage: sbatch submit-jenga-jobs.sh

#SBATCH --mem=4G ## 4GB of RAM
#SBATCH --array=1-1000

python play-genome-jenga.py
