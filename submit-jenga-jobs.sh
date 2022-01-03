#!/bin/bash

#SBATCH -p scavenger
#SBATCH --mem=800 ## 800MB of RAM
#SBATCH --array=1-2 ## change to 1000 for final run.

python play-genome-jenga.py
