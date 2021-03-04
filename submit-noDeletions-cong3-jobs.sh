#!/bin/bash

#SBATCH -p scavenger
#SBATCH --mem=100 ## 100MB of RAM
#SBATCH --array=1-100

python ppi-resilience-analysis.py --dataset cong --analysis 3 --noDeletions True
