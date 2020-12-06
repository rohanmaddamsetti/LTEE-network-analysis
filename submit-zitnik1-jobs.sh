#!/bin/bash

#SBATCH -p scavenger
#SBATCH --mem=100 ## 100MB of RAM
#SBATCH --array=1-100

python ppi-resilience-analysis2.py --dataset zitnik --analysis 1
