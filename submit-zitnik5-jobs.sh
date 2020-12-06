#!/bin/bash

#SBATCH -p scavenger
#SBATCH --mem=100 ## 100MB of RAM
#SBATCH --array=1-3

python ppi-resilience-analysis2.py --dataset zitnik --analysis 5
