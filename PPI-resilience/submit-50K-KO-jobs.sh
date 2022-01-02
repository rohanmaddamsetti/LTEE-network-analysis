#!/bin/bash

## This script submits the 50K single knockout analysis scripts to the Duke Compute Cluster.

sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset cong --analysis 6 --population Ara+1"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset cong --analysis 6 --population Ara+2"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset cong --analysis 6 --population Ara+3"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset cong --analysis 6 --population Ara+4"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset cong --analysis 6 --population Ara+5"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset cong --analysis 6 --population Ara+6"

sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset cong --analysis 6 --population Ara-1"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset cong --analysis 6 --population Ara-2"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset cong --analysis 6 --population Ara-3"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset cong --analysis 6 --population Ara-4"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset cong --analysis 6 --population Ara-5"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset cong --analysis 6 --population Ara-6"

sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset zitnik --analysis 6 --population Ara+1"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset zitnik --analysis 6 --population Ara+2"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset zitnik --analysis 6 --population Ara+3"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset zitnik --analysis 6 --population Ara+4"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset zitnik --analysis 6 --population Ara+5"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset zitnik --analysis 6 --population Ara+6"

sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset zitnik --analysis 6 --population Ara-1"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset zitnik --analysis 6 --population Ara-2"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset zitnik --analysis 6 --population Ara-3"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset zitnik --analysis 6 --population Ara-4"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset zitnik --analysis 6 --population Ara-5"
sbatch -p scavenger --mem=100 --wrap="python ppi-resilience-analysis.py --dataset zitnik --analysis 6 --population Ara-6"
