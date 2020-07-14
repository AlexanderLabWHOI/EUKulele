#!/bin/bash

#SBATCH --qos=unlim
#SBATCH --time=5000
#SBATCH -N 3
#SBATCH --cpus-per-task 15
#SBATCH --mem 180gb
#SBATCH -p scavenger

python EUKulele-config.py --config config.yaml
