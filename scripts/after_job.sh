#!/bin/bash

#SBATCH --time=1
#SBATCH -N 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 1gb

echo "hello"

# read this https://hpc.nih.gov/docs/job_dependencies.html