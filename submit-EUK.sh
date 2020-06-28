#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH -p compute

jobname=$(cat config.yaml | grep jobname | cut -d ":" -f 2 | cut -d " " -f 2)

snakemake   \
        --jobs 200 --use-conda -s EUKulele \
        --cluster-config EUKulele.yaml --cluster "sbatch --parsable --qos=unlim --partition={cluster.queue} --job-name=$jobname.{rule}.{wildcards} --mem={cluster.mem}gb --time={cluster.time} --ntasks={cluster.threads} --nodes={cluster.nodes}"

