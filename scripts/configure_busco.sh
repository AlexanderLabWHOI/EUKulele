#!/bin/bash

# A script that takes 7 positional arguments to run BUSCO on a dataset:
# 1 - the descriptive name for the sample (1 word)
# 2 - the output directory that the results should be written to
# 3 - the location of the new configuration file to be written
# 4 - the input fasta file on which the analysis should be done
# 5 - the number of CPUs available to the user
# 6 - the BUSCO database to be used

SAMPLENAME=$1
OUTPUTDIR=$2
CONFIG_LOC=$3
INPUT_FASTA=$4
CPUS=$5
BUSCO_DB=$6

URL="https://busco-data.ezlab.org/v4/data"

mkdir -p busco_downloads/lineages
if [ ! -d busco_downloads/lineages/$BUSCO_DB ]; then
    wget -nd -r --no-parent -A $BUSCO_DB.*.tar.gz $URL/lineages/ > /dev/null 2>&1
    if [ $? -eq 0 ]; then
        tar -xvzf $BUSCO_DB.*.tar.gz* || true
        rm -f $BUSCO_DB.*.tar.gz* || true # possible race condition
        if [ ! -d busco_downloads/lineages/$BUSCO_DB ]; then
            mv -f $BUSCO_DB busco_downloads/lineages/$BUSCO_DB || true
        fi
    else
        echo "Download unsuccessful"
    fi
fi

if [ -d busco_downloads/lineages/$BUSCO_DB ]; then
    exit 0
exit 1
