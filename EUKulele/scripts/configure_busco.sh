#!/bin/bash

# A script that takes 7 positional arguments to run BUSCO on a dataset:
# 1 - the descriptive name for the sample (1 word)
# 2 - the output directory that the results should be written to
# 3 - the location of the static/template configuration file
# 4 - the location of the new configuration file to be written
# 5 - the input fasta file on which the analysis should be done
# 6 - the number of CPUs available to the user
# 7 - the BUSCO database to be used

SAMPLENAME=$1
OUTPUTDIR=$2
STATIC_CONFIG=$3
CONFIG_LOC=$4
INPUT_FASTA=$5
CPUS=$6
BUSCO_DB=$7

URL="https://busco-data.ezlab.org/v4/data"

mkdir -p busco_downloads
wget -nd -r --no-parent -A $BUSCO_DB.*.tar.gz $URL/lineages/
tar -xvzf $BUSCO_DB.*.tar.gz
rm $BUSCO_DB.*.tar.gz
mv $BUSCO_DB busco_downloads/lineages/$BUSCO_DB