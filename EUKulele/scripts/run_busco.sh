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

mkdir -p $OUTPUTDIR
        
cp euk-env/config/config.ini $CONFIG_LOC
sed -i '/out = /c\out = '$SAMPLENAME $CONFIG_LOC # the name of the output files
sed -i '/out_path = /c\out_path = '$OUTPUTDIR $CONFIG_LOC # what directory the output will be stored in
sed -i '/download_path = /c\download_path = ./busco_downloads/' $CONFIG_LOC
busco -i $INPUT_FASTA -l $BUSCO_DB -m proteins --cpu $CPUS --config $CONFIG_LOC -o $SAMPLENAME -f --offline
mv $OUTPUTDIR/$SAMPLENAME/*/* $OUTPUTDIR/$SAMPLENAME