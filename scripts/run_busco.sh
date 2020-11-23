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
BUSCO_MODE=$7 # proteins vs. transcriptome

mkdir -p $OUTPUTDIR
       
configini="EUKulele/static/config.ini"

BUSCO_DIR="$(dirname $(which busco))"
BUSCO_CONFIG_FILE="$BUSCO_DIR"/../config/config.ini
CONFIG_DIR="$(dirname $CONFIG_LOC)"
if [ ! -f "$BUSCO_CONFIG_FILE" ]; then
    BUSCO_CONFIG_FILE="$(find "$CONDA_PREFIX" | grep busco_config.ini)" #"$CONDA_PREFIX"/config/config.ini
fi
busco_configurator.py $BUSCO_CONFIG_FILE $CONFIG_LOC
echo "python3 busco_configurator.py $BUSCO_CONFIG_FILE $CONFIG_LOC"
sed -i '/out = /c\out = '$SAMPLENAME $CONFIG_LOC # the name of the output files
sed -i '/out_path = /c\out_path = '$OUTPUTDIR $CONFIG_LOC # what directory the output will be stored in
sed -i '/download_path = /c\download_path = ./busco_downloads/' $CONFIG_LOC
# ./references_bins/busco/bin/busco

mkdir -p $OUTPUTDIR/$SAMPLENAME
busco -i $INPUT_FASTA -l $BUSCO_DB -m $BUSCO_MODE --cpu $CPUS --config $CONFIG_LOC -o $SAMPLENAME -f --offline

# need to change back to putting everything back into the base directory to avoid conflicts w new DB names...?

#if [ -f "$OUTPUTDIR/$SAMPLENAME/*/*" ]; then
#    mv $OUTPUTDIR/$SAMPLENAME/*/* $OUTPUTDIR/$SAMPLENAME
#for dir in *; do if [ -d $dir ]; then echo $dir/*; fi; done