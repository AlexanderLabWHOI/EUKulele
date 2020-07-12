#!/bin/bash

# A script that takes 7 positional arguments to run BUSCO on a dataset:
# 0 - the descriptive name for the sample (1 word)
# 1 - the output directory that the results should be written to
# 2 - the location of the static/template configuration file
# 3 - the location of the new configuration file to be written
# 4 - the input fasta file on which the analysis should be done
# 5 - the number of CPUs available to the user
# 6 - the BUSCO database to be used

SAMPLENAME=$0
OUTPUTDIR=$1
STATIC_CONFIG=$2
CONFIG_LOC=$3
INPUT_FASTA=$4
CPUS=$5
BUSCO_DB=$6

mkdir -p $OUTPUTDIR
./busco_configurator.py $BUSCO_CONFIG_FILE {params.static_config}
cp $STATIC_CONFIG $CONFIG_LOC
sed -i '/out = /c\out = '$SAMPLENAME $CONFIG_LOC # the name of the output files
sed -i '/out_path = /c\out_path = '$OUTPUTDIR $CONFIG_LOC # what directory the output will be stored in
busco -i $INPUT_FASTA -l $BUSCO_DB -m proteins --cpu $CPUS --config $CONFIG_LOC -f
mv $OUTPUTDIR/$SAMPLENAME/*/* $OUTPUTDIR/$SAMPLENAME