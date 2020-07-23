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

mkdir -p $OUTPUTDIR
       
configini="EUKulele/static/config.ini"

#rm -rf busco
#git clone https://gitlab.com/ezlab/busco.git
#cp busco/config/config.ini $configini

# euk-env/config/config.ini
#cp euk-env/config/config.ini $CONFIG_LOC
#cp $configini $CONFIG_LOC
BUSCO_DIR="$(dirname $(which busco))"
BUSCO_CONFIG_FILE="$BUSCO_DIR"/../config/config.ini
busco_configurator.py $BUSCO_CONFIG_FILE $CONFIG_LOC
echo "python3 busco_configurator.py $BUSCO_CONFIG_FILE $CONFIG_LOC"
sed -i '/out = /c\out = '$SAMPLENAME $CONFIG_LOC # the name of the output files
sed -i '/out_path = /c\out_path = '$OUTPUTDIR $CONFIG_LOC # what directory the output will be stored in
sed -i '/download_path = /c\download_path = ./busco_downloads/' $CONFIG_LOC
# ./references_bins/busco/bin/busco
busco -i $INPUT_FASTA -l $BUSCO_DB -m proteins --cpu $CPUS --config $CONFIG_LOC -o $SAMPLENAME -f --offline
mv $OUTPUTDIR/$SAMPLENAME/*/* $OUTPUTDIR/$SAMPLENAME