#!/bin/bash

ALLEXITS=0

DATABASE=$1
REF_FASTA="reference.pep.fa"
REF_TABLE="taxonomy-table.txt"
REF_FASTA_URL=$2
REF_TABLE_URL=$3

mkdir -p ${PWD}/$DATABASE

if [[ $DATABASE == "mmetsp" ]]; then
    # Download MMETSP reference FASTA
    wget -O ${PWD}/$DATABASE/$REF_FASTA $REF_FASTA_URL
    ALLEXITS=$(($ALLEXITS + $?))
    
    # Download MMETSP reference taxonomy table
    wget -O ${PWD}/$DATABASE/$REF_TABLE $REF_TABLE_URL
    ALLEXITS=$(($ALLEXITS + $?))
    
    echo "All reference files for MMETSP downloaded to ${PWD}/$DATABASE"
elif [[ $DATABASE == "eukprot" ]]; then
    # Download tar of all EukProt files 
    wget -O ${PWD}/$DATABASE.tgz $REF_FASTA_URL
    ALLEXITS=$(($ALLEXITS + $?))
    
    # Unzip to proteins folder
    tar zxvf ${PWD}/$DATABASE/$DATABASE.tgz -C ${PWD}/$DATABASE
    ALLEXITS=$(($ALLEXITS + $?))
    
    # Download EukProt taxonomy file
    wget -O ${PWD}/$DATABASE/$REF_TABLE $REF_TABLE_URL
    ALLEXITS=$(($ALLEXITS + $?))
    
    ALLFILES=""
    for entry in "${PWD}/$DATABASE/$DATABASE/proteins"/*
    do
      ALLFILES=$ALLFILES" "$entry
    done
    
    for currfile in $ALLFILES
    do 
        ((cat $currfile | sed 's/\./N/g'); echo; echo) >> ${PWD}/$DATABASE/$REF_FASTA
    done
    ALLEXITS=$(($ALLEXITS + $?))
    
    echo "All reference files for EukProt downloaded to ${PWD}/$DATABASE"
elif [[ $DATABASE == "phylodb" ]]; then
    # Download PhyloDB files from Google Drive
    wget --load-cookies /tmp/cookies.txt https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id='$REF_FASTA_URL -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=$REF_FASTA_URL -O ${PWD}/$DATABASE/$DATABASE.tgz && rm -rf /tmp/cookies.txt
    gunzip -C ${PWD}/$DATABASE/$DATABASE.tgz > ${PWD}/$DATABASE/$REF_FASTA
    
    # Download PhyloDB taxonomy table
    wget --load-cookies /tmp/cookies.txt https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id='$REF_FASTA_URL -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=$REF_FASTA_URL -O ${PWD}/$DATABASE/$DATABASE.table.tgz && rm -rf /tmp/cookies.txt
    gunzip -C ${PWD}/$DATABASE/$DATABASE.table.tgz > ${PWD}/$DATABASE/$REF_TABLE_URL
    
    echo "All reference files for PhyloDB downloaded to ${PWD}/$DATABASE"
else
    exit 1
fi