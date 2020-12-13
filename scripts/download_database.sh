#!/bin/bash

ALLEXITS=0

DATABASE=$1
REF_FASTA="reference.pep.fa"
REF_TABLE="taxonomy-table.txt"
REF_FASTA_URL=$2
REF_TABLE_URL=$3
REFERENCE_DIR=$4

#mkdir -p ${PWD}/$REFERENCE_DIR/$DATABASE
mkdir -p $REFERENCE_DIR/$DATABASE

if [[ $DATABASE == "marmmetsp" ]]; then
    # Download MMETSP reference FASTA
    wget -O $REFERENCE_DIR/$DATABASE/$REF_FASTA $REF_FASTA_URL
    ALLEXITS=$(($ALLEXITS + $?))
    
    # Download MMETSP reference taxonomy table
    wget -O $REFERENCE_DIR/$DATABASE/$REF_TABLE $REF_TABLE_URL
    ALLEXITS=$(($ALLEXITS + $?))
    
    echo "All reference files for MarRef-MMETSP downloaded to $REFERENCE_DIR/$DATABASE"
elif [[ $DATABASE == "mmetsp" ]]; then
    # Download MMETSP reference FASTA
    wget -O $REFERENCE_DIR/$DATABASE/$REF_FASTA $REF_FASTA_URL
    ALLEXITS=$(($ALLEXITS + $?))
    
    # Download MMETSP reference taxonomy table
    wget -O $REFERENCE_DIR/$DATABASE/$REF_TABLE $REF_TABLE_URL
    ALLEXITS=$(($ALLEXITS + $?))
    
    echo "All reference files for MMETSP downloaded to $REFERENCE_DIR/$DATABASE"
elif [[ $DATABASE == "eukprot" ]]; then
    # Download tar of all EukProt files 
    wget -O $REFERENCE_DIR/$DATABASE/$DATABASE.tgz $REF_FASTA_URL
    ALLEXITS=$(($ALLEXITS + $?))
    
    # Unzip to proteins folder
    tar zxvf $REFERENCE_DIR/$DATABASE/$DATABASE.tgz -C $REFERENCE_DIR/$DATABASE
    ALLEXITS=$(($ALLEXITS + $?))
    
    # Download EukProt taxonomy file
    wget -O $REFERENCE_DIR/$DATABASE/$REF_TABLE $REF_TABLE_URL
    ALLEXITS=$(($ALLEXITS + $?))
    
    ALLFILES=""
    for entry in "$REFERENCE_DIR/$DATABASE/proteins"/*
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
    # Download PhyloDB reference FASTA
    wget -O $REFERENCE_DIR/$DATABASE/$REF_FASTA.gz $REF_FASTA_URL
    gunzip -f $REFERENCE_DIR/$DATABASE/$REF_FASTA.gz
    #sed -i -e 's/>* .*$//' $REFERENCE_DIR/$DATABASE/$REF_FASTA
    #sed -i $'s/\t/    /g' $REFERENCE_DIR/$DATABASE/$REF_FASTA
    ALLEXITS=$(($ALLEXITS + $?))
    
    # Download PhyloDB reference taxonomy table
    wget -O $REFERENCE_DIR/$DATABASE/$REF_TABLE.gz $REF_TABLE_URL
    gunzip -f $REFERENCE_DIR/$DATABASE/$REF_TABLE.gz
    ALLEXITS=$(($ALLEXITS + $?))
    
    # Download PhyloDB files from Google Drive
    #wget --load-cookies /tmp/cookies.txt https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id='$REF_FASTA_URL -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=$REF_FASTA_URL -O ${PWD}/$DATABASE/$DATABASE.tgz && rm -rf /tmp/cookies.txt
    #gunzip -c ${PWD}/$DATABASE/$DATABASE.tgz > ${PWD}/$DATABASE/$REF_FASTA
    
    # Download PhyloDB taxonomy table
    #wget --load-cookies /tmp/cookies.txt https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id='$REF_FASTA_URL -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=$REF_FASTA_URL -O ${PWD}/$DATABASE/$DATABASE.table.tgz && rm -rf /tmp/cookies.txt
    #gunzip -c ${PWD}/$DATABASE/$DATABASE.table.tgz > ${PWD}/$DATABASE/$REF_TABLE_URL
    
    echo "All reference files for PhyloDB downloaded to ${PWD}/$DATABASE"
elif [[ $DATABASE == "eukzoo" ]]; then
    zenodo_get 1476236
    mv EukZoo_taxonomy_table_v_0.2.tsv $REFERENCE_DIR/$DATABASE/$REF_TABLE
    mv EukZoo_v_0.2.faa $REFERENCE_DIR/$DATABASE/$REF_FASTA
    rm -f EukZoo_creation_and_cleanup.docx
    rm -f EukZoo_KEGG_annotation_v_0.2.tsv
else
    echo "Specified database not found."
    exit 1
fi