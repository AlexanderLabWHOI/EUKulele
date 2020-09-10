#!/bin/bash

SAMPLE_NAME=$1
FASTA_OUT=$2
OUT_DIR=$3

SEQUENCE_FILES=$(ls -d $OUT_DIR/busco/$SAMPLE_NAME/run_eukaryota_odb10/busco_sequences/*/* || true)
SEQUENCE_FILES=$(echo $SEQUENCE_FILES | tr "\n" " ")
len=${#SEQUENCE_FILES[@]}
#if [ "$SEQUENCE_FILES" != ""]; then
if [ "$SEQUENCE_FILES" != " " ] && [ "$SEQUENCE_FILES" != "" ]; then
    cat $SEQUENCE_FILES > $FASTA_OUT
fi