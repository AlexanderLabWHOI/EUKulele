#!/bin/bash

ALLEXITS=0
# still need to check whether software already installed. XX
# still need to remove installs if they're not needed and figure out how the frick BUSCO will work.

mkdir -p references_bins/
export PATH=$PATH:references_bins/ >> ~/.bashrc

# INSTALL DIAMOND

diamond --version
if [ $? -ne 0 ]; then
    wget http://github.com/bbuchfink/diamond/releases/download/v0.9.36/diamond-linux64.tar.gz
    tar xzf diamond-linux64.tar.gz
    # instead make last command be diamond --help
    #export PATH=$PATH:.
    mv -f diamond references_bins/
    rm -rf diamond-linux64.tar.gz*
    diamond --version
fi
ALLEXITS=$(($ALLEXITS + $?))

# INSTALL BLAST
blastp -help
if [ $? -ne 0 ]; then
    wget -r --no-parent -A 'ncbi-blast-*+-x64-linux.tar.gz' ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
    tar -zxvf ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-*+-x64-linux.tar.gz
    shopt -s extglob
    BLASTDIR=ncbi-blast-*/bin
    cp ncbi-blast-*/bin/blastp references_bins/
    mv -f ncbi-blast-* references_bins/
    rm -rf ftp.ncbi.nlm.nih.gov
    #export PATH=ncbi-blast-2.10.1+/bin:$PATH
    #export PATH=$BLASTDIR:$PATH # I still have no idea why this won't save the wildcard expansion
    blastp -help
fi
ALLEXITS=$(($ALLEXITS + $?))

# INSTALL BUSCO=
busco --version
if [ $? -ne 0 ]; then
    rm -rf busco
    git clone https://gitlab.com/ezlab/busco.git
    cp busco/config/config.ini EUKulele/static/config.ini
    mv -vn busco references_bins/
    alias busco="references_bins/busco/bin/busco"
    export PATH=$PATH:references_bins/busco >> ~/.bashrc
    rm -rf busco
    busco --version
fi
#ALLEXITS=$(($ALLEXITS + $?))

# INSTALL TRANSDECODER
TransDecoder.Predict --version
if [ $? -ne 0 ]; then
    wget https://github.com/TransDecoder/TransDecoder/archive/TransDecoder-v5.5.0.tar.gz
    tar -zxvf TransDecoder-v5.5.0.tar.gz
    mv -vn TransDecoder-TransDecoder-v5.5.0 references_bins/TransDecoder
    export PATH=$PATH:references_bins/TransDecoder >> ~/.bashrc
    rm -rf *TransDecoder-v5.5.0*
    TransDecoder.Predict --version
fi
ALLEXITS=$(($ALLEXITS + $?))

echo $ALLEXITS

exit $ALLEXITS
