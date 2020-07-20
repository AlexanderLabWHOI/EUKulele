#!/bin/bash

# One argument is the directory in which the software will be stored.

ALLEXITS=0
# still need to check whether software already installed. XX
# still need to remove installs if they're not needed and figure out how the frick BUSCO will work.

DEST_DIR=$1
echo $DEST_DIR
mkdir -p $DEST_DIR #references_bins/
export PATH=$PATH:$DEST_DIR
echo "export PATH=$PATH:$DEST_DIR" >> ~/.bashrc

# INSTALL DIAMOND

diamond --version
if [ $? -ne 0 ]; then
    wget http://github.com/bbuchfink/diamond/releases/download/v0.9.36/diamond-linux64.tar.gz
    tar xzf diamond-linux64.tar.gz
    mv -f diamond $DEST_DIR
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
    cp ncbi-blast-*/bin/blastp $DEST_DIR
    mv -f ncbi-blast-* $DEST_DIR
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
    alias busco=""$DEST_DIR"busco/bin/busco"
    export PATH=$PATH:"$DEST_DIR"busco
    export PATH=$PATH:"$DEST_DIR"busco >> ~/.bashrc
    rm -rf busco
    rm -rf busco*.log
    busco --version
fi
ALLEXITS=$(($ALLEXITS + $?))

# INSTALL TRANSDECODER
TransDecoder.Predict --version
if [ $? -ne 0 ]; then
    wget https://github.com/TransDecoder/TransDecoder/archive/TransDecoder-v5.5.0.tar.gz
    tar -zxvf TransDecoder-v5.5.0.tar.gz
    mv -vn TransDecoder-TransDecoder-v5.5.0 "$DEST_DIR"TransDecoder
    export PATH=$PATH:"$DEST_DIR"TransDecoder >> ~/.bashrc
    rm -rf *TransDecoder-v5.5.0*
    TransDecoder.Predict --version
fi
ALLEXITS=$(($ALLEXITS + $?))

if [[ "$0" != "$BASH_SOURCE" ]] ; then
  # this script is executed via `source`!
  # An `exit` will close the user's console!
  EXIT=return
else
  # this script is not `source`-d, it's safe to exit via `exit`
  EXIT=exit
fi

echo $ALLEXITS

$EXIT $ALLEXITS