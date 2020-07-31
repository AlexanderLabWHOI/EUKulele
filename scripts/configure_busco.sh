#!/bin/bash

# A script that takes 1 positional argument to configure BUSCO for a dataset:
# 1 - the BUSCO database to be used

BUSCO_DB=$1

URL="https://busco-data.ezlab.org/v4/data"

mkdir -p busco_downloads/lineages
if [ ! -d busco_downloads/lineages/"$BUSCO_DB" ]; then
    wget -nd -r --no-parent -A "$BUSCO_DB".*.tar.gz $URL/lineages/ > /dev/null 2>&1
    if [ $? -eq 0 ]; then
        tar -xvzf "$BUSCO_DB".*.tar.gz* || true
        rm -f "$BUSCO_DB".*.tar.gz* || true
        if [ ! -d busco_downloads/lineages/"$BUSCO_DB" ]; then
            mv -f "$BUSCO_DB" busco_downloads/lineages/"$BUSCO_DB" || true
        fi
    else
        echo "Download unsuccessful"
    fi
fi

if [ -d busco_downloads/lineages/"$BUSCO_DB" ]; then
    exit 0;
else
    exit 1;
fi