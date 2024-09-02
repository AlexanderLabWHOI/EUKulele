#!/bin/bash

awk '{ if (NR > 1 && index(prev, ">") == 0 && index($0, ">") == 0) printf "%s%s", prev, $0; else if (NR > 1) print prev; prev = $0 } END { print prev }' $1
