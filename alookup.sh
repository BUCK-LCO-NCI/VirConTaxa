#!/bin/sh
LD_LIBRARY_PATH=/usr/local/MySQL/5.5.41/lib/
echo "$1 $2 $3"
#/data/BUCKLAB/SCRIPTS/alllookup -f $1 -r $2 $3
/data/starrettgj/scripts/RCAcontigSnakemake/alllookup -f $1 -r $2 $3
