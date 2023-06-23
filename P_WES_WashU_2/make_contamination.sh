#!/bin/bash

module load gatk/4.0.1.2

usage="$0 tumor_pileuptable normal_pileuptable outname outdir"

if [[ $# -ne 4 ]] || [[ ! -f $1 ]] || [[ ! -f $2 ]]
then
    echo -e $usage
    exit
fi

mkdir -p $4


gatk CalculateContamination \
-I $1 \
-O $4/$3_nomatched.table

gatk CalculateContamination \
-I $1 \
-matched $2 \
-O $4/$3.table
