#!/bin/bash

module load gatk/4.0.1.2

##Configuration variables

usage="$0 variants bed contamination_table output"

if [[ $# -ne 4 ]] || [[ ! -f $3 ]] || [[ ! -f $2 ]] || [[ ! -f $1 ]]
then
    echo -e $usage
    exit
fi

gatk FilterMutectCalls \
-V $1 \
-L $2 \
--contamination-table $3 \
-O $4
