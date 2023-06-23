#!/bin/bash

module load gatk/4.0.1.2

##Configuration variables

usage="$0 variants contamination_table output [bed]"

if [[ $# -lt 3 ]] || [[ ! -f $1 ]] || [[ ! -f $2 ]]
then
    echo -e $usage
    exit
fi

limits=""
if [[ $# -gt 3 ]] && [[ $4 != "" ]]
then
    limits="-L "$(readlink -e $4)
fi

gatk FilterMutectCalls \
-V $1 \
--contamination-table $2 \
-O $3 ${limits}
