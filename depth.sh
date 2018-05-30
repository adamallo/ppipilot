#!/bin/bash

usage="$0 bamfile outputfile"

if [[ $# -ne 2 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ ! -s $1 ]]
then
    echo -e $usage
    exit 1
fi

echo "Filtering at QUAL >29 and MAPQ>49"

samtools depth -q 29 -Q 49 $1 | gzip > $2
