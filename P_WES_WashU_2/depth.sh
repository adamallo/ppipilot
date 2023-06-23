#!/bin/bash

usage="$0 bamfile outputfile [bedfile]"

if [[ $# -lt 2 ]] || [[ $# -gt 3 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ ! -s $1 ]]
then
    echo -e $usage
    exit 1
fi

if [[ $# -eq 3 ]] && [[ ! -s $3 ]]
then
    echo "Error opening the bedfile $3"
    echo -e $usage
    exit 1
fi

echo "Filtering at QUAL >29 and MAPQ>49"

if [[ $# -eq 3 ]]
then
    echo "Calculating the depth of positions selected by the BED file $bedfile"
    samtools depth -q 29 -Q 49 $1 -b $3 | gzip > $2
else
    samtools depth -q 29 -Q 49 $1 | gzip > $2
fi
