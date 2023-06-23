#!/bin/bash

module load fastqc/0.11.3

usage="$0 bamfile"

if [[ $# -ne 1 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ ! -f $1 ]]
then
    echo -e $usage
    exit 1
fi

dir=$(dirname  $1)
name=$(basename $1 | sed "s/.bam//g")

mkdir -p $dir/qc_$name/fastqc

fastqc $1 -o $dir/qc_$name/fastqc
