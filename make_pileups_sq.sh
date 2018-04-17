#!/bin/bash
#SBATCH --mem 10000

usage="Usage: $0 bamfile bedfile"

if [[ ! -f $1 ]] || [[ ! -f $2 ]] || [[ $1 -eq "-h" ]] || [[ $1 -eq "--help" ]]
then
    echo -e $usage
    exit 1
fi

dir=$(dirname $1)
name=$(basename $1 | sed "s/.bam//g")
samtools mpileup -f $HUMAN_GENOME -l $2 -q 10 -Q 20 $1 | gzip > $dir/$name.pileup.gz
