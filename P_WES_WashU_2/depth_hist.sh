#!/bin/bash
module load bedtools2/2.24.0

usage="$0 bamfile bedfile nodups"

if [[ $# -ne 3 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ ! -s $1 ]] || [[ ! -s $2 ]]
then
    echo -e $usage
    exit 1
fi

dir=$(dirname  $1)
name=$(basename $1 | sed "s/.bam//g")
bed=$2
bam=$1
nodups=$3

filter=""
dupname="dup"

if [[ $nodups -eq 1 ]]
then
    filter="-F 1024 "
    dupname="nodup"
fi

mkdir -p $dir/qc_$name

samtools view -u -L $bed ${filter} $bam | bedtools coverage -hist -sorted -b stdin -a $bed | grep ^all > $dir/qc_$name/$name.$dupname.depthHist.tsv
