#!/bin/bash
scriptdir=/home/dmalload/ppipilot/scripts
usage="$0 dir [bed] \ndir: directory with directories with directory .bam files inside\n"

if [[ $# -lt 1 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ ! -d $1 ]]
then
    echo -e $usage
    exit
fi

bed=""

if [[ $# -gt 1 ]] && [[ $2 != "" ]]
then
    bed=$2
fi

for dir in $1/*
do
    file=$dir/$(basename $dir).bam
    if [[ ! -s $file ]]
    then
        echo "Skipping $dir"
    else
        name=$(dirname $file)
        echo "Submitting case $name"
        submit $scriptdir/make_pileupsummaries.sh $file $bed
    fi
done
