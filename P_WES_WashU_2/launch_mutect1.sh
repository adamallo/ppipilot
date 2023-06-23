#!/bin/bash
scriptdir=/home/dmalload/ppipilot/scripts
usage="$0 cases.txt outputdir [bed]\ncases.txt: space-separated file with casename file_T file_N.\n"

if [[ $# -lt 2 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ ! -f $1 ]]
then
    echo -e $usage
    exit
fi

outputdir=$2

echo "Outputdir $outputdir"
mkdir -p $outputdir

bed=""

if [[ $# -gt 2 ]] && [[ $3 != "" ]]
then
    bed=$3
fi

while read name tumor normal
do
    normal=$(readlink -e $normal)
    tumor=$(readlink -e $tumor)
    mkdir -p $outputdir/$name
    dir=$(readlink -e $outputdir/$name)
    echo "Submitting case $name. Results in $dir"
    submit $scriptdir/mutect1.sh $tumor $normal $dir $bed
done < $1
