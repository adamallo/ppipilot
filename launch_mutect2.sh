#!/bin/bash
scriptdir=/home/dmalload/ppipilot/scripts
usage="$0 cases.txt outputdir pon \ncases.txt: space-separated file with casename file_T file_N. PON: vcf of the panel of normals for this experiment\n"

if [[ $# -ne 3 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ ! -f $1 ]] || [[ ! -f $3 ]]
then
    echo -e $usage
    exit
fi

outputdir=$2

echo "Outputdir $outputdir"
mkdir -p $outputdir

while read name tumor normal
do
    normal=$(readlink -e $normal)
    tumor=$(readlink -e $tumor)
    mkdir -p $outputdir/$name
    dir=$(readlink -e $outputdir/$name)
    echo "Submitting case $name. Results in $dir"
    submit $scriptdir/mutect2.sh $tumor $normal $3 $dir
done < $1
