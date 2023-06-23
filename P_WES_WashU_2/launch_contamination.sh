#!/bin/bash
scriptdir=/home/dmalload/ppipilot/scripts
usage="$0 cases.txt outputdir \ncases.txt: space-separated file with casename file_T file_N\n"

if [[ $# -ne 2 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ ! -f $1 ]]
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
    normal=$normal"_pileup.table"
    tumor=$(readlink -e $tumor)
    tumor=$tumor"_pileup.table"
    mkdir -p $outputdir/$name
    dir=$(readlink -e $outputdir/$name)
    echo "Submitting case $name. Results in $dir"
    submit $scriptdir/make_contamination.sh $tumor $normal $name $dir
done < $1
