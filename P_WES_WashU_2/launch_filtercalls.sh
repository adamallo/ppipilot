#!/bin/bash
scriptdir=/home/dmalload/ppipilot/scripts
usage="$0 cases.txt outputdir [bed] \ncases.txt: space-separated file with casename file_T file_N. Outputdir: directory with directories with the results of variant calling directory.vcf and contamination estimates directory.table\n"

if [[ $# -lt 2 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ ! -f $1 ]] || [[ ! -d $2 ]]
then
    echo -e $usage
    exit
fi

outputdir=$2

bed=""

if [[ $# -gt 2 ]] && [[ $3 != "" ]]
then
    bed=$3
fi

while read name tumor normal
do
    variants=$(readlink -e "$outputdir/$name/$name.vcf.gz")
    contamination=$(readlink -e "$outputdir/$name/${name}.table")
    echo "Submitting case $name. Results in $outputdir"
    submit $scriptdir/filtercalls.sh $variants $contamination "$outputdir/$name/${name}_filtered.vcf.gz" $bed
done < $1
