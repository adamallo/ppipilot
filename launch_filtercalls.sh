#!/bin/bash
scriptdir=/home/dmalload/ppipilot/scripts
bed=/home/dmalload/ppipilot/beds/xgen-exome-research-panel-targets.bed
usage="$0 cases.txt outputdir \ncases.txt: space-separated file with casename file_T file_N. Outputdir: directory with directories with the results of variant calling directory.vcf and contamination estimates directory.table\n"

if [[ $# -ne 2 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ ! -f $1 ]] || [[ ! -d $2 ]]
then
    echo -e $usage
    exit
fi

outputdir=$2

while read name tumor normal
do
    variants=$(readlink -e "$outputdir/$name/$name.vcf")
    contamination=$(readlink -e "$outputdir/$name/${name}_nomatched.table") #TODO: remove _nomatched
    echo "Submitting case $name. Results in $outputdir"
    submit $scriptdir/filtercalls.sh $variants $bed $contamination "$outputdir/$name/${name}_filtered.vcf" #TODO: use gz next time
done < $1
