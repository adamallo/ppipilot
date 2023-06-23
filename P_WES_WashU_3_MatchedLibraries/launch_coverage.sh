#!/bin/bash
scriptdir=/home/dmalload/ppistudy/scripts
usage="$0 dir \ndir: directory with directories with directory_mdups.bam files inside\n"

if [[ $# -ne 1 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ ! -d $1 ]]
then
    echo -e $usage
    exit
fi

for file in $1/*/*.bam.bai
do
    name=$(dirname $file)
    echo "Submitting case $name"
    submit $scriptdir/gatk_coverage_exome.sbatch $(echo $file | sed "s/\.bai$//")
    #submit $scriptdir/gatk_coverage.sbatch $file
done
