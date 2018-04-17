#!/bin/bash
module load gatk/4.0.1.2

usage="$0 dir\noutputdir must have subdirectories, each with a normal vcf.gz file\n"

if [[ $# -ne 1 ]]
then
    echo -e $usage
    exit
fi

cd $1

vcfs=""

for i in */*.vcf
do
    vcfs=$vcfs" -vcfs "$i
done

eval "gatk CreateSomaticPanelOfNormals$vcfs -O pon.vcf.gz" #this is supposed to be insecure
