#!/bin/bash
#SBATCH --mem-per-cpu 16000M

usage="$0 bamfile intervals\n This script generates an outputfile adding .realigned.bam\n Intervals: intervals generated with RealignerTargetCreator"

if [[ $# -ne 2 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ ! -s $1 ]] || [[ ! -s $2 ]]
then
    echo -e $usage
    exit 1
fi

name=$(basename $1 | sed "s/.bam//g")
dir=$(dirname $1)

gatk -T IndelRealigner -R $HUMAN_GENOME -I $1 -known $(dirname $HUMAN_GENOME)/1000G_phase1.indels.*.vcf.gz -known $(dirname $HUMAN_GENOME)/Mills_and_1000G_gold_standard.indels.*.vcf.gz -targetIntervals $2 -o $dir/$name.realigned.bam
