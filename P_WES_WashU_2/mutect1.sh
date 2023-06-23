#!/bin/bash
#SBATCH -c 8

##Configuration variables

usage="$0 tumor_bam normal_bam output_dir [bed]"

if [[ $# -lt 3 ]] || [[ ! -d $3 ]] || [[ ! -f $2 ]] || [[ ! -f $1 ]]
then
    echo -e $usage
    exit
fi

limits=""
if [[ $# -gt 3 ]] && [[ $4 != "" ]]
then
    limits="--intervals "$(readlink -e $4)
fi

if [[ "$SLURM_CPUS_PER_TASK" == "" ]]
then
    ncpu=1
else
    ncpu=$SLURM_CPUS_PER_TASK
fi


normal=$(readlink -e $2)
tumor=$(readlink -e $1)
output=$(basename $3)
outputdir=$(readlink -e $3)
cd $3

tumorname=$(samtools view -H $tumor | sed -n '/@RG/p' | sed "s/.*SM:\([^ \t]*\).*/\1/g")

mutect --reference_sequence $HUMAN_GENOME --tumor_sample_name tumor --normal_sample_name normal --bam_tumor_sample_name $tumorname --cosmic $(dirname $HUMAN_GENOME)/cosmic.*.vcf --dbsnp $(dirname $HUMAN_GENOME)/dbsnp_138.*.vcf --input_file:normal $normal --input_file:tumor $tumor --out $output.out --vcf $output.vcf --coverage_file $output.wig.txt -nt $ncpu ${limits}
