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
    limits="-l "$(readlink -e $4)
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

lofreq somatic -f $HUMAN_GENOME -n $normal -t $tumor --dbsnp $(dirname $HUMAN_GENOME)/dbsnp_138.*.nosomatic.vcf.gz  -o $output --threads $ncpu ${limits} --germline
