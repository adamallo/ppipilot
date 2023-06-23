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
    limits=" --callRegions="$(readlink -e $4)
fi

normal=$(readlink -e $2)
tumor=$(readlink -e $1)
output=$(basename $3)
outputdir=$(readlink -e $3)

strelkaS --normalBam=$normal --tumorBam=$tumor --referenceFasta=$HUMAN_GENOME --exome --outputCallableRegions --runDir=$outputdir${limits}

$outputdir/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK --memGb $(($SLURM_MEM_PER_CPU*$SLURM_CPUS_PER_TASK/1024))
