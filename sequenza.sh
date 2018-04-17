#!/bin/bash
#SBATCH -c 8

usage="Usage: $0 seqzfile ncores"
scriptdir=/home/dmalload/ppipilot/scripts

if [[ ! -f $1 ]]
then
    echo $usage
    exit 1
fi


Rscript $scriptdir/sequenza.R $1 $SLURM_CPUS_PER_TASK
