#!/bin/bash
#SBATCH -c 8

usage="Usage: $0 seqzfile chr female\nchr: string to append to the chr number to identify chromosomes (i.e., hg19 vs b37)\nfemale: TRUE if female, FALSE if male"
scriptdir=/home/dmalload/ppipilot/scripts

if [[ $# -ne 3 ]] || [[ ! -f $1 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
    echo -e $usage
    exit 1
fi


Rscript $scriptdir/sequenza.R $1 $2 $3 $SLURM_CPUS_PER_TASK
