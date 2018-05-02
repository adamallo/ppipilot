#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu 8000M

usage="$0 inputfile outputpref(dir) description rmbam\brmbam: binary flag to indicate if the original bam should be deleted or not"

if [[ $# -ne 4 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ ! -f $1 ]]
then
    echo -e $usage
    exit 1
fi

mkdir -p $(dirname $2)

fgbio.sh CollectDuplexSeqMetrics --input=$1 --output=$2 -u=T --description=$3

if [[ $4 -eq 1 ]]
then
    rm -f $1
fi
