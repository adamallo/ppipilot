#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu 8000M

usage="$0 inputfile outputpref(dir) rmbam\brmbam: binary flag to indicate if the original bam should be deleted or not"

if [[ $# -ne 3 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help"]] || [[ ! -f $1 ]]
then
    echo -e $usage
    exit 1
fi

mkdir -p $(dirname $2)

fgbio.sh CollectDuplexSeqMetrics --input=$1 --output=$2 -u=T --description=$(basename $2)

if [[ $3 -eq 1 ]]
then
    rm -f $1
fi
