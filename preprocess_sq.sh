#!/bin/bash
#SBATCH --mem 10000

module unload python
module load pypy27/5.6

scriptdir=/home/dmalload/ppipilot/scripts

usage="Usage: $0 tumor normal outname"

if [[ $# -ne 3 ]] || [[ ! -f $1 ]] || [[ ! -f $2 ]]
then
    echo -e $usage
    exit 1
fi

mkdir -p $3

gcfile="$(dirname $HUMAN_GENOME)/b37.gc50Base.txt.gz"
dir=$3
outname=$(basename $dir)

pypy $scriptdir/sequenza-utils.py pileup2seqz -gc $gcfile -t $1 -n $2 | gzip > $dir/${outname}_temp.seqz.gz
#pypy `which sequenza-utils.py` bam2seqz -gc $gcfile --fasta $HUMAN_GENOME -t $1 -n $2 | gzip > $dir/$3_temp.seqz.gz
pypy $scriptdir/sequenza-utils.py seqz-binning -w 50 -s $dir/${outname}_temp.seqz.gz | gzip > $dir/$outname.seqz.gz
rm -f $dir/${outname}_temp.seqz.gz
