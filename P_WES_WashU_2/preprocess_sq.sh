#!/bin/bash
#SBATCH --mem 10000

module unload python
module load pypy27/5.6

scriptdir=/home/dmalload/ppipilot/scripts

usage="Usage: $0 tumor normal normal2 outname\n Normal2: is only used to estimate depth and depth.ratio. It can be a normal form other patient from the same batch in order to reduce batch effects in the exome capture procedure"

if [[ $# -ne 4 ]] || [[ ! -f $1 ]] || [[ ! -f $2 ]] || [[ ! -f $3 ]]
then
    echo -e $usage
    exit 1
fi

mkdir -p $4

#gcfile=$(echo $(dirname $HUMAN_GENOME)/*gc*Base.txt.gz | xargs ls )
#gcfile=$(echo $(dirname $HUMAN_GENOME)/*gc*Base.lite.oorder.txt.gz | xargs ls )

##TODO: delete this
gcfile="/home/dmalload/my_storage/hg19/hg19.gc50Base.lite.oorder.txt.gz"
HUMAN_GENOME="/home/dmalload/my_storage/hg19/hg19.fa"

dir=$4
outname=$(basename $dir)

#pypy $scriptdir/sequenza-utils.py pileup2seqz -gc $gcfile -t $1 -n $2 | gzip > $dir/${outname}_temp.seqz.gz
#pypy $scriptdir/sequenza-utils.py seqz-binning -w 50 -s $dir/${outname}_temp.seqz.gz | gzip > $dir/$outname.seqz.gz

##I am filtering by mapq and q at the pileup level
##I may need to modify the default filter -N 20 that filters positions with less than 20 reads in total (the two samples together) when using consensus reads

pypy `which sequenza-utils` bam2seqz -p -f illumina -gc $gcfile -q 0 --fasta $HUMAN_GENOME -t $1 -n $2 -n2 $3 -N 2 -o $dir/${outname}_temp.seqz.gz
pypy `which sequenza-utils` seqz_binning -w 50 -s $dir/${outname}_temp.seqz.gz -o $dir/$outname.seqz.gz

rm -f $dir/${outname}_temp.seqz.gz
