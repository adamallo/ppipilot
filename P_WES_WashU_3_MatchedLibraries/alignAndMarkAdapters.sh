#!/bin/bash
#SBATCH -c 8

scriptdir=/home/dmalload/ppistudy/scripts

module load picard/2.9.2
set -o pipefail

usage="$0 bamfile umi \n umi: binary flag, 1 includes a 8bp UMI, 0 no UMI "

if [[ $# -lt 2 ]] || [[ ! -f $1 ]]
then
    echo -e $usage
    exit 1
fi

umi=1
umi=$2
read_structure="8M+T 8M+T" ##Only change for a different data origin

bamfile=$(readlink -e $1)
runname=$(samtools view -H $1 | grep "@RG" | sed 's/.*LB:\([^\t]*\).*/\1/')
samplename=$(samtools view -H $1 | grep "@RG" | sed 's/.*SM:\([^\t]*\).*/\1/')

if [[ "$SLURM_CPUS_PER_TASK" == "" ]]
then
    SLURM_CPUS_PER_TASK=1
fi

##Temp dir
TMP_DIR=/scratch/dmalload/${runname}
mkdir -p $TMP_DIR

dir=$(dirname $1)
mkdir -p $dir/$samplename 
cd $dir/$samplename

##Settings with a highest probability of being changed
UMIminMAPQ=50 ###Filtering for UMI calling (In the IDT manual they use 10, the default is 30, looking at the empirical distribution of our data basically all of them are at 60). Here I am more stringent since this info cannot be recovered from the consensus and may influence strongly their estimation
UMIminPHRED=20 ###Filtering for UMI calling (In the IDT manual they use 30, the default is 10). This info cannot be recovered from consensus reads. However, it is less important than UMIminMAPQ since even if we include bases with low quality they should be "consensed out" if we have enough number of reads comming from the same molecule and these phreds are taking into consideration to calculate the phred score of the consensus position.
UMIminreads=1 ###I can filter them a posteriori, doing it here just improves computation time, but I do not think it is worth losing data considering low coverage in some samples

##############

if [[ ! -s ${runname}_markadapters.bam ]] && [[ ! -s ${runname}_unaligned.bam ]]
then
    if [[ $umi -eq 1 ]]
    then
        fgbio.sh --tmp-dir=$TMP_DIR --compression=0 ExtractUmisFromBam --input=$bamfile --output=/dev/stdout --read-structure=$read_structure --molecular-index-tags=ZA ZB --single-tag=RX | sambamba view -f bam --nthreads=$SLURM_CPUS_PER_TASK -o ${runname}_unaligned.bam /dev/stdin
        if [[ -s ${runname}_unaligned.bam ]]
        then
            rm -f $bamfile
        fi
    else    
        mv $bamfile ${runname}_unaligned.bam
    fi
fi

#Marks illumina adapters in a new bam. This is not standard, and does tricks like artificially lowering the quality scores of bases that are marked as adapters. This bam is an intermediate file that will be deleted

if [[ ! -s ${runname}.bam ]] && [[ ! -s ${runname}_markadapters.bam ]]
then
    mkdir -p qc_$runname/markadapters
    java -Xms512m -Xmx8G -jar $PICARDJAR MarkIlluminaAdapters I=${runname}_unaligned.bam O=${runname}_markadapters.bam M=qc_$runname/markadapters/${runname}_markadapters.metrics TMP_DIR=$TMP_DIR #TMP_DIR optional to process large files
fi

#Piped: generates modified fastq from the previous bam, aligns them to the reference genome and adds back the lost information from the original bam file (also modifies some alignments, making primary those that are not chimeric if there are several options). The resulting bam is standard

if [[ ! -s ${runname}.bam ]]
then
    java -Xmx8G -jar $PICARDJAR SamToFastq \
    I=${runname}_markadapters.bam \
    FASTQ=/dev/stdout \
    CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
    TMP_DIR=$TMP_DIR VALIDATION_STRINGENCY=SILENT | \
    bwa mem -M -t $SLURM_CPUS_PER_TASK -p $HUMAN_GENOME /dev/stdin | \
    java -Xmx16G -jar $PICARDJAR MergeBamAlignment \
    ALIGNED_BAM=/dev/stdin \
    UNMAPPED_BAM=${runname}_unaligned.bam \
    OUTPUT=${runname}.bam \
    R=$HUMAN_GENOME CREATE_INDEX=true ADD_MATE_CIGAR=true \
    CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
    INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
    TMP_DIR=$TMP_DIR EXPECTED_ORIENTATIONS=FR SO=queryname \
    VALIDATION_STRINGENCY=SILENT

    if [[ -s ${runname}.bam ]]
    then
        rm -f ${runname}_unaligned.bam ${runname}_markadapters.bam
    fi
fi
