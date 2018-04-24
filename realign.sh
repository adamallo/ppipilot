#!/bin/bash
#SBATCH -c 8

module load picard/2.9.2

usage="$0 bam"

if [[ $# -ne 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-h" ]]
then
    echo -e $usage
    exit 1
fi

sample=$(echo $1 | sed "s/.bam//g")

TMP_DIR=/scratch/dmalload/$(basename $1 | sed "s/.bam//g")
mkdir $TMP_DIR

java -Xmx8G -jar $PICARDJAR RevertSam I=${sample}.bam O=${sample}_unaligned.bam RESTORE_ORIGINAL_QUALITIES=false REMOVE_DUPLICATE_INFORMATION=false

java -Xmx8G -jar $PICARDJAR SamToFastq \
    I=${sample}_unaligned.bam \
    FASTQ=/dev/stdout \
    CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
    TMP_DIR=$TMP_DIR | \
    bwa mem -M -t $SLURM_CPUS_PER_TASK -p $HUMAN_GENOME /dev/stdin | \
    java -Xmx16G -jar $PICARDJAR MergeBamAlignment \
    ALIGNED_BAM=/dev/stdin \
    UNMAPPED_BAM=${sample}_unaligned.bam \
    OUTPUT=${sample}_new.bam \
    R=$HUMAN_GENOME CREATE_INDEX=true ADD_MATE_CIGAR=true \
    CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
    INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
    TMP_DIR=$TMP_DIR EXPECTED_ORIENTATIONS=FR SO=coordinate \
    VALIDATION_STRINGENCY=SILENT
