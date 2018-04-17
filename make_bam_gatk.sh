#!/bin/bash
#SBATCH -c 8

module load picard/2.9.2
set -o pipefail

#@K00274:115:HMJNCBBXX:6:1101:2960:1209 1:N:0:ATTACTCG+AGGCTATA
#@machine:run:flowcell:lane:tile:posx:posy[:UMI] read:isfiltered:controlnumber(0 almost always):index

##Temp dir
sample=$1
TMP_DIR=/scratch/dmalload/${sample}
mkdir $TMP_DIR

cd $sample

#Makes unaligned bam with readgroup
if [[ ! -f ${sample}_unaligned.bam ]]
then
    pu=$(gunzip -c *_R1_*fastq.gz | head -1 | sed "s/^[^:]*:[^:]*:\([^:]*\):\([^:]*\).*/\1.\2/g")
    pu="$pu.$sample"
    f1=$(ls *_R1_*fastq.gz)
    f2=$(ls *_R2_*fastq.gz)
    java -Xms512m -Xmx8G -jar $PICARDJAR FastqToSam F1=$f1 F2=$f2 O=${sample}_unaligned.bam RG=${sample}_1 SM=$sample PU=$pu PL=illumina CN=uiowa
    if [[ -f ${sample}_unaligned.bam ]]
    then
        rm -rf *fastq.gz ##we have backups in the nas
    fi
fi

#Marks illumina adapters in a new bam. This is not standard, and does tricks like artificially lowering the quality scores of bases that are marked as adapters. This bam is an intermediate file that will be deleted

if [[ ! -f ${sample}_markadapters.bam ]]
then
    java -Xms512m -Xmx8G -jar $PICARDJAR MarkIlluminaAdapters I=${sample}_unaligned.bam O=${sample}_markadapters.bam M=${sample}_markadapters.metrics TMP_DIR=$TMP_DIR #TMP_DIR optional to process large files
fi

#Piped: generates modified fastq from the previous bam, aligns them to the reference genome and adds back the lost information from the original bam file (also modifies some alignments, making primary those that are not chimeric if there are several options). The resulting bam is standard

if [[ ! -f ${sample}.bam ]]
then
    java -Xmx8G -jar $PICARDJAR SamToFastq \
    I=${sample}_markadapters.bam \
    FASTQ=/dev/stdout \
    CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
    TMP_DIR=$TMP_DIR | \
    bwa mem -M -t $SLURM_CPUS_PER_TASK -p $HUMAN_GENOME /dev/stdin | \
    java -Xmx16G -jar $PICARDJAR MergeBamAlignment \
    ALIGNED_BAM=/dev/stdin \
    UNMAPPED_BAM=${sample}_unaligned.bam \
    OUTPUT=${sample}.bam \
    R=$HUMAN_GENOME CREATE_INDEX=true ADD_MATE_CIGAR=true \
    CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
    INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
    TMP_DIR=$TMP_DIR 
fi

#Markduplicates
#Not piped since apparently it generates memory problems if done
if [[ ! -f ${sample}_mdups.bam ]]
then
    java -Xmx8G -jar $PICARDJAR MarkDuplicates I=${sample}.bam O=${sample}_mdups.bam ASSUME_SORTED=true M=${sample}_markduplicates.metrics CREATE_INDEX=true TMP_DIR=$TMP_DIR
fi

if [[ -f ${sample}_mdups.bam ]]
then
    rm -f ${sample}_unaligned.bam
    rm -f ${sample}_markadapters.bam
#    rm -f ${sample}.bam ##Just in case. #TODO: remove
fi

rm -rf $TMP_DIR
