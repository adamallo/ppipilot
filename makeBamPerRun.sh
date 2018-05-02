#!/bin/bash
#SBATCH -c 8

module load picard/2.9.2
set -o pipefail

usage="$0 R1.bam R2.bam umi\n umi: binary flag, 1 includes a 8bp UMI, 0 no UMI "
umi=1

##Settings with a highest probability of being changed
read_structure="8M+T 8M+T" ##Only change for a different data origin
pl="illumina"
cn="washu"

#Example of fastq file naming format
#H_ZT-591-3-lib1_S14_L002_I1_001.fastq.gz
#H_([^_]*)_([^_]*)_([^_]*)_([^_]*)_([^_]*).fastq.gz \1: libraryname \2: sample \3: lane \4: IDreads [I1,I2,R1,R2] \5: SomeKindOfId

#Example of Fastq seq identifier and description with legend
#@K00274:115:HMJNCBBXX:6:1101:2960:1209 1:N:0:ATTACTCG+AGGCTATA
#@machine:run:flowcell:lane:tile:posx:posy[:UMI] read:isfiltered:controlnumber(0 almost always):index

#if [[ $# -ne 3 ]] || [[ ! -s $1 ]] || [[ ! -s $2 ]]
#then
#    echo -e $usage
#    exit 1
#fi

R1=$(readlink -e $1)
R2=$(readlink -e $2)
runname=$(basename $1 | sed "s/_.._[0-9]\+.fastq.gz$//")
samplelib=$(echo $runname | sed "s/\([^_]*\)_\([^_]*\)_\([^_]*\)_\([^_]*\)$/\2/")
samplename=$(echo $samplelib | sed "s/.lib[0-9]\+//g")
umi=$3

##Temp dir
TMP_DIR=/scratch/dmalload/${runname}
mkdir $TMP_DIR

dir=$(dirname $1)
cd $dir

#Makes unaligned bam with readgroup
if [[ ! -s ${runname}_markadapters.bam ]] && [[ ! -s ${runname}_unaligned.bam ]] && [[ ! -s ${runname}_unaligned_temp.bam ]] ##Working in the ifs
then
    rg=$(gunzip -c *_R1_*fastq.gz | head -1 | sed "s/^[^:]*:[^:]*:\([^:]*\):\([^:]*\).*/\1.\2/g") #Flowcell.lane
    pu="$rg.$samplelib" #Flowcell.lane.samplelib
    java -Xms512m -Xmx8G -jar $PICARDJAR FastqToSam F1=$R1 F2=$R2 O=${runname}_unaligned_temp.bam RG=$rg SM=$samplename PU=$pu LB=$samplelib PL=$pl CN=$cn
    if [[ -s ${runname}_unaligned_temp.bam ]]
    then
        rm -rf ${runname}*fastq.gz ##we have backups in the nas
    fi
fi

if [[ ! -s ${runname}_markadapters.bam ]] && [[ ! -s ${runname}_unaligned.bam ]] && [[ -s ${runname}_unaligned_temp.bam ]]
then
    if [[ $umi -eq 1 ]]
    then
        fgbio.sh --tmp-dir=$TMP_DIR --compression=0 ExtractUmisFromBam --input=${runname}_unaligned_temp.bam --output=/dev/stdout --read-structure=$read_structure --molecular-index-tags=ZA ZB --single-tag=RX | sambamba view -f bam --nthreads=$SLURM_CPUS_PER_TASK -o ${runname}_unaligned.bam /dev/stdin
        if [[ -s ${runname}_unaligned.bam ]]
        then
            rm -f ${runname}_unaligned_temp.bam
        fi
    else    
        mv ${runname}_unaligned_temp.bam ${runname}_unaligned.bam
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

rm -rf $TMP_DIR
