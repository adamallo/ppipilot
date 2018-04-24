#!/bin/bash
#SBATCH -c 8
scriptdir=/home/dmalload/ppipilot/scripts

module load picard/2.9.2
set -o pipefail

submit="submit"

usage="$0 sampledir umi\n umi: binary flag, 1 includes a 8bp UMI, 0 no UMI "
umi=1


##Settings with a highest probability of being changed
UMIminMAPQ=30 ###Filtering for UMI calling (In the IDT manual they use 10, the default is 30). Here I am more stringent since this info cannot be recovered from the consensus and may influence strongly their estimation
UMIminPHRED=10 ###Filtering for UMI calling (In the IDT manual they use 30, the default is 10). This info cannot be recovered from consensus reads. However, it is less important than UMIminMAPQ since even if we include bases with low quality they should be "consensed out" if we have enough number of reads comming from the same molecule
UMIminreads=1 ###I can filter them a posteriori, doing it here just improves computation time, but I do not think it is worth losing data considering low coverage in some samples

if [[ $# -ne 2 ]] || [[ ! -d $1 ]]
then
    echo -e $usage
    exit 1
fi

sampledir=$1
sample=$(basename $sampledir)
umi=$2

##Temp dir
TMP_DIR=/scratch/dmalload/$(basename $sampledir)
mkdir $TMP_DIR

cd $sampledir

##Only for umis
if [[ $umi -eq 1 ]]
then

    ###Merge all bams (from different runs/lanes and/or libraries)
    inputlist="" 
    for i in *.bam
    do
        inputlist="${inputlist}I=$i "
    done
    
    if [[ $inputlist == "" ]]
    then
        echo "No bam files detected in $sampledir"
        exit 1
    fi
    
    ##Reads pairs are sorted by 5' position of the two reads, library, UMI and read name.
    ##Umis can accumulate mutations that are then accumulated following a directed graph (pcr cycles). This is taken into consideration in the paired strategy. It is based on making clusters of reads with a given edit distance (1 seems optimal) and directed connections that explain groups with that mutation and others in less number of reads. If this wasn't done estimates would be biased since there would be a number of consensus reads from the same molecule
    java -Xmx8G -jar $PICARDJAR MergeSamFiles ${inputlist}O=${sample}_merged.bam USE_THREADING=True SORT_ORDER=coordinate

    ##For some reason I cannot pipe the output of MergeSamFiles to fgbio

    fgbio.sh --compression=0 --tmp-dir=$TMP_DIR GroupReadsByUmi --input=${sample}_merged.bam --output=/dev/stdout --strategy=paired --raw-tag=RX --assign-tag=MI --min-map-q=$UMIminMAPQ --edits=1 | sambamba view -f bam --nthreads=$SLURM_CPUS_PER_TASK -o ${sample}_grouped.bam /dev/stdin
 
    #I could pipe these two, but I need the output for CollectDuplexSeqMetrics since it does not work well if mates are in different chromosomes (I guess from secondary alignments
    
    fgbio.sh --compression=0 --tmp-dir=$TMP_DIR CallDuplexConsensusReads --input=grouped.bam --output=/dev/stdout --min-input-base-quality=$UMIminPHRED --min-reads=$UMIminreads --sort-order=queryname | sambamba view -f bam --nthreads=$SLURM_CPUS_PER_TASK -o ${sample}_unalignedConsensus.bam /dev/stdin

    #Submit aditional job for 
    
    $submit $scriptdir/duplexmetrics.sh ${sample}_grouped.bam duplex_metrics/${sample} 1
 
    ##Sorting order must be queryname here, otherwise SamToFastq fails due to lack of ram (with any amount of RAM really)

    ##No clipping this time since it was done previously (right after marking adapters)    
    java -Xmx4G -jar $PICARDJAR SamToFastq \
    I=${sample}_unalignedConsensus.bam \
    FASTQ=/dev/stdout \
    INTERLEAVE=true NON_PF=true \
    TMP_DIR=$TMP_DIR VALIDATION_STRINGENCY=SILENT | \
    bwa mem -M -t $SLURM_CPUS_PER_TASK -p $HUMAN_GENOME /dev/stdin | \
    java -Xmx8G -jar $PICARDJAR MergeBamAlignment \
    ALIGNED_BAM=/dev/stdin \
    UNMAPPED_BAM=${sample}_unalignedConsensus.bam \
    OUTPUT=${sample}.bam \
    R=$HUMAN_GENOME CREATE_INDEX=true ADD_MATE_CIGAR=true \
    CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
    INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
    ATTRIBUTES_TO_RETAIN=RG \
    ATTRIBUTES_TO_RETAIN=RX \
    ATTRIBUTES_TO_RETAIN=MI \
    ATTRIBUTES_TO_RETAIN=aD \
    ATTRIBUTES_TO_RETAIN=bD \
    ATTRIBUTES_TO_RETAIN=cD \
    ATTRIBUTES_TO_RETAIN=aM \
    ATTRIBUTES_TO_RETAIN=bM \
    ATTRIBUTES_TO_RETAIN=cM \
    ATTRIBUTES_TO_RETAIN=aE \
    ATTRIBUTES_TO_RETAIN=bE \
    ATTRIBUTES_TO_RETAIN=cE \
    ATTRIBUTES_TO_REVERSE=ad \
    ATTRIBUTES_TO_REVERSE=bd \
    ATTRIBUTES_TO_REVERSE=ae \
    ATTRIBUTES_TO_REVERSE=be \
    ATTRIBUTES_TO_REVERSE=ac \
    ATTRIBUTES_TO_REVERSE=bc \
    ATTRIBUTES_TO_REVERSE=aq \
    ATTRIBUTES_TO_REVERSE=bq \
    TMP_DIR=$TMP_DIR EXPECTED_ORIENTATIONS=FR SO=coordinate \
    VALIDATION_STRINGENCY=SILENT
    ##I am clipping overlapping reads here. The IDT tutorial recomends to hard-clip them using fgbio's ClipBam. I do not see the point on that and actually they are doing it twice (since the default here is true). I may want to come back to this in the future and see if it changes something. I suspect they are just hard-clipping them for VarDict
    ##I am including secondary alignments like GATK/Mutect recommend. Not in IDT tutorial
    ##I am retaining all statistics from consensus calling + MI (to get UMI statistics from this bam file). Reversing (if needed) all per-base attributes
    
    if [[ -s $sample.bam ]]
    do
        rm -f $(echo $inputlist | sed "s/I=//g")
        rm -f ${sample}_unalignedConsensus.bam
        ##I only need the merged bam (to use without UMI consensuation) and the $sample.bam
    done
else 
    #Markduplicates ##Not necessary for umis since they are being consensuated 
    if [[ ! -f ${sample}_mdups.bam ]]
    then

        ###Merge all bams (from different runs/lanes and/or libraries)
        inputlist=""
        
        for i in *.bam
        do
            inputlist="${inputlist}I=$i "
        done
        
        if [[ $inputlist == "" ]]
        then
            echo "No bam files detected in $sampledir"
            exit 1
        fi
#NOT tested
        ##Merging bams with sorting queryname for markduplicates, markduplicates and sorting back to coordinate for Mutect
        java -Xmx8G -jar $PICARDJAR MergeSamFiles ${inputlist}O=/dev/stdout USE_THREADING=True SORT_ORDER=queryname | java -Xmx8G -jar $PICARDJAR MarkDuplicates I=/dev/stdin O=/dev/stdout ASSUME_SORTED=true M=${sample}_markduplicates.metrics CREATE_INDEX=false TMP_DIR=$TMP_DIR | sambamba sort -t $SLURM_CPUS_PER_TASK -m ${SLURM_MEM_PER_NODE}M -o ${sample}.bam /dev/stdin
        
        if [[ -s $sample.bam ]]
        then
            rm -f $(echo $inputlist | sed "s/I=//g")
        fi 
fi

rm -rf $TMP_DIR
