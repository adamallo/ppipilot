#!/bin/bash
#SBATCH -c 8
scriptdir=/home/dmalload/ppistudy/scripts

module load picard/2.9.2
set -o pipefail

submit="submit"

usage="$0 sampledir umi [bed]\n umi: binary flag, 1 includes a 8bp UMI, 0 no UMI "
umi=1

#############################
##TODO:TMP_DIR makes markduplicates extremely slow. I have to doublecheck if this a common problem for other programs used in this pipeline
#############################

##Settings with a highest probability of being changed
UMIminMAPQ=50 ###Filtering for UMI calling (In the IDT manual they use 10, the default is 30, looking at the empirical distribution of our data basically all of them are at 60). Here I am more stringent since this info cannot be recovered from the consensus and may influence strongly their estimation
UMIminPHRED=20 ###Filtering for UMI calling (In the IDT manual they use 30, the default is 10). This info cannot be recovered from consensus reads. However, it is less important than UMIminMAPQ since even if we include bases with low quality they should be "consensed out" if we have enough number of reads comming from the same molecule and these phreds are taking into consideration to calculate the phred score of the consensus position.
UMIminreads=1 ###I can filter them a posteriori, doing it here just improves computation time, but I do not think it is worth losing data considering low coverage in some samples

if [[ $# -lt 2 ]] || [[ ! -d $1 ]]
then
    echo -e $usage
    exit 1
fi

sampledir=$1
sample=$(basename $sampledir)
umi=$2

bed=""

if [[ $# -gt 2 ]] && [[ -s $3 ]]
then
    bed=$(readlink -e $3)
fi

if [[ "$SLURM_CPUS_PER_TASK" == "" ]]
then
    SLURM_CPUS_PER_TASK=1
fi

##Temp dir
TMP_DIR=/scratch/dmalload/$(basename $sampledir)
mkdir -p $TMP_DIR

cd $sampledir

##Only for umis
if [[ $umi -eq 1 ]]
then
    inputlist="" 
    
    ##Reads pairs are sorted by 5' position of the two reads, library, UMI and read name.
    ##Umis can accumulate mutations that are then accumulated following a directed graph (pcr cycles). This is taken into consideration in the paired strategy. It is based on making clusters of reads with a given edit distance (1 seems optimal) and directed connections that explain groups with that mutation and others in less number of reads. If this wasn't done estimates would be biased since there would be a number of consensus reads from the same molecule

    if [[ ! -s ${sample}_merged.bam ]]
    then

        ###Merge all bams (from different runs/lanes and/or libraries)
        if [[ $(ls *.bam | wc -l) -eq 1 ]]
        then
            filename=$(ls *.bam)
            mv $filename ${sample}_merged.bam
        else
            for i in *.bam
            do
                inputlist="${inputlist}I=$i "
            done
            
            if [[ $inputlist == "" ]]
            then
                echo "No bam files detected in $sampledir"
                exit 1
            fi
            
            java -Xmx8G -jar $PICARDJAR MergeSamFiles ${inputlist}O=${sample}_merged.bam USE_THREADING=True SORT_ORDER=coordinate
        fi

        #QC
        $submit $scriptdir/bamqc.sh ${sample}_merged.bam 
        $submit $scriptdir/fastqc.sh ${sample}_merged.bam
        $submit $scriptdir/flagstat.sh ${sample}_merged.bam
    fi

    if [[ ! -s ${sample}_grouped.bam ]]
    then
        ##For some reason I cannot pipe the output of MergeSamFiles to fgbio and I cannot use sambamba to make the bam or otherwise I need to re-sort the output?
        fgbio.sh --tmp-dir=$TMP_DIR GroupReadsByUmi --input=${sample}_merged.bam --output=${sample}_grouped.bam --strategy=paired --raw-tag=RX --assign-tag=MI --min-map-q=$UMIminMAPQ --edits=1
        #fgbio.sh --compression=0 --tmp-dir=$TMP_DIR GroupReadsByUmi --input=${sample}_merged.bam --output=/dev/stdout --strategy=paired --raw-tag=RX --assign-tag=MI --min-map-q=$UMIminMAPQ --edits=1 | sambamba view -f bam --nthreads=$SLURM_CPUS_PER_TASK -o ${sample}_grouped.bam /dev/stdin
    fi

    if [[ ! -s ${sample}.bam ]]
    then 
        #I could pipe with GroupReadsByUmi, but I need the output for CollectDuplexSeqMetrics since it does not work well if mates are in different chromosomes (I guess from secondary alignments
    
        fgbio.sh --compression=0 --tmp-dir=$TMP_DIR CallDuplexConsensusReads --input=${sample}_grouped.bam --output=/dev/stdout --min-input-base-quality=$UMIminPHRED --min-reads=$UMIminreads --sort-order=queryname | sambamba view -f bam --nthreads=$SLURM_CPUS_PER_TASK -o ${sample}_unalignedConsensus.bam /dev/stdin     
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
        #QC
        $submit $scriptdir/bamqc.sh ${sample}.bam 
        $submit $scriptdir/fastqc.sh ${sample}.bam
        $submit $scriptdir/flagstat.sh ${sample}.bam
        $submit $scriptdir/depth_hist.sh ${sample}.bam $bed 0
        $submit $scriptdir/duplexmetrics.sh ${sample}_grouped.bam qc_$sample/duplex_metric $sample 0 ##I may want to switch this to 1, but I may also want the input file?
        
        if [[ -s "$sample.bam" ]]
        then
            rm -f $(echo $inputlist | sed "s/I=//g")
            rm -f ${sample}_unalignedConsensus.bam
            ##I only need the merged bam (to use without UMI consensuation) and the $sample.bam
        fi
    fi
else 
    #Markduplicates ##Not necessary for umis since they are being consensuated 
    if [[ ! -f ${sample}.bam ]]
    then
        if [[ ! -s ${sample}_merged.bam ]]
        then
            inputlist=""
            
            ###Merge all bams (from different runs/lanes and/or libraries)
            if [[ $(ls *.bam | wc -l) -eq 1 ]]
            then
                filename=$(ls *.bam)
                mv $filename ${sample}_merged.bam
            else
                for i in *.bam
                do
                    inputlist="${inputlist}I=$i "
                done
                
                if [[ $inputlist == "" ]]
                then
                    echo "No bam files detected in $sampledir"
                    exit 1
                fi
                ##Merging bams with sorting queryname for markduplicates, markduplicates and sorting back to coordinate for Mutect
                java -Xmx8G -jar $PICARDJAR MergeSamFiles ${inputlist}O=${sample}_merged.bam USE_THREADING=True SORT_ORDER=queryname
            fi
        fi

        if [[ ! -s ${sample}.bam ]]
        then
            #java -Xmx8G -jar $PICARDJAR MarkDuplicates I=${sample}_merged.bam O=/dev/stdout M=${sample}_markduplicates.metrics CREATE_INDEX=false TMP_DIR=$TMP_DIR VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=0 | sambamba sort -t $SLURM_CPUS_PER_TASK -m $(($SLURM_MEM_PER_CPU*$SLURM_CPUS_PER_TASK-8192))M -o ${sample}.bam /dev/stdin ##Adding the TMP_DIR makes markduplicates extremely slow for no apparent reason

            java -Xmx8G -jar $PICARDJAR MarkDuplicates I=${sample}_merged.bam O=/dev/stdout M=${sample}_markduplicates.metrics CREATE_INDEX=false VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=0 | sambamba sort -t $SLURM_CPUS_PER_TASK -m $(($SLURM_MEM_PER_CPU*$SLURM_CPUS_PER_TASK-8192))M -o ${sample}.bam /dev/stdin

            $submit $scriptdir/bamqc.sh ${sample}.bam 
            $submit $scriptdir/fastqc.sh ${sample}.bam
            $submit $scriptdir/flagstat.sh ${sample}.bam
            $submit $scriptdir/depth_hist.sh ${sample}.bam $bed 0
            $submit $scriptdir/depth_hist.sh ${sample}.bam $bed 1
        fi
        
        if [[ -s $sample.bam ]]
        then
            rm -f $(echo $inputlist | sed "s/I=//g")
        fi
    fi 
fi

rm -rf $TMP_DIR

