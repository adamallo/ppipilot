#!/bin/bash
scriptdir=/home/dmalload/ppipilot/scripts
usage="$0 parentdir umi\n parentdir: directory with fastq files, two per sample with _R1_ and _R2_ in their name\n umi: binary flag indicating whether reads have 8-bp UMIs or not\n"

if [[ $# -ne 2 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ ! -d $1 ]]
then
    echo -e $usage
    exit
fi

dir=$1
umi=$2
dirbase=$(dirname $dir)
depsep=":"
dependency="--dependency=afterok:"

#Example of fastq file naming format
#H_ZT-591-3-lib1_S14_L002_I1_001.fastq.gz
#H_([^_]*)_([^_]*)_([^_]*)_([^_]*)_([^_]*).fastq.gz \1: sample-libraryname \2: sample? \3: lane \4: IDreads [I1,I2,R1,R2] \5: SomeKindOfId

##File reorganization

casedirs=($(ls $dir/* | xargs -L 1 basename | sed "s/\([^_]*\)_\([^_]*\)_\([^_]*\)_\([^_]*\)_\([^_]*\)_\([^_]*\).fastq.gz/\2/g" | sed "s/.lib[0-9]\+//g" | sort | uniq))

for casedir in "${casedirs[@]}"
do
    mkdir $casedir
    cp $dir/*$casedir* $casedir
done

declare -A jobs

for singlerun in $dir/*/*_R1_*
do
    mate=$(echo $singlerun | sed "s/_R1_/_R2_/")
    if [[ -f $singlerun ]] && [[ -f $mate ]]
    then
        jobid=$(submit $scriptdir/makeBamPerRun.sh $singlerun $mate $umi | sed "s/Queued job \([0-9]\+\)/\1/")
        jobs+=(["$singlerun"]="$jobid") ##This has to be for each casedir, not singlerun
        echo "Submitting makeBamPerRun job for $singlerun and $mate with job_id $jobid\n"
    else
        echo "The mate $mate for $singlerun cannot be found!"
        exit 1
    fi
done


for casedir in "${casedirs[@]}"
do
    ##Merge and either markduplicates or umi stuff

    if [[ $(ls $casedir/*.fastq.gz | wc -l) -eq 2 ]]
    then
       
        echo "Submitting make_bam_gatk job for sample $sample"
        #sbatch -p private $scriptdir/make_bam_gatk.sh $sample
        submit $scriptdir/make_bam_gatk.sh $(basename $casedir)
    else
        echo "Skipping bam preparation for $casedir"
    fi
done


for sdir in $dir/
do
    sample=$(basename $sdir)
    if [[ $(ls $sdir*.fastq.gz | wc -l) -eq 2 ]]
    then
       
        echo "Submitting make_bam_gatk job for sample $sample"
        #sbatch -p private $scriptdir/make_bam_gatk.sh $sample
        submit $scriptdir/make_bam_gatk.sh $sample
    else
        echo "Skipping bam preparation for $sdir "
    fi

    ##I need to add dependencies here!
 
    if [[ -f $sdir/${sample}_mdups.bam ]] && [[ ! -d $sdir/qc ]]
    then 
        echo "Submitting qc job for sample $sample"
        #sbatch -p private $scriptdir/qc.sh $sample
        submit $scriptdir/qc.sh $sample
    else
        echo "Skipping qc for $sdir "
    fi
done
