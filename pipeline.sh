#!/bin/bash
scriptdir=/home/dmalload/ppipilot/scripts
usage="$0 fastqdir parentdir umi [bed]\n fastqdir: directory with fastqfiles\nparentdir: directory in which new directories for each case will be generated\n umi: binary flag indicating whether reads have 8-bp UMIs or not\n"

if [[ $# -lt 3 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ ! -d $1 ]]
then
    echo -e $usage
    exit
fi

mkdir -p $2

basedir=$2
fastqdir=$1
umi=$3
depsep=":"
dependency="--dependency=afterok:"

bed=""

if [[ $# -gt 3 ]] && [[ -s $4 ]]
then
    bed=$(readlink -e $4)
fi

#Example of fastq file naming format
#H_ZT-591-3-lib1_S14_L002_I1_001.fastq.gz
#H_([^_]*)_([^_]*)_([^_]*)_([^_]*)_([^_]*).fastq.gz \1: sample-libraryname \2: sample? \3: lane \4: IDreads [I1,I2,R1,R2] \5: SomeKindOfId

##File reorganization

casedirs=($(ls $fastqdir/* | xargs -L 1 basename | sed "s/\([^_]*\)_\([^_]*\)_\([^_]*\)_\([^_]*\)_\([^_]*\)_\([^_]*\).fastq.gz/\2/g" | sed "s/.lib[0-9]\+//g" | sort | uniq))

for casedir in "${casedirs[@]}"
do
    mkdir -p $basedir/$casedir
    cp -n $fastqdir/*$casedir* $basedir/$casedir
done

declare -A jobs

for singlerun in $basedir/*/*_R1_*
do
    mate=$(echo $singlerun | sed "s/_R1_/_R2_/")
    if [[ -f $singlerun ]] && [[ -f $mate ]]
    then
        jobid=$(submit $scriptdir/makeBamPerRun.sh $singlerun $mate $umi | sed "s/Queued job \([0-9]\+\)/\1/")
        this_casedir=$(dirname $singlerun | xargs basename)
        if [[ ${jobs["$this_casedir"]} == "" ]]
        then
            jobs["$this_casedir"]=$dependency$jobid
        else
            jobs["$this_casedir"]=${jobs["$this_casedir"]}$depsep$jobid
        fi
        echo "Submitting makeBamPerRun job for $singlerun and $mate for case $this_casedir with job_id $jobid\n"
    else
        echo "The mate $mate for $singlerun cannot be found!"
        exit 1
    fi
done

for casename in "${casedirs[@]}"
do
    casedir=$basedir/$casename
    echo "Submitting combineBamPerSample job for sample $casename with dependency" ${jobs["$casename"]}
    job=$(submit -c 8 ${jobs["$casename"]} $scriptdir/combineBamPerSample.sh $casedir $umi $bed | sed "s/Queued job \([0-9]\+\)/\1/") 
    
    #launchpileupsummaries depends on job
    #need to store this jobs by pair

done

#PON depends on all jobs of casenames that are normal

##SNV
#foreach pair
#makecontamination depends on launchpileupsummaries of pair
#makemutect2 depends on PON and launchcontamination
#filtercalls depends on makemutect2


##CNV
