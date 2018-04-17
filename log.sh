#!/bin/bash
scriptdir=/home/dmalload/ppipilot/scripts
usage="$0 parentdir\n parentdir: directory with one directory per sample, with the two fastq files inside, with _R1_ and _R2_ in their name\n"

if [[ $# -ne 1 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
    echo -e $usage
    exit
fi

cd $1

for sdir in */
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
