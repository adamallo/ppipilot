#!/bin/bash
#SBATCH --mem 8000
#SBATCH -c 4

module load gatk/4.0.1.2

##Configuration variables

af="0.00000812" #1/123137 Exomes in genpop
genpop="/home/dmalload/ppipilot/gpop/af_gnomad2.0.2.vcf.gz"

usage="$0 tumor_bam normal_bam pon output_dir"

if [[ $# -ne 4 ]] || [[ ! -d $4 ]] || [[ ! -f $3 ]] || [[ ! -f $2 ]] || [[ ! -f $1 ]]
then
    echo -e $usage
    exit
fi

cd $3

normal=$2
tumor=$1

output=$(basename $4)

tumorname=$(samtools view -H $tumor | sed -n '/@RG/p' | sed "s/.*SM:\([^ \t]*\).*/\1/g");
gatk --java-options "-Xmx8G" Mutect2 -R $HUMAN_GENOME -I $normal -I $tumor --tumor-sample $tumorname \
-pon $3 --germline-resource $genpop --af-of-alleles-not-in-resource $af \
--output $outputdir/$output.vcf.gz -bamout $outputdir/$output.bam
