#!/bin/bash
#SBATCH --mem 8000
#SBATCH -c 4

module load gatk/4.0.1.2

##Configuration variables

af="0.00000812" #1/123137 Exomes in genpop
genpop="/home/dmalload/ppipilot/gpop/af_gnomad2.0.2.vcf.gz"

usage="$0 tumor_bam output_dir"

if [[ $# -ne 2 ]] || [[ ! -d $2 ]] || [[ ! -f $1 ]]
then
    echo -e $usage
    exit
fi

cd $2
tumor=$1
output=$(basename $2)

tumorname=$(samtools view -H $tumor | sed -n '/@RG/p' | sed "s/.*SM:\([^ \t]*\).*/\1/g");
gatk --java-options "-Xmx8G" Mutect2 -R $HUMAN_GENOME -I $tumor --tumor-sample $tumorname \
--germline-resource $genpop --af-of-alleles-not-in-resource $af \
--output $output.vcf
