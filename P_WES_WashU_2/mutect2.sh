#!/bin/bash
#SBATCH --mem 8000
#SBATCH -c 4

module load gatk/4.0.1.2

##Configuration variables

af="0.00000812" #1/123137 Exomes in genpop
genpop="/home/dmalload/ppipilot/gpop/af_gnomad2.0.2.vcf.gz"

usage="$0 tumor_bam normal_bam pon output_dir [bed]"

if [[ $# -lt 4 ]] || [[ ! -d $4 ]] || [[ ! -f $3 ]] || [[ ! -f $2 ]] || [[ ! -f $1 ]]
then
    echo -e $usage
    exit
fi

limits=""
if [[ $# -gt 4 ]] && [[ $5 != "" ]]
then
    limits="-L "$(readlink -e $5)
fi

pon=$(readlink -e $3)
normal=$(readlink -e $2)
tumor=$(readlink -e $1)
output=$(basename $4)
outputdir=$(readlink -e $4)
cd $4

tumorname=$(samtools view -H $tumor | sed -n '/@RG/p' | sed "s/.*SM:\([^ \t]*\).*/\1/g");
gatk --java-options "-Xmx8G" Mutect2 -R $HUMAN_GENOME -I $normal -I $tumor --tumor-sample $tumorname \
-pon $pon --germline-resource $genpop --af-of-alleles-not-in-resource $af \
--output $outputdir/$output.vcf.gz -bamout $outputdir/$output.bam ${limits}
