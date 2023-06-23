#!/bin/bash
#SBATCH --mem-per-cpu=8192

module load gatk/4.0.1.2

##Configuration variables
minimum_allele_frequency="0.05"

#Adapted from https://github.com/broadinstitute/gatk/blob/master/scripts/mutect2_wdl/mutect_resources.wdl

cd $1
name=$2

if [[ $(ls *.vcf* | wc -l) -ne 1 ]]
then
    wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz
fi

input_vcf=$(ls *.vcf.bgz)

if [[ ! -f header ]]
then
    gunzip -c $input_vcf | grep '#' > header
fi

if [[ ! -f body ]]
then 
    gunzip -c $input_vcf | grep -v '#' | grep -P "\tPASS\t" > body
fi

# delete any INFO fields before AF (PASS<TAB><other info fields>;AF=___ --> PASS<TAB>AF=____)
# delete any INFO fields after AF (;<other info fields><end of line> --> nothing)

if [[ ! -f simplified_info ]]
then
    sed -e 's/PASS\t.*AF=/PASS\tAF=/g' -e 's/[;].*$//g' body > simplified_info
    #rm -f body
fi

# replace ID (3rd) and QUAL (6th) columns with '.' (empty)

if [[ ! -f simplified_body ]]
then
    while read contig pos id ref alt qual filter info; do
        printf "$contig\t$pos\t.\t$ref\t$alt\t.\t$filter\t$info\n"
    done < simplified_info > simplified_body
    #rm -f simplified_info
fi

if [[ ! -f simplified.vcf ]]
then
    cat header simplified_body > simplified.vcf
    #rm -f simplified_body
    #rm -f header
fi


gatk --java-options -Xmx8G IndexFeatureFile -F simplified.vcf
gatk --java-options -Xmx8G SelectVariants -V simplified.vcf -O $name.vcf.gz

##variants_for_contamination

gatk --java-options -Xmx8G SelectVariants -V $name.vcf.gz \
-select-type SNP -restrict-alleles-to BIALLELIC \
-select "AF > ${minimum_allele_frequency}" \
-O ${name}_vcontaminant_minAF${minimum_allele_frequency}.vcf.gz \
--lenient

rm -f simplified.vcf simplified.vcf.idx
