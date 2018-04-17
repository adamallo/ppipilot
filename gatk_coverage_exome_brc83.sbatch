#!/bin/bash
#SBATCH --mem 10000

module load gatk/3.5.0
name=$(dirname $1)/coverage_exome_brc83

java -Xms512m -Xmx8G -jar /packages/6x/gatk/3.5.0/GenomeAnalysisTK.jar -T DepthOfCoverage -R $HUMAN_GENOME -L /home/dmalload/ppipilot/beds/xgen-exome-research-panel-targets.bed -L /home/dmalload/ppipilot/beds/BRC83_probes_3769b74109d6438f9140353bec5ac2f7.bed -I $1 -o $name --summaryCoverageThreshold 50 --summaryCoverageThreshold 40 --summaryCoverageThreshold 30 --summaryCoverageThreshold 20 --summaryCoverageThreshold 10
