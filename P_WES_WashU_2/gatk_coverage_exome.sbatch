#!/bin/bash
#SBATCH --mem 10000

module load gatk/3.5.0
name=$(dirname $1)/coverage_exome

java -Xms512m -Xmx8G -jar /packages/6x/gatk/3.5.0/GenomeAnalysisTK.jar -T DepthOfCoverage -R $HUMAN_GENOME -L /home/dmalload/ppipilot/beds/xgen-exome-research-panel-targets.bed -I $1 -o $name --summaryCoverageThreshold 50 --summaryCoverageThreshold 40 --summaryCoverageThreshold 30 --summaryCoverageThreshold 20 --summaryCoverageThreshold 10
