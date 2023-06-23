#!/bin/bash

module load fastqc/0.11.3

cd $1

name=$(basename $1)

mkdir -p qc/fastqc
mkdir -p qc/bamqc
mkdir -p qc/metrics
mv *.metrics qc/metrics

fastqc ${name}_mdups.bam -o qc/fastqc
bamqc ${name}_mdups.bam -o qc/bamqc
samtools flagstat ${name}_mdups.bam > qc/samtools_flagstat.txt
