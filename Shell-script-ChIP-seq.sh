#!/bin/bash
# Assembly preparation
conda create -n ChIP \
-c bioconda -c conda-forge \
bowtie2 deeptools samtools macs2 trimmomatic=0.35 fastqc igvtools

# Data cleaning
trimmomatic-0.35.jar PE -threads 8 \
sample_1.clean.fq.gz sample_2.clean.fq.gz \
sample.R1.fq.gz sample_1.unmap.fq.gz sample.R2.fq.gz sample_2.unmap.fq.gz \
ILLUMINACLIP:adapters/TruSeq3-PE.fa:3:20:10:1:true \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:13 MINLEN:50 2> sample.PE.trim.log

trimmomatic-0.35.jar SE -threads 8 \
sample.fastq.gz \
sample.SE.fq.gz \
ILLUMINACLIP:adapters/TruSeq3-PE.fa:3:20:10:1:true \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:13 MINLEN:50 2> sample.SE.trim.log

# Quality cheak
fastqc -t 12 -o out_path sample_1.fq.gz sample_2.fq.gz
fastqc -t 12 -o out_path sample.fq.gz

# Mapping
bowtie2-build -f Ref/P_INF_CHROMOSOMES.fasta Ref/Pi
bowtie2-build -f Ref/Pm3010_flye.fasta Ref/Pm

bowtie2 -p 6 -3 5 --local -x Ref/Pi \
-U clean/Pisample-H3K27me3-IP.fq.gz \
-S results/Pisample-H3K27me3-IP.sam

bowtie2 -p 6 -3 5 --local -x Ref/Pi \
-U clean/Pisample-input.fq.gz \
-S results/Pisample-input.sam

bowtie2 -p 6 -3 5 --local -x Ref/Pm \
-U clean/Pmsample-H3K27me3-IP.fq.gz \
-S results/Pmsample-H3K27me3-IP.sam

bowtie2 -p 6 -3 5 --local -x Ref/Pm \
-U clean/Pmsample-input.fq.gz \
-S results/Pmsample-input.sam

# Converting to bam files and sorting
samtools sort results/name-H3K27me3-IP.sam \
-o results/name-H3K27me3-IP.sorted.bam

samtools index results/name-H3K27me3-IP.sorted.bam

samtools sort results/name-input.sam \
-o results/name-input.sorted.bam

samtools index results/name-input.sorted.bam

# Stat of peaks
macs2 callpeak -t results/Pisample-H3K27me3-IP.sorted.bam -c results/Pisample-input.sorted.bam -f BAM -g 227158000 --keep-dup all --outdir callpeak -name PiH3K27me3 --broad

macs2 callpeak -t results/Pmsample-H3K27me3-IP.sorted.bam -c results/Pmsample-input.sorted.bam -f BAM -g 156096338 --keep-dup all --outdir callpeak -name PmH3K27me3 --broad

intersectBed -a HSE1_gene.bed -b callpeak/name-H3K27me3.broadPeak -f 0.5 -r -wo

# Visualization of data
bamCompare -b1 results/name-H3K27me3-IP.sorted.bam -b2 results/name-input.sorted.bam -o log2ratio.bw
