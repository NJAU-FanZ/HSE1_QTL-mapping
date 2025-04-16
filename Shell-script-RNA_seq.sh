#!/bin/bash
# Assembly preparation
conda create -n RNA \
-c bioconda -c conda-forge \
hisat2 samtools htseq python=3.7 trimmomatic=0.35 fastqc igvtools

# Data cleaning
trimmomatic-0.35.jar PE -threads 8 \
sample_1.clean.fq.gz sample_2.clean.fq.gz \
sample.R1.fq.gz sample_1.unmap.fq.gz sample.R2.fq.gz sample_2.unmap.fq.gz \
ILLUMINACLIP:adapters/TruSeq3-PE.fa:3:20:10:1:true \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:13 MINLEN:50 2> sample.PE.trim.log

# Quality cheak
fastqc -t 12 -o out_path sample_1.fq.gz sample_2.fq.gz

# Mapping
extract_exons.py RefSeq.gtf > genome.exon
extract_splice_sites.py RefSeq.gtf > genome.ss

hisat2-build Ref.fasta --ss genome.ss --exon genome.exon Ref_tran

hisat2 -p 10 -x Ref_tran -1 sample.R1.fq.gz -2 sample.R2.fq.gz -S sample.hisat.sam 2>sample.hisat.align.log

# Converting to bam files and sorting
samtools sort sample.hisat.sam \
-o sample.sorted.bam

samtools index sample.sorted.bam

# Calculate the reads of each gene
htseq-count -f bam -r name -s no -t exon \
-i transcript_id -m intersection-nonempty \
sample.sorted.bam RefSeq.gtf \
>sample.counts_out.txt
