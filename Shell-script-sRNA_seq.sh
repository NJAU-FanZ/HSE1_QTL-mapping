#!/bin/bash
# Assembly preparation
conda create -n sRNA \
-c bioconda -c conda-forge \
fastx_toolkit samtools bowtie deeptools

# Quality filtering
fastx_trimmer -f 6 -i input.fastq -o output.fastq

# Quality cheak
fastqc -t 12 -o out_path sample-sRNA.fq

# Mapping
bowtie2-build -f Ref/P_INF_CHROMOSOMES.fasta Ref/Pi
bowtie2-build -f Ref/Pm3010_flye.fasta Ref/Pm

bowtie2 -p 6 -3 5 --local -x Ref/Pi \
-U clean/Pisample-sRNA.fq.gz \
-S results/Pisample-sRNA-Pi.sam

bowtie2 -p 6 -3 5 --local -x Ref/Pm \
-U clean/Pmsample-sRNA.fq.gz \
-S results/Pmsample-sRNA-Pm.sam

# Converting to bam files and sorting
samtools view -bS sample_sRNA.sam > sample_sRNA.bam

samtools sort sample_sRNA.bam -o sample_sRNA.sorted.bam

samtools index sample_sRNA.sorted.bam

# Stat of peaks
ShortStack --bamfile sample_sRNA.sorted.bam \
--genomefile ref.fasta \
--locifile Interested_loci.txt \
--nohp --mincov 0.5rpmm

# Visualization of data
BamCoverage -p 16 -b sample_sRNA.sorted.bam \
--outFileName sample_coverage.sRNA.bw \
--normalizeUsing CPM
