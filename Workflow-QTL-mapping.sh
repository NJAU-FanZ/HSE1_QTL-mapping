#!/bin/bash
# Assembly preparation  
conda create -n WGS \
-c bioconda -c conda-forge \
gatk4 bwa plink samtools=1.9 vcftools snpeff python=3.7 trimmomatic=0.35 fastqc

# Data cleaning 
trimmomatic-0.35.jar PE -threads 8 \
sample_1.clean.fq.gz sample_2.clean.fq.gz \
sample.R1.fq.gz sample_1.unmap.fq.gz sample.R2.fq.gz sample_2.unmap.fq.gz \
ILLUMINACLIP:adapters/TruSeq3-PE.fa:3:20:10:1:true \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:13 MINLEN:50 2> sample.PE.trim.log

# Quality cheak
fastqc -t 12 -o out_path sample_1.fq.gz sample_2.fq.gz

# Mapping  
bwa index Ref/P_INF_CHROMOSOMES.fasta

gatk CreateSequenceDictionary -R Ref/P_INF_CHROMOSOMES.fasta -O Ref/P_INF_CHROMOSOMES.dict  
samtools faidx Ref/P_INF_CHROMOSOMES.fasta

bwa mem -t 16 -M -B 3 \
-R '@RG\tID:name\tSM:name\tLB:name\tPL:ILLUMINA' Ref/P_INF_CHROMOSOMES.fasta \
CleanDATA/name_1.fq.gz CleanDATA/name_2.fq.gz \
-f bam/name.sam

samtools sort -o bam/name.bam bam/name.sam

gatk MarkDuplicates -I bam/name.bam \
-O bam/name.MarkDup.bam \
-M bam/name.markdup_metrics.txt

samtools view -h -q 1 -F 256 -F 2048 \
-b bam/name.MarkDup.bam > bam/name_uniq.bam

samtools index bam/name_uniq.bam

# Variant calling and filtering  
gatk HaplotypeCaller \
-R Ref/P_INF_CHROMOSOMES.fasta \
-I bam/name_uniq.bam \
-O bam/name.gvcf.gz --emit-ref-confidence GVCF -stand-call-conf 30 \
>bam/name.HaplotypeCaller.log  

gatk CombineGVCFs \
-R Ref/P_INF_CHROMOSOMES.fasta \
-V bam/name1.gvcf.gz \
...  \
-V bam/name115.gvcf.gz \
-O results/combine-samples.gvcf \
>results/combine-samples.gvcf.log  

gatk GenotypeGVCFs \
-R Ref/P_INF_CHROMOSOMES.fasta \
-V results/combine-samples.gvcf \
-O results/combine-samples.vcf \
>results/combine-samples.GenotypeGvcf.log

# SNP filtering and GWAS  
gatk CombineGVCFs \
  -R Ref/P_INF_CHROMOSOMES.fasta \
  -V bam/name1.gcvf \
  -V bam/name2.gcvf \
...
  -V bam/name115.gcvf \
  -O results/combine-samples.gvcf \
  >results/combine-samples.gvcf.log
gatk GenotypeGVCFs \
  -R Ref/P_INF_CHROMOSOMES.fasta  \
  -V results/combine-samples.gvcf \
  -O results/combine-samples.vcf \
  >results/combine-samples.GenotypeGvcf.log
  
#SNP filting
gatk SelectVariants -select-type SNP \
  -V results/combine-samples.vcf \
  -O results/combine-samples.snp.vcf \
  >results/combine-samples.snp.log
  
gatk SelectVariants -select-type INDEL \
  -V results/combine-samples.vcf \
  -O results/combine-samples.indel.vcf \
  >results/combine-samples.indel.log
  
gatk VariantFiltration -V results/combine-samples.snp.vcf \
  -O results/combine-samples-filt.snp.vcf \
  -filter "QD < 2.0" --filter-name "QD2"  \
  -filter "QUAL < 30.0" --filter-name "QUAL30"  \
  -filter "SOR > 3.0" --filter-name "SOR3"  \
  -filter "FS > 60.0" --filter-name "FS60"  \
  -filter "MQ < 40.0" --filter-name "MQ40"  \
  >results/combine-snp-samples.filter.log
  
gatk VariantFiltration -V results/combine-samples.indel.vcf \
  -O results/combine-samples-filt.indel.vcf \
  -filter "QD < 2.0" --filter-name "QD2"  \
  -filter "QUAL < 30.0" --filter-name "QUAL30"  \
  -filter "SOR > 3.0" --filter-name "SOR3"  \
  -filter "FS > 60.0" --filter-name "FS60"  \
  -filter "MQ < 40.0" --filter-name "MQ40"  \
  >results/combine-indel-samples.filter.log
  
gatk MergeVcfs -I results/combine-samples-filt.indel.vcf \
  -I results/combine-samples-filt.snp.vcf \
  -O results/Filter-combine-samples.vcf \
  >results/Filter-combine-samples.log
  
vcftools --vcf results/Filter-combine-samples.vcf \
  --max-missing 0.95 \
  --maf 0.05 \
  --min-meanDP 5 \
  --recode \
  --recode-INFO-all \
  --remove-filtered-all \
  --out results/Filter-samples
  
vcftools --vcf results/Filter-samples.recode.vcf \
  --plink --out results/Filter-samples
  
plink --vcf results/Filter-samples.recode.vcf \
--maf 0.05 --geno 0.1 --recode vcf-iid \
--out results/Filter-samples-1 --allow-extra-chr

plink --vcf results/Filter-samples-1.vcf \
--recode 12 transpose --output-missing-genotype 0 \
--out results/Filter-samples-2 --autosome-num 90 --allow-extra-chr

emmax-kin-intel64 results/Filter-samples-2 \
-v -d 10 -o results/Filter-samples.BN.kinf

# Parameters for linkage disequilibrium analysis
LDBlockShow -InVCF Filter-samples.vcf \
-OutPut out-20k \
-Region chr:loci-10k:loci+10k -OutPng -SeleVar 2

ShowLDSVG -InPreFix out-20k \
-OutPut chrN-20k \
-InGWAS gwaschrN.pvalue -Cutline 7 -ShowNum 100 -PointSize 30 -OutPng
