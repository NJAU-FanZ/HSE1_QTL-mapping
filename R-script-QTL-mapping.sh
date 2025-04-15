
# R package rMVP analysis
library(rMVP)
MVP.Data(filePhe="pheno.txt",fileVCF="Filter-samples.vcf",out="pop")
MVP.Data.Kin("name.BN.kinf", out="pop", maxLine=1e4, sep='\t')
genotype <- attach.big.matrix("pop.geno.desc")
phenotype <- read.table("pop.phe",head=TRUE)
map <- read.table("pop.geno.map" , head = TRUE)
Kinship <- attach.big.matrix("pop.kin.desc")
imMVP <- MVP(
    phe=phenotype,
    geno=genotype,
    map=map,
    K=Kinship,
    #CV.GLM=Covariates,  ##if you have additional covariates, please keep there open.
    #CV.MLM=Covariates,
    #CV.FarmCPU=Covariates,
    nPC.GLM=5,   ##if you have added PCs into covariates, please keep there closed.
    nPC.MLM=3,  ##if you don't want to add PCs as covariates, please comment out the parameters instead of setting the nPC to 0.
    nPC.FarmCPU=3,
    maxLine=10000,   #smaller value would reduce the memory cost
    #ncpus=10,
    vc.method="BRENT",  ##only works for MLM
    maxLoop=10,
    method.bin="static",   ## "FaST-LMM", "static" (#only works for FarmCPU)
    #permutation.threshold=TRUE,
    #permutation.rep=100,
    threshold=0.05,
    method=c("GLM", "MLM", "FarmCPU")
)
str(imMVP)
