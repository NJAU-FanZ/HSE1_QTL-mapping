# Assembly preparation
devtools::install_github("xiaolei-lab/rMVP")

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

## Visualization of mhtplot from NJAU-Cheng He
######################################################
###############The variables will be used##################
######################################################
chr.length <- c(22923084,21262383,16897382,16741615,14697857,14531184,13769070,13450935,
                13431482,13206068,12896474,12552869,11760846,10674240,10077107)
#use a 'for' cycle to produce several useful vectors.
for(i.chr in 1:15){
    if(i.chr==1){
        pos.elment1 <- sum(chr.length[1:i.chr])
        pos.elment2 <- pos.elment1/2
        pos <- pos.elment1
        pos.chr <- pos.elment2
    }else{
        pos.elment2 <- pos.elment1+chr.length[i.chr]/2
        pos.elment1 <- sum(chr.length[1:i.chr])
        pos <- c(pos,pos.elment1)
        pos.chr <- c(pos.chr,pos.elment2)}
}

######################################################
#############Using 'for' cycles to draw mahtplot#############
######################################################

   gwas.data <- read.table("FarmCPU.txt",sep="\t",header=TRUE)
   p <- gwas.data$p
   ps <- sort(p)
   order <- c(1:length(ps))/length(ps)*0.05
   bhcutoff <- ps[which((ps-order) > 0)[1]]
   mht.data.pro <- gwas.data[,c(1,2,3)]
   dimnames(mht.data.pro)[[2]] <- c("chr","pos","p")
   

   sig.data <- read.table("FarmCPU_signal.txt",sep="\t",header=TRUE)
   A <- sig.data$A
   As <- sort(A)
   order <- c(1:length(As))/length(As)*0.05
   bhcutoff1 <- As[which((As-order) > 0)[1]]
   signal.data.pro <- sig.data[,c(1,2,3)]
   dimnames(signal.data.pro)[[2]] <- c("chr","pos","p")
#use a second 'for' cycle to modify data structure.
   for(j in 2:15){
       if(j==2){
           chr.data <- mht.data.pro[mht.data.pro[,1]==j,]
           col2.data <- chr.data[,2]+pos[j-1]
           chr.data[,2] <-  col2.data
           mht.data <- rbind(mht.data.pro[mht.data.pro[,1]==1,],chr.data)
       }else{
           chr.data <- mht.data.pro[mht.data.pro[,1]==j,]
           col2.data <- chr.data[,2]+pos[j-1]
           chr.data[,2] <-  col2.data
           mht.data <- rbind(mht.data,chr.data)
            }
   }

   for(j in 2:15){
       if(j==2){
           chr.data <- signal.data.pro[signal.data.pro[,1]==j,]
           col2.data <- chr.data[,2]+pos[j-1]
           chr.data[,2] <-  col2.data
           signal.data <- rbind(signal.data.pro[signal.data.pro[,1]==1,],chr.data)
       }else{
           chr.data <- signal.data.pro[signal.data.pro[,1]==j,]
           col2.data <- chr.data[,2]+pos[j-1]
           chr.data[,2] <-  col2.data
           signal.data <- rbind(signal.data,chr.data)
            }
   }
#set some variables for mhtplot.
   pngfile.names <- paste("FarmCPU_213s_GWAS-t2.png",sep="")
   title.name <- paste("FarmCPU_213s_GWAS_mhtplot",sep="")
   threshold <- -log10(bhcutoff) ## FDR cutoff line
   threshold_l <- 6.25 ## cutoff line
   threshold_2 <- 6.50 ## Bonferroni矫正
   threshold_3 <- 6.32 ## Benjamini-Hochberg矫正
   label.name<-c(paste("chr.",1:15,sep=""))
   color.array<-rep(c("#BEBEBE","#BEBEBE"), 8)
  
#open a png file for plotting.
   png(pngfile.names,res=200,height=1500,width=3200)
   y <- -log10(mht.data[,3])
   x <- mht.data[,2]
   plot(x,y,type="p",cex=0,xlab="Chromosome",ylab="-log10(p-value)",xlim=c(0,pos[15]),
        ylim=c(0,10),xaxs="i",yaxs="i",xaxt="n",family="serif",font=2)       
#use a 'for' cycle to draw points with differnt colors.
   for(k in 1:15){
       x.zc <- mht.data[mht.data[,1]==k,2]
       y.zc <- -log10(mht.data[mht.data[,1]==k,3])
       points(x.zc,y.zc,col=color.array[k],pch=19,cex=1)
   }      
   for(m in 1:15){
       x.zc <- signal.data[signal.data[,1]==m,2]
       y.zc <- -log10(signal.data[signal.data[,1]==m,3])
       points(x.zc,y.zc,col="#8B2323",pch=19,cex=1.5)
   }      

# Modify and embellish the mhtplot.
   axis(1,at=pos.chr,labels=paste(1:15,sep=""),font=2)
   lines(x=c(0,pos[15]),y=c(threshold_3,threshold_3),lty=2,type="l",col="black",lwd=3)
   lines(x=c(0,pos[15]),y=c(0,0),type="l",col="black",lwd=1)
   title(title.name,cex.main=2)
   par(font.lab=2,font.axis=2)
   dev.off()

