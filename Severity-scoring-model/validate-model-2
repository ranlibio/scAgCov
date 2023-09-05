#!/usr/bin/Rscript
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(ggsci)
library(ggpubr)
library(ggrepel)
library(forcats)

##combine count
files <- list.files("./")

geneInfor <- read.table(files[1],row.names = 1,header = TRUE,check.names = FALSE)
geneInfor <- geneInfor[,-c(1,2,3,4,7)]

countMatrix <- matrix(NA,nrow = nrow(geneInfor),ncol = length(files))
rownames(countMatrix) <- rownames(geneInfor)
colnames(countMatrix) <- sapply(files,function(x) unlist(strsplit(x,"_"))[1])

i <- 1
for (i in 1:length(files)) {
  tmp <- read.table(files[i],row.names = 1,header = TRUE,check.names = FALSE)
  countMatrix[rownames(tmp),i] <- tmp[,7]
}

countMatrix <- as.data.frame(countMatrix)
###
countMatrix$geneID <- geneInfor[rownames(countMatrix),2]
countMatrix <- aggregate(countMatrix[-c(ncol(countMatrix))],by=list(countMatrix$geneID),
                         FUN=sum,na.rm=TRUE)
rownames(countMatrix) <- countMatrix$Group.1
countMatrix <- countMatrix[,-1]
countMatrix <- countMatrix[rowSums(countMatrix>5)>ncol(countMatrix)/2,]

##calculate RPKM
expFPKM <- matrix(0,nrow = nrow(countMatrix),ncol = ncol(countMatrix))
rownames(expFPKM ) <- rownames(countMatrix)
colnames(expFPKM ) <- colnames(countMatrix)

i <- 1
for (i in 1:nrow(expFPKM)) {
  maxLen <- max(geneInfor$Length[which(geneInfor$gene_name %in% rownames(expFPKM)[i])])
  expFPKM[i,] <- as.numeric(countMatrix[i,])*10^9/as.numeric((maxLen*colSums(countMatrix)))
}


sampleInfor <- read.table("GSE161777_sampleInfor.txt",header = TRUE,row.names = 1)

##ROC analysis
finalGene <- read.csv("final.csv")
upGene <- finalGene$geneID[which(finalGene$class=="up")]
downGene <- finalGene$geneID[which(finalGene$class=="down")]

cValue <- c()
rValue <- c()
i <- 1
for (i in 1:nrow(sampleInfor)) {
  upexp <- mean(as.numeric(expFPKM[upGene,sampleInfor$geo[i]]))
  downexp <- mean(as.numeric(expFPKM[downGene,sampleInfor$geo[i]]))
  rValue <- c(rValue,upexp/downexp)
}

sampleInfor$ratio <- rValue
sampleInfor$group <- ""
sampleInfor$group[which(sampleInfor$pseudotime==0)] <- "s1"
sampleInfor$group[which(sampleInfor$pseudotime>0 & sampleInfor$pseudotime<5)] <- "s2"
sampleInfor$group[which(sampleInfor$pseudotime>4)] <- "s3"

ggplot(sampleInfor,aes(x=group,y=ratio,color=group))+
  geom_boxplot(outlier.color = NA,width=0.5)+
  geom_point(pch=20,color=NA)+
  geom_jitter(width = 0.1,size=0.5)+
  scale_y_continuous(breaks = seq(0,60,20),limits = c(0,60))+
  scale_color_nejm()+theme_pubr()+theme(legend.position="none")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5))

