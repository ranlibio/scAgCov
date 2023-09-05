#!/usr/bin/Rscript
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(ggsci)
library(ggpubr)
library(ggrepel)
library(GSVA)
library(forcats)

##calculate RPKM
sampleInfor <- read.table("sampleInfor.txt",row.names = 1,header = TRUE,sep = "\t",check.names = FALSE)
sampleInfor <- sampleInfor[which(sampleInfor$Sample_type=="fresh PBMC" | sampleInfor$Sample_type=="frozen PBMC"),]

geneLen <- read.table("gencode_v35.txt",row.names = 1)
psuedoCount <- read.table("psuedo-bulk.txt",header = TRUE,row.names = 1,check.names = FALSE)
psuedoCount <- psuedoCount[rowSums(psuedoCount>5)>ncol(psuedoCount)/2,rownames(sampleInfor)]
psuedoCount <- psuedoCount[,-which(colSums(psuedoCount)<20*10^6)]
psuedoCount <- psuedoCount[intersect(rownames(psuedoCount),unique(geneLen$V2)),]
expFPKM <- matrix(0,nrow = nrow(psuedoCount),ncol = ncol(psuedoCount))
rownames(expFPKM) <- rownames(psuedoCount)
colnames(expFPKM) <- colnames(psuedoCount)

i <- 1
for (i in 1:nrow(expFPKM)) {
  maxLen <- max(geneLen$V4[which(geneLen$V2 %in% rownames(expFPKM)[i])])
  expFPKM[i,] <- as.numeric(psuedoCount[i,])*10^9/as.numeric((maxLen*colSums(psuedoCount)))
}

##calculate ssGSEA score
ageGene <- read.table("degreport-p0.05-cov-genes-table.upC.downC.txt",row.names = 1,header = TRUE)
gomt <- data.frame(term=rep(c("up","down"),time=c(479,455)),gene=c(rownames(ageGene)[which(ageGene$cluster==1)],rownames(ageGene)[which(ageGene$cluster==2)]))
gomt <- split(gomt$gene,gomt$term)
ssgsea <- gsva(as.matrix(expFPKM),gomt,method="ssgsea")
ssgsea <- ssgsea[,intersect(colnames(ssgsea),rownames(sampleInfor))]
ssgsea <- ssgsea[c(2,1),]

##adjusted ssGSEA score
ssgsea_tmp <- ssgsea
colGroup <- data.frame(row.names = rownames(sampleInfor),Sex=factor(sampleInfor$Sex),
                       Sampletype=factor(sampleInfor$Sampletype),Stage=factor(sampleInfor$Stage),
                       CoVIDseverity=factor(sampleInfor$CoVIDseverity),AgeBin=factor(sampleInfor$AgeBin),
                       ST=factor(sampleInfor$ST))

tt <- removeBatchEffect(ssgsea_tmp,covariates=model.matrix(~Sex+Sampletype,data=colGroup),
                        design = model.matrix(~Stage,data=colGroup))

plotData <- data.frame(x=sampleInfor[colnames(tt),"Stage"],y=tt[2,])
plotData$y <- scale(plotData$y)
plotData$x <- fct_inorder(plotData$x)

ggplot(data=plotData,aes(x=x,y=y,color=x))+
  geom_boxplot(outlier.color = NA,width=0.5)+
  geom_point(pch=20,color=NA)+
  geom_jitter(width = 0.1,size=0.5)+
  scale_y_continuous(breaks = seq(-3,3,3),limits = c(-3,3))+
  scale_color_nejm()+theme_pubr()+theme(legend.position="none")+stat_compare_means(method = "kruskal.test")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5))

