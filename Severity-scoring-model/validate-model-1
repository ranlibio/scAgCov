#!/usr/bin/Rscript
library(ggplot2)
library(tidyr)
library(ggsci)
library(ggpubr)
library(ggrepel)
library(forcats)

##Inputs
sampleInfor <- read.table("sampleInfor.txt",row.names = 1,header = TRUE,sep = "\t",check.names = FALSE)

psuedoCount <- read.table("psuedo-bulk.txt",header = TRUE,row.names = 1,check.names = FALSE)
psuedoCount <- psuedoCount[,rownames(sampleInfor)]
psuedoCount <- psuedoCount[rowSums(psuedoCount>5)>ncol(psuedoCount)/2,]
psuedoCount <- psuedoCount[,which(colSums(psuedoCount)>10*10^6)]

##calculate RPKM
geneLen <- read.table("D:/liran_paper/covid19/data/trainModel/gencode_v35.txt",row.names = 1)
expFPKM <- matrix(0,nrow = nrow(psuedoCount),ncol = ncol(psuedoCount))
rownames(expFPKM) <- rownames(psuedoCount)
colnames(expFPKM) <- colnames(psuedoCount)

i <- 1
for (i in 1:nrow(expFPKM)) {
  maxLen <- max(geneLen$V4[which(geneLen$V2 %in% rownames(expFPKM)[i])])
  expFPKM[i,] <- as.numeric(psuedoCount[i,])*10^9/as.numeric((maxLen*colSums(psuedoCount)))
}

expFPKM <- expFPKM[which(!is.na(expFPKM[,1])),]

##ROC analysis
finalGene <- read.csv("final.csv")
upGene <- finalGene$geneID[which(finalGene$class=="up")]
downGene <- finalGene$geneID[which(finalGene$class=="down")]

tmp_infor <- sampleInfor[colnames(expFPKM),]
cValue <- c()
rValue <- c()
i <- 1
for (i in 1:nrow(tmp_infor)) {
  upexp <- mean(as.numeric(expFPKM[upGene,rownames(tmp_infor)[i]]))
  downexp <- mean(as.numeric(expFPKM[downGene,rownames(tmp_infor)[i]]))
  rValue <- c(rValue,upexp/downexp)
}

cValue <- as.character(tmp_infor$Stage)
ROCdata <- data.frame(outcome=cValue,value=rValue)
rownames(ROCdata) <- rownames(tmp_infor)
ROCdata$time <- sampleInfor[rownames(ROCdata),"Samplingday"]

##Boxplot
ggplot(ROCdata,aes(x=outcome,y=value,color=outcome))+
  geom_boxplot(outlier.color = NA,width=0.5)+
  geom_point(pch=20,color=NA)+
  geom_jitter(width = 0.1,size=0.5)+
  scale_y_continuous(breaks = seq(0,100,20),limits = c(0,100))+
  scale_color_nejm()+theme_pubr()+theme(legend.position="none")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5))

##Correlation
corData <- ROCdata[which(ROCdata$outcome=="mc" | ROCdata$outcome=="sc"),]
corData$time <- as.numeric(corData$time)
ggplot(corData,aes(x=time,y=value))+geom_point(color="#808080",cex=0.5)+
  scale_x_continuous(breaks = seq(0,150,30),limits = c(0,140))+
  scale_y_continuous(breaks = seq(0,50,10),limits = c(0,50))+
  theme_classic() +theme_pubr()+ theme(legend.position="none")+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)

##plot ROC curve
library(pROC)
library(plotROC)
ROCcurve <- roc(response=ROCdata$outcome,predictor=ROCdata$value,levels = c("sp","mp"),
                ci=TRUE,SP=TRUE,SE=TRUE,plot = TRUE,percent=TRUE,auc = TRUE)

plot.roc(ROCcurve, print.auc = TRUE, print.thres = "best")

