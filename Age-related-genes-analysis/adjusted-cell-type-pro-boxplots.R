#!/usr/bin/Rscript
library(ggplot2)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(stringr)
library(plyr)

##Only the cell types with non-zero estimated proportions in > 50% of samples were retained
cellTypes <- read.table("PBMC_5sigMarker_Results.txt",row.names = 1,header = TRUE,sep = "\t")
cellTypes <- cellTypes[,colSums(cellTypes>0)/nrow(cellTypes)>0.5]
cellTypes <- cellTypes[rowSums(cellTypes>0)/ncol(cellTypes)>0.5,]

##Clinical information
sampleAnnot <- read.table("GTEx_Portal.txt",header = TRUE,sep = "\t",row.names = 1)
sampleAnnot <- unique(sampleAnnot[,c(2:6)])
rownames(sampleAnnot) <- sampleAnnot$Subject.ID

sampleAnnot$Hardy[which(sampleAnnot$Hardy!="Ventilator case")] <- "no"
sampleAnnot$Hardy[which(sampleAnnot$Hardy=="Ventilator case")] <- "yes"
sampleAnnot$AgeBin <- str_replace(sampleAnnot$AgeBin,"-","_")

mydf <- cbind(cellTypes,sampleAnnot[rownames(cellTypes),"AgeBin"])
colnames(mydf)[ncol(mydf)] <- "AgeBin"

##adjusted boxplots
mydf1 <- as.data.frame(matrix(nrow=nrow(mydf),ncol=(ncol(mydf)-1)))
colnames(mydf1) <- colnames(mydf)[1:ncol(mydf1)]
rownames(mydf1) <- rownames(mydf)

sampleAnnot <- sampleAnnot[rownames(mydf),]
colGroup <- data.frame(row.names = rownames(sampleAnnot),sex=factor(sampleAnnot$Sex),hardy=factor(sampleAnnot$Hardy),AgeBin=factor(sampleAnnot$AgeBin))

for (i in 1:(ncol(mydf1))){
  c <- colnames(mydf1)[i]
  mycell <- c
  mydata <- as.matrix(mydf[,colnames(mydf)==mycell])
  cres <- glm(mydata[,1] ~ (sex+hardy),data=colGroup)
  mydf1[,i] <- as.vector(cres$residuals)
}
mydf1$AgeBin <- colGroup$AgeBin

##Plotting
plot_list <- list()
i <- 1
for (i in 1:(ncol(mydf1)-1)){
  mycell <- colnames(mydf1)[i]
  mydata <- mydf1[,c(mycell,"AgeBin")]
  g <- ggplot(data=mydata,aes_string(x="AgeBin",y=mycell))+geom_boxplot(aes(fill=AgeBin),outlier.size = 0.5)+
    scale_fill_nejm()+theme_pubr()+theme(legend.position="none")+ggtitle(label=mycell)+
    theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),plot.title = element_text(hjust = 0.5))
  plot_list[[i]]= g
}

pdf("adjusted-CIBERSORTx-boxplots-agebin5.pdf",height=6,width=6,useDingbats=FALSE)
do.call(grid.arrange,plot_list)
dev.off()
