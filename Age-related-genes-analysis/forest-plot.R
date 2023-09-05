#!/usr/bin/Rscript
library(MASS)
library(ordinal)
library(ggplot2)
library(ggpubr)
library(forcats)

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

sampleAnnot <- sampleAnnot[rownames(cellTypes),]
colGroup <- data.frame(row.names = rownames(sampleAnnot),sex=factor(sampleAnnot$Sex),age=sampleAnnot$Age,hardy=factor(sampleAnnot$Hardy))

##ordinal logistic regression model
##Cell-type proportions were adjusted using residuals from a generalized linear model
mydf <- matrix(nrow=ncol(cellTypes),ncol=5)
i <- 1
for (i in 1:ncol(cellTypes)){
  c <- colnames(cellTypes)[i]
  myCell <- c
  mydata <- as.numeric(cellTypes[,colnames(cellTypes)==myCell])
  mydataO <- ordered(mydata,levels=c(sort(unique(mydata))))
  
  cres <- clm(mydataO ~ (age+sex+hardy),data=colGroup,link="logit",Hess=TRUE)
  csum <- summary(cres)$coef
  mydf[i,1] <- c
  mydf[i,2:5] <- csum["age",]
}
colnames(mydf) <- c("Cell","OLR_estimate","StdError","z","pval")

mydf <- as.data.frame(mydf)
mydf$OLR_estimate <- as.numeric(as.character(mydf$OLR_estimate))
mydf$StdError <- as.numeric(as.character(mydf$StdError))
mydf$pval <- as.numeric(as.character(mydf$pval))

mydf$lowerCI <- mydf$OLR_estimate - 1.96*mydf$StdError
mydf$upperCI <- mydf$OLR_estimate + 1.96*mydf$StdError
mydf$adjp <- p.adjust(mydf$pval,method="BH")

rownames(mydf) <- mydf$Cell
mydf <- mydf[order(mydf$adjp),]
mydf$nlp <- log10(mydf$adjp)*-1
updata <- mydf[mydf$OLR_estimate > 0 & mydf$adjp < 0.05,]
downdata <- mydf[mydf$OLR_estimate < 0 & mydf$adjp < 0.05,]
otherdata <- mydf[mydf$OLR_estimate ==0 | mydf$adjp > 0.05,]
mydf$color <- rep("black",nrow(mydf))

mydf[mydf$Cell %in% updata$Cell,"color"] <- "#00AFBB"
mydf[mydf$Cell %in% downdata$Cell,"color"] <- "#FC4E07"

mydf <- mydf[nrow(mydf):1,]
mydf$Cell <- fct_inorder(mydf$Cell)

pdf("OLR-CIBERSORTx-blood.age5.pdf",height=2.5,width=4,useDingbats=FALSE)
ggplot(mydf,aes(x=OLR_estimate,y=Cell))+geom_vline(xintercept=0,linetype=2,color="gray44")+
  geom_errorbar(aes(xmin=lowerCI,xmax=upperCI),width=0.2)+geom_point(aes(color=color,size=nlp))+
  theme_bw()+scale_color_manual(values=c("#FC4E07","#00AFBB","black")) + 
  scale_size_continuous(name = "-log pvalue",range=c(1,4))
dev.off()

