#!/usr/bin/Rscript
library(DESeq2)
library(stringr)
library(limma)
library(DEGreport)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(ggsci)

gtexExp <- read.table("gene_reads_whole_blood.gct",header = TRUE,row.names = 1,check.names = FALSE)##Expression matrix
gtexExp <- aggregate(gtexExp[-c(1,2)],by=list(gtexExp$Description),FUN=sum,na.rm=TRUE)
rownames(gtexExp) <- gtexExp$Group.1
gtexExp <- gtexExp[,-1]

i <- 1
for (i in 1:ncol(gtexExp)) {
  tmp <- unlist(strsplit(colnames(gtexExp)[i],"-"))
  colnames(gtexExp)[i] <- paste(tmp[1],tmp[2],sep = "-")
}## change the sample name

sampleAnnot <- read.table("GTEx_Portal.txt",header = TRUE,sep = "\t",row.names = 1) ##Clinical information
sampleAnnot <- unique(sampleAnnot[,c(2:6)])
rownames(sampleAnnot) <- sampleAnnot$Subject.ID
sampleAnnot <- sampleAnnot[colnames(gtexExp),]

sampleAnnot$Hardy[which(sampleAnnot$Hardy!="Ventilator case")] <- "no"
sampleAnnot$Hardy[which(sampleAnnot$Hardy=="Ventilator case")] <- "yes"
sampleAnnot$AgeBin <- str_replace(sampleAnnot$AgeBin,"-","_")##"20-29" "30-39" "40-49" "50-59" "60-69" "70-79"

geneAnnot <- read.table("geneInfo.tab",row.names = 1)
geneAnnot <- geneAnnot[which(geneAnnot$V3=="protein_coding"),]
gtexExp <- gtexExp[which(rownames(gtexExp) %in% geneAnnot$V2),rownames(sampleAnnot)]

##Differential expression analysis
colGroup <- data.frame(row.names = rownames(sampleAnnot),sex=factor(sampleAnnot$Sex),
                       age=factor(sampleAnnot$AgeBin),
                       hardy=factor(sampleAnnot$Hardy))

dds <- DESeqDataSetFromMatrix(countData = gtexExp, colData = colGroup, 
                              design= ~sex + hardy + age,tidy=FALSE)
dds <- dds[rowSums(counts(dds) >= 5) >= ncol(gtexExp)/2, ]
dds <- DESeq(dds,test = "LRT",reduced= ~sex+hardy)
saveRDS(dds,"DESeq2.LRT.agebin.rds")

dds <- readRDS("DESeq2.LRT.agebin.rds")
res <- results(dds,contrast=c("age","70_79","60_69"),independentFiltering=TRUE,alpha=0.05,pAdjustMethod="BH")
res <- res[order(res$padj),]
head(res)
summary(res)
write.table(res,"GTEx-blood-DEseq2-results.LRT.timeBin.7-6.Sexhardy.txt",sep="\t",row.names=TRUE,quote = FALSE)

##"batch correction" for downstream visualization of adjusted expression values
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
plotPCA(vsd,"sex")
plotPCA(vsd,"hardy")
plotPCA(vsd,"age")

assay(vsd) <- removeBatchEffect(assay(vsd),covariates=model.matrix(~sex+hardy,data=colGroup),design = model.matrix(~age,data=colGroup))
write.table(format(as.data.frame(assay(vsd)),digits=3),"GTEX-vst.adjusted.txt",sep="\t",row.names=TRUE,quote = FALSE)

##Gene clusters
adjCount <- read.table("GTEX-vst.adjusted.txt",row.names=1,check.names = FALSE)
difResult <- read.table("GTEx-blood-DEseq2-results.LRT.timeBin.combine.Sexhardy.txt",sep="\t",header=TRUE,check.names = FALSE,row.names = 1)
difResult <- difResult[difResult$padj<0.05,]
adjCount <- adjCount[rownames(adjCount) %in% rownames(difResult),]

clusters <- degPatterns(adjCount, metadata = colGroup, time = "age")
filtclusters <- clusters$normalized
upc <- c("5")
downc <- c("6")
filtclusters <- filtclusters[filtclusters$cluster %in% upc | filtclusters$cluster %in% downc,]
filtclusters[filtclusters$cluster %in% upc,"cluster"] = 1
filtclusters[filtclusters$cluster %in% downc,"cluster"] = 2

pdf("degreport-p0.05-cov-genes-color-boxplot.lines.pdf",height=4,width=9,useDingbats=FALSE)
degPlotCluster(filtclusters,time="age",color="age",points = FALSE,boxes = FALSE)+
  geom_boxplot(outlier.color = NA,fill=NA,width=0.6)+
  geom_point(pch=20,color=NA)+
  geom_jitter(width = 0.06,size=0.5)+
  geom_line(aes_string(group="genes"),alpha=0.05)+
  scale_color_nejm()+theme_bw()
dev.off()

mydf <- clusters$df
write.table(mydf,"degreport-p0.05-cov-genes-table.all.txt",sep="\t",row.names=FALSE,quote = FALSE)

mydf3 <- mydf[mydf$cluster %in% upc | mydf$cluster %in% downc,]
mydf3[mydf3$cluster %in% upc,2] = "1"
mydf3[mydf3$cluster %in% downc,2] = "2"

write.table(mydf3,"degreport-p0.05-cov-genes-table.upC.downC.txt",sep="\t",row.names=FALSE,quote = FALSE)
