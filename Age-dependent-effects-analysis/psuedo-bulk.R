#!/usr/bin/Rscript
library(Seurat)
library(Matrix)

##psuedo-bulk
geneID <- read.table("GSE158055_covid19_features.tsv")
files <- list.files("./")

rawCount <- matrix(0,nrow = length(geneID$V1),ncol = length(files))
rownames(rawCount) <- geneID$V1
colnames(rawCount) <- files

i <- 1
for (i in 1:length(files)) {
  
  print(i)
  print(Sys.time())
  scData <- Read10X(files[i],gene.column = 1)
  scData <- CreateSeuratObject(counts = scData,min.cells = 3, min.features = 200)
  scCount <- GetAssayData(object = scData, slot = "counts")
  rawCount[rownames(scCount),files[i]] <- rowSums(scCount)
}

write.table(rawCount,file = "psuedo-bulk.txt",row.names = TRUE,quote = FALSE,col.names = TRUE,sep = "\t")
