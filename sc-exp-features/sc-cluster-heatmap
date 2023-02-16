#!/usr/bin/Rscript
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggrepel)
library(tidyr)

##calculate the proportion of a given cell-type expressing a age-related gene
control_expMatrix <- read.csv("control_countExpMatrix.csv",row.names = 1,check.names = FALSE)
control_rawMatrix <- read.csv("control_countRawMatrix.csv",row.names = 1,check.names = FALSE)
moderate_expMatrix <- read.csv("moderate_countExpMatrix.csv",row.names = 1,check.names = FALSE)
moderate_rawMatrix <- read.csv("moderate_countRawMatrix.csv",row.names = 1,check.names = FALSE)
severe_expMatrix <- read.csv("severe_countExpMatrix.csv",row.names = 1,check.names = FALSE)
severe_rawMatrix <- read.csv("severe_countRawMatrix.csv",row.names = 1,check.names = FALSE)

filter_process <- function(dataMatrix,n){
  i <- 1
  retain_cols<- c()
  for (i in 1:ncol(dataMatrix)) {
    lines <- which(dataMatrix[,i]>n*20)
    if(length(lines)>nrow(dataMatrix)/2){
      retain_cols <- c(retain_cols,colnames(dataMatrix)[i])
    }
  }
  
  i <- 1
  retain_rows<- c()
  for (i in 1:nrow(dataMatrix)) {
    lines <- which(dataMatrix[i,]>n*20)
    if(length(lines)>0){
      retain_rows <- c(retain_rows,rownames(dataMatrix)[i])
    }
  }
  
  return(dataMatrix[retain_rows,retain_cols])
  
}

control_f <- filter_process(control_rawMatrix,15)
moderate_f <- filter_process(moderate_rawMatrix,15)
severe_f <- filter_process(severe_rawMatrix,12)

cellTypes <- sort(unique(c(colnames(control_f),colnames(moderate_f),colnames(severe_f))))
allGene <- sort(unique(c(rownames(control_f),rownames(moderate_f),rownames(severe_f))))

count_ratio <- matrix(0,nrow = length(allGene),ncol = length(cellTypes))
rownames(count_ratio) <- allGene
colnames(count_ratio) <- cellTypes

calculate_ratio <- function(count_ratio,rawMatrix,expMatrix){
  i <- 1
  for (i in 1:nrow(count_ratio)) {
    k <- 1
    for (k in 1:ncol(count_ratio)) {
      row_num <- which(rownames(rawMatrix)==rownames(count_ratio)[i])
      col_num <- which(colnames(rawMatrix)==colnames(count_ratio)[k])
      if(length(row_num)>0 & length(col_num)>0){
        if(rawMatrix[row_num,col_num]>0){
          count_ratio[i,k] <- expMatrix[row_num,col_num]/rawMatrix[row_num,col_num]*100
        }
      }
    }
  }
  
  return(count_ratio)
}

control_ratio <- calculate_ratio(count_ratio,control_rawMatrix,control_expMatrix)
moderate_ratio <- calculate_ratio(count_ratio,moderate_rawMatrix,moderate_expMatrix)
severe_ratio <- calculate_ratio(count_ratio,severe_rawMatrix,severe_expMatrix)

##Cluster of genes based on expression in some cell-types
diffResult <- read.table("degreport-p0.05-cov-genes-table.upC.downC.txt",row.names = 1,header = TRUE)
upGene <- rownames(diffResult)[which(diffResult$cluster==1)]
downGene <- rownames(diffResult)[which(diffResult$cluster==2)]

control_ratio_up <- control_ratio[intersect(rownames(control_ratio),upGene),]
control_ratio_down <- control_ratio[intersect(rownames(control_ratio),downGene),]
moderate_ratio_up <- moderate_ratio[intersect(rownames(moderate_ratio),upGene),]
moderate_ratio_down <- moderate_ratio[intersect(rownames(moderate_ratio),downGene),]
severe_ratio_up <- severe_ratio[intersect(rownames(severe_ratio),upGene),]
severe_ratio_down <- severe_ratio[intersect(rownames(severe_ratio),downGene),]

bk <- seq(0,100,by=1)
plotData <- pheatmap(t(control_ratio_up),scale = "none",show_rownames = TRUE,
                     show_colnames = FALSE,cluster_rows = FALSE,cluster_cols = TRUE,
                     treeheight_row = 0,treeheight_col = 0,border_color = NA,na_col = "white",
                     color = colorRampPalette(brewer.pal(n = 8, name ="RdPu"))(101),
                     legend_breaks=seq(0,100,20),breaks=bk,cutree_cols = 4)

cluster_cols <- rownames(control_ratio_up)[plotData$tree_col$order]
cluster_rows <- colnames(control_ratio_up)
cols_cluster <- cutree(plotData$tree_col,k=4)

geneSet <- names(cols_cluster)[which(cols_cluster==2)]
write.table(geneSet,file = "geneSet.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)

bk <- seq(0,100,by=1)
pheatmap(t(con_moderate_down[cluster_cols,cluster_rows]),scale = "none",show_rownames = TRUE,
         show_colnames = FALSE,cluster_rows = FALSE,cluster_cols = FALSE,
         treeheight_row = 0,treeheight_col = 0,border_color = NA,na_col = "white",
         color = colorRampPalette(brewer.pal(n = 8, name ="RdPu"))(101),
         legend_breaks=seq(0,100,20),breaks=bk)
