#!/usr/bin/Rscript
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(ggsci)
library(ggpubr)
library(ggrepel)
library(Seurat)
library(Matrix)

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

##selecte genes whose expression was progressively increased or decreased from healthy controls to mild-moderate to severe patients
diffResult <- read.table("degreport-p0.05-cov-genes-table.upC.downC.txt",row.names = 1,header = TRUE)

upGene <- rownames(diffResult)[which(diffResult$cluster==1)]
downGene <- rownames(diffResult)[which(diffResult$cluster==2)]

control_ratio_up <- control_ratio[intersect(rownames(control_ratio),upGene),]
control_ratio_down <- control_ratio[intersect(rownames(control_ratio),downGene),]
moderate_ratio_up <- moderate_ratio[intersect(rownames(moderate_ratio),upGene),]
moderate_ratio_down <- moderate_ratio[intersect(rownames(moderate_ratio),downGene),]
severe_ratio_up <- severe_ratio[intersect(rownames(severe_ratio),upGene),]
severe_ratio_down <- severe_ratio[intersect(rownames(severe_ratio),downGene),]

selected_genes <- matrix(0,nrow = nrow(control_ratio_up),ncol = ncol(control_ratio_up))
rownames(selected_genes) <- rownames(control_ratio_up)
colnames(selected_genes) <- colnames(control_ratio_up)

i <- 1
for (i in 1:nrow(selected_genes)) {
  k <- 1
  for (k in 1:ncol(selected_genes)) {
    a1 <- control_ratio[rownames(selected_genes)[i],colnames(selected_genes)[k]]
    a2 <- moderate_ratio[rownames(selected_genes)[i],colnames(selected_genes)[k]]
    a3 <- severe_ratio[rownames(selected_genes)[i],colnames(selected_genes)[k]]
    if(max(a1,a3)>40){
      tmp <- (a3+10)/(a1+10)
      if((tmp>1.5) | (tmp<0.67))
      {
        if(a1>a2 & a2>a3){
          selected_genes[i,k] <- (-1)
        }
        if(a1<a2 & a2<a3){
          selected_genes[i,k] <- 1
        }
      }
    }
  }
}

selected_genes <- selected_genes[which(rowSums(abs(selected_genes))>0),which(colSums(abs(selected_genes))>0)]

##cell-types with significant difference in proportion of cells expressing certain age-pos genes among patients with different severity during progression stages
sampleInfor <- read.table("sample_infor.txt",row.names = 1,header = TRUE,sep = "\t",check.names = FALSE)

stage_1 <- rownames(sampleInfor)[which(sampleInfor$CoVID_severity=="control" & sampleInfor$Sample_time=="control")]
stage_2 <- rownames(sampleInfor)[which(sampleInfor$CoVID_severity=="mild/moderate" & sampleInfor$Sample_time=="progression")]
stage_3 <- rownames(sampleInfor)[which(sampleInfor$CoVID_severity=="severe/critical" & sampleInfor$Sample_time=="progression")]

files <- c(stage_1,stage_2,stage_3)

sigMatrix <- matrix(1,nrow = nrow(selected_genes),ncol = ncol(selected_genes))
rownames(sigMatrix) <- rownames(selected_genes)
colnames(sigMatrix) <- colnames(selected_genes)

i <- 1
for (i in 1:nrow(sigMatrix)) {
  j <- 1
  for (j in 1:ncol(sigMatrix)) {
    
    tmp_ratio <- c()
    tmp_class <- c()
    
    k <- 1
    for (k in 1:length(files)) {
      
      tmp_file <- read.csv(paste(files[k],".csv",sep = ""),row.names = 1,header = TRUE,check.names = FALSE)
      if(is.na(tmp_file[rownames(sigMatrix)[i],colnames(sigMatrix)[j]])){
        tmp_ratio <- c(tmp_ratio,0)
      }else{
        tmp_ratio <- c(tmp_ratio,tmp_file[rownames(sigMatrix)[i],colnames(sigMatrix)[j]])
      }
      
      if(files[k] %in% stage_1){
        tmp_class <- c(tmp_class,"stage1")
      }
      if(files[k] %in% stage_2){
        tmp_class <- c(tmp_class,"stage2")
      }
      if(files[k] %in% stage_3){
        tmp_class <- c(tmp_class,"stage3")
      }
    }
    
    aovData <- data.frame(ratio=tmp_ratio,class=tmp_class)
    kt_tmp <- kruskal.test(ratio~class, aovData) 
    sigMatrix[i,j] <- kt_tmp$p.value
  }
}

sigMatrix.padjust <- matrix(p.adjust(c(sigMatrix),method = "BH"),nrow=nrow(sigMatrix),ncol=ncol(sigMatrix))
colnames(sigMatrix.padjust) <- colnames(sigMatrix)
rownames(sigMatrix.padjust) <- rownames(sigMatrix)

##plotting
sigMatrix.padjust <- sigMatrix.padjust*abs(selected_genes)
sigMatrix.padjust <- sigMatrix.padjust[rowSums(sigMatrix.padjust>0 & sigMatrix.padjust<0.05)>0,]
sigMatrix.padjust <- sigMatrix.padjust[,colSums(sigMatrix.padjust>0 & sigMatrix.padjust<0.05)>0]

i <- 1
for (i in 1:nrow(sigMatrix.padjust)) {
  sigMatrix.padjust[i,which(as.numeric(sigMatrix.padjust[i,])==0 | 
                              as.numeric(sigMatrix.padjust[i,])>=0.05)] <- ""
  sigMatrix.padjust[i,which(as.numeric(sigMatrix.padjust[i,])<0.001)] <- "***"
  sigMatrix.padjust[i,which(as.numeric(sigMatrix.padjust[i,])<0.01)] <- "**"
  sigMatrix.padjust[i,which(as.numeric(sigMatrix.padjust[i,])<0.05)] <- "*"
}

pheatmap(t(selected_genes[rownames(sigMatrix.padjust),colnames(sigMatrix.padjust)]),
         scale = "none",cluster_cols = TRUE,treeheight_col = 0,number_color="white",
         cluster_rows=FALSE, display_numbers=t(sigMatrix.padjust),
         border_color = "white",
         color = c("#0571B0","white","#CA0020"))

