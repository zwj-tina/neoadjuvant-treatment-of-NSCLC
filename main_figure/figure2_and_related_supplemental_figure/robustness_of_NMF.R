library(NMF)
library(ComplexHeatmap)
library(reshape2)
library(tidyverse)
library(dplyr)
library(readxl)
library(viridis)
###################################################################
count <- 1
module_merge <- list()

for(i in 1:200){
  info <- read.csv("all_sub_cell_type.csv")
  group <- read.csv("NMF_all_group_5.csv")
  
  #remove 20% samples randomly
  random.samples <- sample(group$sampleID, 45)
  print(random.samples)
  group <- group[!group$sampleID %in% random.samples,]
  info <- info[info$sampleID %in% group$sampleID,]
  print(length(unique(info$sampleID)))
  
  df <- table(info$sampleID,info$sub_cell_type)
  ratio <- as.data.frame(df / rowSums(df))
  colnames(ratio) <- c("sampleID","cell.type","Freq")
  print(length(unique(ratio$sampleID)))
  
  ratio <- dcast(ratio, sampleID ~ ratio$cell.type, value.var = "Freq")
  rownames(ratio) <- ratio$sampleID
  
  ratio <- ratio[,-1]
  head(ratio)
  ratio[is.na(ratio)] <- 0
  
  #normalization
  scale_ratio <- apply(ratio, MARGIN = 2, function(x) (x-min(x))/(max(x)-min(x)))
  scale_ratio <- as.data.frame(scale_ratio)
  scale_ratio <- t(scale_ratio)
  
  seed = 2020820
  for(rk in 2:10){
    nmf.rank5 <- nmf(scale_ratio, 
                     rank = rk, 
                     nrun=200,
                     seed = seed, 
                     method = "lee")
    
    index <- extractFeatures(nmf.rank5,"max") 
    for(j in 1:rk){
      part <- scale_ratio[index[[j]],]
      module_merge[[count]] <- rownames(part)
      count <- count + 1
    }
  }
}

module_merge


Mat <- matrix(0, ncol = length(unique(info$sub_cell_type)), nrow = length(unique(info$sub_cell_type)))
rownames(Mat) <- unique(info$sub_cell_type)
colnames(Mat) <- unique(info$sub_cell_type)
head(Mat)

for (i in 1:length(unique(info$sub_cell_type))) {
  for (j in 1:length(unique(info$sub_cell_type))) {
    number <- 0
    for(m in 1:length(module_merge)){
      if((rownames(Mat)[i] %in% module_merge[[m]]) & (rownames(Mat)[j] %in% module_merge[[m]])){
        number <- number + 1
      }
    }
    Mat[i,j] <- number
  }
}

custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
pheatmap(as.matrix(Mat), cluster_cols=T, cluster_rows=T, 
         clustering_distance_rows="euclidean", color=custom_magma, 
         fontsize=12,treeheight_row=0,treeheight_col=30, 
         cellheight = 7,cellwidth = 7,show_rownames=T, 
         show_colnames=F,clustering_method = "ward.D2", border_color = NA)

