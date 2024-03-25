library(NMF)
library(ComplexHeatmap)
library(reshape2)
library(tidyverse)
library(dplyr)
library(readxl)
###################################################################
info <- read.csv("all_sub_cell_type.csv")
head(info)
length(unique(info$sampleID))

df <- table(info$sampleID,info$sub_cell_type)
ratio <- as.data.frame(df / rowSums(df))
head(ratio)
colnames(ratio) <- c("sampleID","cell.type","Freq")
head(ratio)
length(unique(ratio$sampleID))

sample.info <- as.data.frame(read_excel("sample.xlsx"))
head(sample.info)
sample.info <- sample.info[sample.info$sampleID %in% LUAD.group$sampleID,]


pathological_response_level <- c()
for(each in sample.info$pathological_response){
  if(each %in% c("MPR","pCR")){
    pathological_response_level <- c(pathological_response_level, "MPR")
  }else{
    pathological_response_level <- c(pathological_response_level, "non-MPR")
  }
}

sample.info$pathological_response_level <- pathological_response_level

response.meta <- sample.info[, c("sampleID","smoking_history","cancer_type","pre_treatment_staging",
                                 "PDL1_TPS","PD1","chemotherapy","targeted_therapy","cycles",
                                 "pathological_response","pathological_response_level",
                                 "pathological_response_rate","radiological_response",
                                 "RVT_pre_dominant_histology")]
response.meta <- response.meta %>% distinct(sampleID, .keep_all = TRUE)
head(response.meta)
length(unique(response.meta$sampleID))

ratio <- dcast(ratio, sampleID ~ ratio$cell.type, value.var = "Freq")


merge.version <- merge(ratio, response.meta, by = "sampleID", all.x = TRUE)
head(merge.version)
dim(merge.version)

rownames(ratio) <- ratio$sampleID
ratio <- ratio[merge.version$sampleID,]
dim(ratio)

ratio <- ratio[,-1]
head(ratio)
ratio[is.na(ratio)] <- 0


#normalization
scale_ratio <- apply(ratio, MARGIN = 2, function(x) (x-min(x))/(max(x)-min(x)))
head(scale_ratio)
scale_ratio <- as.data.frame(scale_ratio)
head(scale_ratio)
scale_ratio <- t(scale_ratio)
head(scale_ratio)
dim(scale_ratio)

#ranks <- 2:10
#estim.coad <- nmf(scale_ratio, ranks, nrun=100)
#plot(estim.coad)

#consensusmap(estim.coad)


#NMF,rank=5
seed = 2020820

nmf.rank5 <- nmf(scale_ratio, 
                 rank = 5, 
                 nrun=200,
                 seed = seed, 
                 method = "brunet")

Heatmap(as.matrix(nmf.rank5@consensus))

index <- extractFeatures(nmf.rank5,"max") 
#index[[2]] <- c(15,17,24,25)
#index[[4]] <- c(43,46,20,47)

new.index <- list()
new.index[[1]] <- index[[4]]
new.index[[2]] <- index[[5]]
new.index[[3]] <- index[[1]]
new.index[[4]] <- index[[2]]
new.index[[5]] <- index[[3]]

sig.order <- unlist(new.index)
NMF.Exp.rank5 <- scale_ratio[sig.order,]
NMF.Exp.rank5 <- na.omit(NMF.Exp.rank5) 
dim(NMF.Exp.rank5)

group <- predict(nmf.rank5) 

#adjust the position of the module
new.group <- c()
for(each in group){
  if(each %in% c("4")){
    new.group <- c(new.group, "1")
  }
  if(each %in% c("5")){
    new.group <- c(new.group, "2")
  }
  if(each %in% c("1")){
    new.group <- c(new.group, "3")
  }
  if(each %in% c("2")){
    new.group <- c(new.group, "4")
  }
  if(each %in% c("3")){
    new.group <- c(new.group, "5")
  }
}
new.group <- factor(new.group, levels = c("1","2","3","4","5"))

#plot heatmap
z_ratio <- scale(ratio)/4
head(z_ratio)
z_ratio <- as.data.frame(z_ratio)
head(z_ratio)
z_ratio <- t(z_ratio)

plot_matrix <- z_ratio[sig.order,]
plot_matrix <- na.omit(plot_matrix)
dim(plot_matrix)


info.matrix <- as.data.frame(t(NMF.Exp.rank5))
head(info.matrix)
info.matrix$sampleID <- rownames(info.matrix)
info.matrix$group <- new.group
info.matrix <- merge(info.matrix, response.meta, by = "sampleID", all.x = TRUE)
head(info.matrix)

gene.group <- c()
for(each in rownames(NMF.Exp.rank5)){
  if(each %in% rownames(scale_ratio)[new.index[[1]]]){
    gene.group <- c(gene.group, "module1")
  }else if(each %in% rownames(scale_ratio)[new.index[[2]]]){
    gene.group <- c(gene.group, "module2")
  }else if(each %in% rownames(scale_ratio)[new.index[[3]]]){
    gene.group <- c(gene.group, "module3")
  }else if(each %in% rownames(scale_ratio)[new.index[[4]]]){
    gene.group <- c(gene.group, "module4")
  }else if(each %in% rownames(scale_ratio)[new.index[[5]]]){
    gene.group <- c(gene.group, "module5")
  }
}

PDL1_TPS_group <- c()
for(each in info.matrix$PDL1_TPS){
  if (each %in% c("<1%","0")){
    PDL1_TPS_group <- c(PDL1_TPS_group,"<1%") 
  }else if (each %in% c("0.01","0.02","0.05","0.08","0.03","0.3","0.2","0.35","0.4")){
    PDL1_TPS_group <- c(PDL1_TPS_group,"1%-49%") 
  }else if (each %in% c("0.7","0.6","0.9","0.8","0.85","1","0.55","0.75","0.65","0.5")){
    PDL1_TPS_group <- c(PDL1_TPS_group,">=50%") 
  }else{
    PDL1_TPS_group <- c(PDL1_TPS_group,"Not tested") 
  }
}
info.matrix$PDL1_TPS_group <- PDL1_TPS_group


RVT_pre_dominant_histology <- c()
for(each in info.matrix$RVT_pre_dominant_histology){
  if (each %in% c("complex acinar","MIP","acinar","papillary")){
    RVT_pre_dominant_histology <- c(RVT_pre_dominant_histology,"non-solid") 
  }else if (each %in% c("solid")){
    RVT_pre_dominant_histology <- c(RVT_pre_dominant_histology,"solid") 
  }else{
    RVT_pre_dominant_histology <- c(RVT_pre_dominant_histology,"unknown") 
  }
}
info.matrix$RVT_pre_dominant_histology <- RVT_pre_dominant_histology


info.matrix[is.na(info.matrix)] <- "unknown"
ha = HeatmapAnnotation(smokingHistory = factor(info.matrix$smoking_history, levels = c("Y","N","unknown")),
                       cycles = factor(info.matrix$cycles, levels = c("unknown","2","3","4","5","6")),
                       PDL1_TPS = factor(info.matrix$PDL1_TPS_group, levels = c("Not tested","<1%","1%-49%",">=50%")),
                       RVT_pre_dominant_histology = factor(info.matrix$RVT_pre_dominant_histology, levels = c("unknown","non-solid","solid")),
                       pathologicalResponse = factor(info.matrix$pathological_response, levels = c("non-MPR","nPR","pPR","MPR","pCR")),
                       pathologicalResponseLevel = factor(info.matrix$pathological_response_level, levels = c("MPR","non-MPR")),
                       group = factor(info.matrix$group, levels = c("1","2","3","4","5")),
                       col = list(smokingHistory = c("Y" = "#F6E382","N" = "#B8DCC5","unknown" = "#82B0D2"),
                                  cycles = c("unknown" = "#E8E8D0","2" = "#DEDEBE","3" = "#CDCD9A","4"="#B9B973","5"="#AFAF61","6"="#949449"),
                                  pathologicalResponseLevel = c("MPR" = "#2868A6", "non-MPR" = "#B1161C"),
                                  pathologicalResponse = c("nPR" = "#EDB752","pPR" = "#FFC1B3","MPR" = "#AFD483","pCR" = "#469832","non-MPR" = "#DFE0DF"),
                                  PDL1_TPS = c("Not tested" = "#D1E9E9","<1%"="#B3D9D9","1%-49%"="#6FB7B7",
                                               ">=50%"="#4F9D9D"),
                                  RVT_pre_dominant_histology = c("non-solid" = "#63BAAB",
                                                                 "solid" = "#D87635",
                                                                 "unknown" = "#C9CACB"),
                                  group = c("1"="#E84C35","2"="#4FBAD6",
                                            "3"="#00A289","4"="#3C5487",
                                            "5"="#F29B80")),
                       simple_anno_size = unit(0.5, "cm"))


a <- Heatmap(plot_matrix, name = "ratio",
             top_annotation = ha,
             row_split = gene.group,
             column_split = new.group,
             row_gap = unit(2, "mm"),
             column_gap = unit(2, "mm"),
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             column_order = order(factor(info.matrix$pathological_response, levels = c("non-MPR","nPR","pPR","MPR","pCR"))),
             row_names_gp = grid::gpar(fontsize = 10),
             column_names_gp = grid::gpar(fontsize = 5)) 
a
