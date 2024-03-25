library(tidyverse)
library(dbplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(readxl)
########################################################################################
#each group major cell type
#################################################################################################
info <- read.csv("/home/zhangwj/data_yi/neoadjuvant/data/scRNA_h5ad_S/all_sub_cell_type.csv")
head(info)
info <- info[info$major_cell_type %in% c("T/NK cells"),]

df <- table(info$sampleID,info$sub_cell_type)
ratio <- as.data.frame(df / rowSums(df))
head(ratio)
colnames(ratio) <- c("sampleID","cell.type","Freq")
head(ratio)

sample.info <- as.data.frame(read_excel("sample.xlsx"))
head(sample.info)
pathological_response <- c()
for(each in sample.info$pathological_response){
  if(each %in% c("MPR","pCR")){
    pathological_response <- c(pathological_response, "MPR")
  }else{
    pathological_response <- c(pathological_response, "non-MPR")
  }
}

sample.info$pathological_response <- pathological_response

response.meta <- sample.info[, c("sampleID","smoking_history","cancer_type","pre_treatment_staging",
                                 "PDL1_TPS","PD1","chemotherapy","targeted_therapy","cycles",
                                 "pathological_response","pathological_response_rate","radiological_response",
                                 "RVT_pre_dominant_histology")]
response.meta <- response.meta %>% distinct(sampleID, .keep_all = TRUE)
head(response.meta)

ratio <- merge(ratio, response.meta, by = "sampleID", all.x = TRUE)
head(ratio)

ratio$pathological_response <- factor(ratio$pathological_response, levels = c("MPR","non-MPR"))

cluster.info <- read.csv("NMF_LUSC_group_5.csv")
cluster.info <- cluster.info[,-1]
head(cluster.info)

ratio <- merge(ratio, cluster.info, by = "sampleID", all.x = TRUE)
head(ratio)
ratio <- ratio[ratio$group %in% c("1","2","3","4","5"),]
head(ratio)
ratio$group <- paste0("group",ratio$group)
ratio$group <- factor(ratio$group, levels = c("group1","group2","group3","group4","group5"))


compaired <- list(c("group1","group2"),c("group2","group3"),
                  c("group2","group4"),c("group2","group5"))
ggboxplot(ratio, x = "group", y = "Freq",
          fill = "group",
          x.text.angle=90,size = 0.3,pt.size = 0,facet.by = "cell.type",
          dot.size = 0) +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = T,test = t.test) + 
  scale_fill_manual(values=c("group1"="#E84C35","group2"="#4FBAD6",
                             "group3"="#00A289","group4"="#3C5487",
                             "group5"="#F29B80")) +
  theme(legend.position="none")

#######################################################################################
########################################################################################
#only myeloid cell
#################################################################################################
info <- read.csv("all_sub_cell_type.csv")
head(info)

info <- info[info$major_cell_type %in% c("Myeloid cell"),]

df <- table(info$sampleID,info$sub_cell_type)
ratio <- as.data.frame(df / rowSums(df))
head(ratio)
colnames(ratio) <- c("sampleID","cell.type","Freq")
head(ratio)

sample.info <- as.data.frame(read_excel("sample.xlsx"))
head(sample.info)
pathological_response <- c()
for(each in sample.info$pathological_response){
  if(each %in% c("MPR","pCR")){
    pathological_response <- c(pathological_response, "MPR")
  }else{
    pathological_response <- c(pathological_response, "non-MPR")
  }
}

sample.info$pathological_response <- pathological_response

response.meta <- sample.info[, c("sampleID","smoking_history","cancer_type","pre_treatment_staging",
                                 "PDL1_TPS","PD1","chemotherapy","targeted_therapy","cycles",
                                 "pathological_response","pathological_response_rate","radiological_response",
                                 "RVT_pre_dominant_histology")]
response.meta <- response.meta %>% distinct(sampleID, .keep_all = TRUE)
head(response.meta)

ratio <- merge(ratio, response.meta, by = "sampleID", all.x = TRUE)
head(ratio)

ratio$pathological_response <- factor(ratio$pathological_response, levels = c("MPR","non-MPR"))

cluster.info <- read.csv("NMF_LUSC_group_5.csv")
cluster.info <- cluster.info[,-1]
head(cluster.info)

ratio <- merge(ratio, cluster.info, by = "sampleID", all.x = TRUE)
head(ratio)
ratio <- ratio[ratio$group %in% c("1","2","3","4","5"),]
head(ratio)
ratio$group <- paste0("group",ratio$group)
ratio$group <- factor(ratio$group, levels = c("group1","group2","group3","group4","group5"))

NN <- ratio[ratio$cell.type == "mDC_LAMP3",]
head(NN)
NN <- NN[,!colnames(NN) %in% c("Freq","cell.type")]
head(NN)
for(each in unique(ratio$cell.type)){
  NN[paste0(each,"_Freq")] <- ratio[ratio$cell.type == each,]$Freq
}
head(NN)

NN$DC_Freq <- NN$cDC1_CLEC9A_Freq + NN$cDC2_CD1C_Freq + NN$mDC_LAMP3_Freq +NN$pDC_LILRA4_Freq
compaired <- list(c("group1", "group2"),c("group2", "group3"),
                  c("group2", "group4"),c("group2", "group5"))


ggboxplot(NN, x = "group", y = "mDC_LAMP3_Freq",
          color = "group",add="jitter",
          x.text.angle=0) + labs(x='group', y= 'mDC_LAMP3 / (all myeloid cells)') + 
  theme(legend.position="none") +
  scale_color_manual(values=c("group1"="#E84C35","group2"="#4FBAD6",
                              "group3"="#00A289","group4"="#3C5487",
                              "group5"="#F29B80"),guide = "none") +  
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
