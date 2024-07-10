library(reshape2)
library(tidyverse)
library(dplyr)
library(readxl)
###################################################################
info <- read.csv("all_sub_cell_type.csv")
head(info)
length(unique(info$sampleID))

cluster.info <- read.csv("NMF_all_group_5.csv")
cluster.info <- cluster.info[,-1]
head(cluster.info)
dim(cluster.info)


info <- info[info$sampleID %in% cluster.info$sampleID,]
length(unique(info$sampleID))


df <- table(info$sampleID,info$sub_cell_type)
ratio <- as.data.frame(df / rowSums(df))
head(ratio)
colnames(ratio) <- c("sampleID","cell.type","Freq")
head(ratio)

ratio <- dcast(ratio, sampleID ~ ratio$cell.type, value.var = "Freq")
rownames(ratio) <- ratio$sampleID
head(ratio)

ratio$group1_module <- ratio$`CD8T_NK-like_FGFBP2` + ratio$NK_CD16hi_FGFBP2 + ratio$CD4T_Tm_ANXA1 + ratio$Mφ_FCGR3A
ratio$group2_module <- ratio$Bm_TNFSF9 + ratio$Bm_FCRL4 + ratio$Bm_PDE4D + ratio$Bm_CD74 + ratio$ILC3_KIT + ratio$Bm_TNF + ratio$Bn_TCL1A
ratio$group3_module <- ratio$`CD8T_Tem_GZMK+GZMH+` + ratio$CD8T_Trm_ZNF683 + ratio$`CD8T_Tem_GZMK+NR4A1+` + ratio$CD8T_Tm_IL7R + ratio$CD8T_MAIT_KLRB1
ratio$group4_module <- ratio$CD4T_Treg_FOXP3 + ratio$CD4T_Treg_CCR8 + ratio$CD4T_Tfh_CXCL13 + ratio$`CD4T_Th1-like_CXCL13` + 
  ratio$CD4T_Treg_MKI67 + ratio$CD8T_ISG15 + ratio$CD8T_terminal_Tex_LAYN + ratio$CD8T_Tex_CXCL13
ratio$group5_module <- ratio$Mφ_VCAN + ratio$Mφ_FOLR2 + ratio$cDC2_CD1C + ratio$Mφ_CXCL2 + 
  ratio$Mφ_DNAJB1 + ratio$Mφ_ISG15 + ratio$mDC_LAMP3 + ratio$Mφ_MARCO + ratio$Mφ_CXCL10 + ratio$pDC_LILRA4 + ratio$Mφ_MMP9 + ratio$cDC1_CLEC9A

ratio <- merge(ratio,cluster.info,by = "sampleID")
head(ratio)
ratio$group <- paste0("group",ratio$group)
ratio$group <- factor(ratio$group, levels = c("group1","group2","group3","group4","group5"))

compaired <- list(c("group1","group2"),
                  c("group1","group3"),
                  c("group1","group4"),
                  c("group1","group5"))
ggboxplot(ratio, x = "group", y = "`CD8T_NK-like_FGFBP2`",
          color = "group",add="jitter",add.params=list(size=0.5),
          x.text.angle=0) + labs(y= 'CD8T_NK-like_FGFBP2 / CD45+') + 
  theme(legend.position="none") + 
  scale_color_manual(values=c("group1"="#E84C35","group2"="#4FBAD6",
                              "group3"="#00A289","group4"="#3C5487",
                              "group5"="#F29B80")) +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)




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

a <- table(ratio$pathological_response,ratio$group)
a <- as.data.frame(a / rowSums(a))
colnames(a) <- c("response","group","Freq")
head(a)
ggbarplot(a, x="response", y="Freq", fill = "group",
          x.text.angle=90) + theme(legend.position = "right") + 
  scale_fill_manual(values=c("group1"="#E84C35","group2"="#4FBAD6",
                             "group3"="#00A289","group4"="#3C5487",
                             "group5"="#F29B80")) 

LUSC <- ratio[ratio$cancer_type %in% c("LUSC"),]
a <- table(LUSC$pathological_response,LUSC$group)
a <- as.data.frame(a / rowSums(a))
colnames(a) <- c("response","group","Freq")
head(a)
ggbarplot(a, x="response", y="Freq", fill = "group",
          x.text.angle=90) + theme(legend.position = "right") + 
  scale_fill_manual(values=c("group1"="#E84C35","group2"="#4FBAD6",
                             "group3"="#00A289","group4"="#3C5487",
                             "group5"="#F29B80")) 

LUAD <- ratio[ratio$cancer_type %in% c("LUAD"),]
a <- table(LUAD$pathological_response,LUAD$group)
a <- as.data.frame(a / rowSums(a))
colnames(a) <- c("response","group","Freq")
head(a)
ggbarplot(a, x="response", y="Freq", fill = "group",
          x.text.angle=90) + theme(legend.position = "right") + 
  scale_fill_manual(values=c("group1"="#E84C35","group2"="#4FBAD6",
                             "group3"="#00A289","group4"="#3C5487",
                             "group5"="#F29B80")) 



#PD1
PDL1_TPS_group <- c()
for(each in ratio$PDL1_TPS){
  if (each %in% c("<1%","0")){
    PDL1_TPS_group <- c(PDL1_TPS_group,"<1%") 
  }else if (each %in% c("0.02","0.05","0.08","0.03","0.01","0.3","0.2","0.35","0.4")){
    PDL1_TPS_group <- c(PDL1_TPS_group,"1%-49%") 
  }else if (each %in% c("0.5","0.7","0.6","0.9","0.8","0.85","1","0.55","0.75","0.65")){
    PDL1_TPS_group <- c(PDL1_TPS_group,">=50%") 
  }else{
    PDL1_TPS_group <- c(PDL1_TPS_group,"Not tested") 
  }
}
ratio$PDL1_TPS_group <- PDL1_TPS_group

ratio.part <- ratio[ratio$PDL1_TPS_group != "Not tested",]
ratio.part$PDL1_TPS_group <- factor(as.vector(ratio.part$PDL1_TPS_group), levels = c("<1%","1%-49%",">=50%"))
ratio.part$cluster <- paste0(ratio.part$group,"_",ratio.part$pathological_response)

a <- table(ratio.part$cluster,ratio.part$PDL1_TPS_group)
df <- as.data.frame(a)
head(df)
colnames(df) <- c("cluster","PDL1_TPS","number")
head(df)
df$cluster <- factor(df$cluster,
                     levels = c("group1_MPR","group1_non-MPR",
                                "group2_MPR","group2_non-MPR",
                                "group3_MPR","group3_non-MPR",
                                "group4_MPR","group4_non-MPR",
                                "group5_MPR","group5_non-MPR"))

ggbarplot(df, x="cluster", y="number", fill = "PDL1_TPS",
          x.text.angle=90) + theme(legend.position = "right") + 
  scale_fill_manual(values=c("<1%"="#B3D9D9","1%-49%"="#4F9D9D",
                             ">=50%"="#3D7878"))




group <- c()
response <- c()
for(each in df$cluster){
  group <- c(group, str_split(each, "_")[[1]][1])
  response <- c(response, str_split(each, "_")[[1]][2])
}
df$group <- group
df$response <- response


ggbarplot(df, x="PDL1_TPS", y="number", fill = "group",
          x.text.angle=90,facet.by = "response") + theme(legend.position = "right") + 
  scale_fill_manual(values=c("group1"="#E84C35","group2"="#4FBAD6",
                             "group3"="#00A289","group4"="#3C5487",
                             "group5"="#F29B80"))
ggsave("/home/zhangwj/data_yi/neoadjuvant/revision2/figure/PDL1_2.pdf",width = 5, height = 4)

#alluvial
library(readxl)
sample.info <- as.data.frame(read_excel("sample.xlsx"))
head(sample.info)

cluster.info <- read.csv("NMF_all_group_5.csv")
cluster.info <- cluster.info[,-1]
head(cluster.info)

sample.info <- merge(sample.info, cluster.info, by = "sampleID", all.x = TRUE)
head(sample.info)
sample.info <- sample.info[sample.info$group %in% c("1","2","3","4","5"),]
head(sample.info)

sample.info$group <- paste0("group",sample.info$group)

pathological_response <- c()
for(each in sample.info$pathological_response){
  if(each %in% c("MPR","pCR")){
    pathological_response <- c(pathological_response, "MPR")
  }else{
    pathological_response <- c(pathological_response, "non-MPR")
  }
}

sample.info$pathological_response <- pathological_response
sample.info$cluster <- paste0(sample.info$group,"_",sample.info$pathological_response)

#PD1
PDL1_TPS_group <- c()
for(each in sample.info$PDL1_TPS){
  if (each %in% c("<1%","0")){
    PDL1_TPS_group <- c(PDL1_TPS_group,"<1%") 
  }else if (each %in% c("0.01","0.02","0.05","0.08","0.03","0.3","0.2","0.35","0.4")){
    PDL1_TPS_group <- c(PDL1_TPS_group,"1%-49%") 
  }else if (each %in% c("0.5","0.7","0.6","0.9","0.8","0.85","1","0.55","0.75","0.65")){
    PDL1_TPS_group <- c(PDL1_TPS_group,">=50%") 
  }else{
    PDL1_TPS_group <- c(PDL1_TPS_group,"Not tested") 
  }
}
sample.info$PDL1_TPS_group <- PDL1_TPS_group

sample.info <- sample.info[sample.info$PDL1_TPS_group != "Not tested",]
sample.info$PDL1_TPS_group <- factor(as.vector(sample.info$PDL1_TPS_group), levels = c("<1%","1%-49%",">=50%"))
a <- table(sample.info$cluster,sample.info$PDL1_TPS_group)
df <- as.data.frame(a)
head(df)
colnames(df) <- c("cluster","PDL1_TPS","number")
head(df)
df$cluster <- factor(df$cluster,
                     levels = c("group1_MPR","group1_non-MPR",
                                "group2_MPR","group2_non-MPR",
                                "group3_MPR","group3_non-MPR",
                                "group4_MPR","group4_non-MPR",
                                "group5_MPR","group5_non-MPR"))


df$group <- sapply(as.vector(df$cluster), function(x) strsplit(x,"_")[[1]][1])
head(df)
df$pathological_response <- sapply(as.vector(df$cluster), function(x) strsplit(x,"_")[[1]][2])
head(df)
df <- df[!df$number %in% c(0),]
ggplot(data = df,
       aes(axis1 = PDL1_TPS,   # First variable on the X-axis
           axis2 = group,   # Third variable on the X-axis
           y = number)) +
  geom_alluvium(aes(fill = pathological_response,order = pathological_response)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void() +
  scale_fill_manual(values=c("MPR" = "#2868A6", "non-MPR" = "#B1161C"))

library(ggalluvial)
cluster.info <- read.csv("NMF_all_group_5.csv")
cluster.info <- cluster.info[,-1]
head(cluster.info)

sample.info <- as.data.frame(read_excel("sample.xlsx"))
head(sample.info)
sample.info <- sample.info[sample.info$sampleID %in% cluster.info$sampleID,]
sample.info <- sample.info[,c("sampleID","pathological_response","radiological_response")]

pathological_response_level <- c()
for(each in sample.info$pathological_response){
  if(each %in% c("MPR","pCR")){
    pathological_response_level <- c(pathological_response_level, "MPR")
  }else{
    pathological_response_level <- c(pathological_response_level, "non-MPR")
  }
}
sample.info$pathological_response_level <- pathological_response_level
rownames(sample.info) <- sample.info$sampleID
head(sample.info)
sample.info <- sample.info[cluster.info$sampleID,]
sample.info$group <- cluster.info$group
sample.info$group <- paste0("group",sample.info$group)
sample.info$sub.group <- paste0(sample.info$group,"_",sample.info$pathological_response_level)
sample.info$sub.group <- factor(sample.info$sub.group,
                                levels = c("group1_MPR","group1_non-MPR",
                                           "group2_MPR","group2_non-MPR",
                                           "group3_MPR","group3_non-MPR",
                                           "group4_MPR","group4_non-MPR",
                                           "group5_MPR","group5_non-MPR"))
head(sample.info)

sample.info <- sample.info[sample.info$radiological_response %in% c("SD","PR","CR","PD"),]
head(sample.info)
mm <- as.data.frame(table(sample.info$sub.group,sample.info$radiological_response))
colnames(mm) <- c("sub.group","radiological_response","number")
mm$radiological_response <- factor(mm$radiological_response,levels = c("CR","PR","SD","PD"))
head(mm)

ggplot(data = mm,
       aes(axis1 = sub.group,   # First variable on the X-axis
           axis2 = radiological_response,   # Third variable on the X-axis
           y = number)) +
  geom_alluvium(aes(fill = radiological_response,order = radiological_response)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void() +
  scale_fill_manual(values=c("CR"="#E6A4B4","PR"="#FFD9C0","SD"="#8CC0DE","PD"="#0B60B0"))

