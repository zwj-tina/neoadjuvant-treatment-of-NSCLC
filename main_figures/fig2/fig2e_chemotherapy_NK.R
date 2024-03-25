library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(dplyr)

info <- read.csv("all_sub_cell_type.csv")
head(info)
df <- table(info$sampleID,info$sub_cell_type)
ratio <- as.data.frame(df / rowSums(df))
head(ratio)
colnames(ratio) <- c("sampleID","cell.type","Freq")
head(ratio)

library(readxl)
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

ratio <- dcast(ratio, sampleID ~ ratio$cell.type, value.var = "Freq")

head(ratio)
merge.version <- merge(ratio, response.meta, by = "sampleID", all.x = TRUE)
head(merge.version)
#patients with only chemotherapy
merge.version <- merge.version[merge.version$PD1 %in% c("No"),]
head(merge.version)


ggboxplot(merge.version, x = "pathological_response", y = "`CD8T_NK-like_FGFBP2`",
          color = "pathological_response",add = "jitter",
          x.text.angle=0,size = 0.5,pt.size = 1) +
  stat_compare_means(aes(group = pathological_response)) + 
  scale_color_manual(values = c("MPR" = "#2868A6", "non-MPR" = "#B1161C")) + 
  theme(legend.position="none")
