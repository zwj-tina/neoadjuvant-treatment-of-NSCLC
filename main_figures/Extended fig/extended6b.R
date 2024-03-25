Tex_relevant <- read.csv("CD8Tex_relevant_clonotype_number.csv")
Tex_relevant <- Tex_relevant[,-1]
colnames(Tex_relevant) <- c("sampleID","number")
rownames(Tex_relevant) <- Tex_relevant$sampleID
head(Tex_relevant)
dim(Tex_relevant)

LUSC.group <- read.csv("NMF_LUSC_group_5.csv")
LUSC.group <- LUSC.group[,-1]
rownames(LUSC.group) <- LUSC.group$sampleID
head(LUSC.group)
dim(LUSC.group)

samples <- intersect(Tex_relevant$sampleID, LUSC.group$sampleID)
length(samples)

Tex_relevant <- Tex_relevant[samples,]
LUSC.group <- LUSC.group[samples,]

Tex_relevant$group <- LUSC.group$group
head(Tex_relevant)

library(readxl)
sample.info <- as.data.frame(read_excel("sample.xlsx"))
sample.info <- sample.info[,c("sampleID","pathological_response")]
head(sample.info)
response <- list(
  "pCR" = "MPR",
  "MPR" = "MPR",
  "pPR" = "non-MPR",
  "nPR" = "non-MPR",
  "non-MPR" = "non-MPR"
)
sample.info$pathology <- sapply(as.vector(sample.info$pathological_response), function(x) response[[x]])
head(sample.info)
rownames(sample.info) <- sample.info$sampleID
sample.info <- sample.info[rownames(Tex_relevant),]
head(sample.info)

Tex_relevant$pathological_response <- sample.info$pathology
head(Tex_relevant)
Tex_relevant$group <- paste0("group",Tex_relevant$group)
head(Tex_relevant)



Treg <- read.csv("expanded_CD4Treg_clonotype_number.csv")
Treg <- Treg[,-1]
colnames(Treg) <- c("sampleID","number")
rownames(Treg) <- Treg$sampleID
head(Treg)
dim(Treg)
Treg <- Treg[Tex_relevant$sampleID,]
dim(Treg)
Treg_level <- c()
for(each in Treg$number){
  if(each >= 4){
    Treg_level <- c(Treg_level,"high")
  }else{
    Treg_level <- c(Treg_level,"low")
  }
}
Treg$Treg_level <- Treg_level
head(Treg)

Tex_relevant$Treg_level <- Treg_level

head(Tex_relevant)

Tex_relevant$type <- paste0(Tex_relevant$pathological_response,"_",Tex_relevant$Treg_level)
head(Tex_relevant)

Tex_relevant <- Tex_relevant[Tex_relevant$pathological_response %in% c("non-MPR"),]
head(Tex_relevant)

Tex_relevant$type <- factor(Tex_relevant$type, levels = c("non-MPR_low","non-MPR_high"))
Tex_relevant$number <- as.numeric(Tex_relevant$number)

compaired <- list(c("non-MPR_low","non-MPR_high"))
ggboxplot(Tex_relevant, x = "type", y = "number",
          color = "type",add="jitter",add.params=list(size=0.6),
          x.text.angle=0) + labs(x='group', y= 'number of Tex-relevant clonotypes') + 
  theme(legend.position="none")  + 
  scale_color_manual(values=c("#D9BFAE","#8CB4A3")) + 
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test) 




Treg <- read.csv("expanded_CD4Treg_clonotype_number.csv")
Treg <- Treg[,-1]
colnames(Treg) <- c("sampleID","number")
rownames(Treg) <- Treg$sampleID
head(Treg)
dim(Treg)

LUSC.group <- read.csv("NMF_LUSC_group_5.csv")
LUSC.group <- LUSC.group[,-1]
rownames(LUSC.group) <- LUSC.group$sampleID
head(LUSC.group)
dim(LUSC.group)

samples <- intersect(Treg$sampleID, LUSC.group$sampleID)
length(samples)

Treg <- Treg[samples,]
LUSC.group <- LUSC.group[samples,]

Treg$group <- LUSC.group$group
head(Treg)

library(readxl)
sample.info <- as.data.frame(read_excel("sample.xlsx"))
sample.info <- sample.info[,c("sampleID","pathological_response")]
head(sample.info)
response <- list(
  "pCR" = "MPR",
  "MPR" = "MPR",
  "pPR" = "non-MPR",
  "nPR" = "non-MPR",
  "non-MPR" = "non-MPR"
)
sample.info$pathology <- sapply(as.vector(sample.info$pathological_response), function(x) response[[x]])
head(sample.info)
rownames(sample.info) <- sample.info$sampleID
sample.info <- sample.info[rownames(Treg),]
head(sample.info)

Treg$pathological_response <- sample.info$pathology
head(Treg)
dim(Treg)

Treg_level <- c()
for(each in Treg$number){
  if(each > 3){
    Treg_level <- c(Treg_level,"high")
  }else{
    Treg_level <- c(Treg_level,"low")
  }
}
Treg$Treg_level <- Treg_level
head(Treg)

Treg$type <- paste0(Treg$pathological_response,"_",Treg$Treg_level)
head(Treg)

Treg <- Treg[Treg$pathological_response %in% c("non-MPR"),]
head(Treg)

Treg$type <- factor(Treg$type, levels = c("non-MPR_low","non-MPR_high"))
Treg$number <- as.numeric(Treg$number)

compaired <- list(c("non-MPR_low","non-MPR_high"))
ggboxplot(Treg, x = "type", y = "number",
          color = "type",add="jitter",add.params=list(size=0.6),
          x.text.angle=0) + labs(x='group', y= 'number of Tex-relevant clonotypes') + 
  theme(legend.position="none")  + 
  scale_color_manual(values=c("#D9BFAE","#8CB4A3")) + 
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test) 

