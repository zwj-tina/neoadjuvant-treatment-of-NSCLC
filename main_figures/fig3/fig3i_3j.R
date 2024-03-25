library(ggpubr)

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

Tex_relevant$type <- paste0(Tex_relevant$group,"_",Tex_relevant$pathological_response)
head(Tex_relevant)
Tex_relevant$type <- factor(as.vector(Tex_relevant$type), levels = c("group1_MPR","group1_non-MPR",
                                                               "group2_MPR","group2_non-MPR",
                                                               "group3_MPR","group3_non-MPR",
                                                               "group4_MPR","group4_non-MPR",
                                                               "group5_MPR","group5_non-MPR"))

compaired <- list(c("group1_non-MPR", "group1_MPR"),
                  c("group2_non-MPR", "group2_MPR"),
                  c("group3_non-MPR", "group3_MPR"),
                  c("group4_non-MPR", "group4_MPR"),
                  c("group5_non-MPR", "group5_MPR"),
                  c("group3_non-MPR", "group4_non-MPR"),
                  c("group2_non-MPR", "group4_non-MPR"),
                  c("group5_non-MPR", "group4_non-MPR"),
                  c("group1_non-MPR", "group4_non-MPR"))

ggboxplot(Tex_relevant, x = "type", y = "number",
          color = "type",add="jitter",
          x.text.angle=90) + labs(x='group', y= 'clonotype number of Tex-relevent cells') + 
  theme(legend.position="none") + 
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)


Tex_relevant$pathological_response <- factor(Tex_relevant$pathological_response, levels = c("MPR","non-MPR"))
Tex_relevant$group <- factor(Tex_relevant$group, levels = c("group1","group2",
                                                      "group3","group4",
                                                      "group5"))
ggboxplot(Tex_relevant, x = "group", y = "number",
          color = "pathological_response",add = "jitter",
          x.text.angle=0,size = 0.5,add.params=list(size=0.6)) + 
  scale_color_manual(values=c("MPR" = "#2868A6", "non-MPR" = "#B1161C")) + 
  labs(y= 'number of Tex-relevant clonotypes')




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
Treg$group <- paste0("group",Treg$group)
head(Treg)

Treg$type <- paste0(Treg$group,"_",Treg$pathological_response)
head(Treg)
Treg$type <- factor(as.vector(Treg$type), levels = c("group1_MPR","group1_non-MPR",
                                                                     "group2_MPR","group2_non-MPR",
                                                                     "group3_MPR","group3_non-MPR",
                                                                     "group4_MPR","group4_non-MPR",
                                                                     "group5_MPR","group5_non-MPR"))

compaired <- list(c("group3_non-MPR", "group4_non-MPR"),
                  c("group2_non-MPR", "group4_non-MPR"),
                  c("group5_non-MPR", "group4_non-MPR"),
                  c("group1_non-MPR", "group4_non-MPR"))

ggboxplot(Treg, x = "type", y = "number",
          color = "type",add="jitter",
          x.text.angle=90) + labs(x='group', y= 'clonotype number of activated Treg clonotypes') + 
  theme(legend.position="none") + 
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)


Treg$pathological_response <- factor(Treg$pathological_response, levels = c("MPR","non-MPR"))
Treg$group <- factor(Treg$group, levels = c("group1","group2",
                                                            "group3","group4",
                                                            "group5"))
ggboxplot(Treg, x = "group", y = "number",
          color = "pathological_response",add = "jitter",
          x.text.angle=0,size = 0.5,add.params=list(size=0.6)) + 
  scale_color_manual(values=c("MPR" = "#2868A6", "non-MPR" = "#B1161C")) + 
  labs(y= 'clonotype number of activated Treg clonotypes')
