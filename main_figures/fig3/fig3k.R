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
colnames(Tex_relevant) <- c("sampleID","Tex_number","group","pathological_response")
head(Tex_relevant)

Treg <- read.csv("expanded_CD4Treg_clonotype_number.csv")
Treg <- Treg[,-1]
colnames(Treg) <- c("sampleID","number")
rownames(Treg) <- Treg$sampleID
head(Treg)
dim(Treg)
Treg <- Treg[Tex_relevant$sampleID,]
dim(Treg)

Tex_relevant$Treg_number <- Treg$number
head(Tex_relevant)

ggscatter(Tex_relevant[Tex_relevant$pathological_response %in% c("non-MPR"),], 
          x = "Tex_number",
          y = "Treg_number", color = "group",shape = "group")+
  scale_color_manual(values=c("group1"="#E84C35","group2"="#4FBAD6",
                              "group3"="#00A289","group4"="#3C5487",
                              "group5"="#F29B80")) + 
  xlim(0,60) + ylim(0,60)
