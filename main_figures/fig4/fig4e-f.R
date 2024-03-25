library(ggpubr)
Tex_relevant <- read.csv("CD8Tex_relevant_clonotype_number.csv")
Tex_relevant <- Tex_relevant[,-1]
colnames(Tex_relevant) <- c("sampleID","number")
rownames(Tex_relevant) <- Tex_relevant$sampleID
head(Tex_relevant)
dim(Tex_relevant)

LUAD.group <- read.csv("NMF_LUAD_group_5.csv")
LUAD.group <- LUAD.group[,-1]
rownames(LUAD.group) <- LUAD.group$sampleID
head(LUAD.group)
dim(LUAD.group)

samples <- intersect(Tex_relevant$sampleID, LUAD.group$sampleID)
length(samples)

Tex_relevant <- Tex_relevant[samples,]
LUAD.group <- LUAD.group[samples,]

Tex_relevant$group <- LUAD.group$group
head(Tex_relevant)

library(readxl)
sample.info <- as.data.frame(read_excel("sample.xlsx"))
sample.info <- sample.info[,c("sampleID","pathological_response")]
sample.info <- sample.info[sample.info$pathological_response %in% c("pPR","pCR","MPR","nPR"),]
sample.info$pathological_response <- factor(sample.info$pathological_response,
                                            levels = c("pPR","pCR","MPR","nPR"))
response <- list(
  "pCR" = "MPR",
  "MPR" = "MPR",
  "pPR" = "pPR",
  "nPR" = "nPR"
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

ggscatter(Tex_relevant[Tex_relevant$pathological_response %in% c("nPR"),], 
          x = "Tex_number",
          y = "Treg_number", color = "group",shape = "group")+
  scale_color_manual(values=c("group1"="#E84C35","group2"="#4FBAD6",
                              "group3"="#00A289","group4"="#3C5487",
                              "group5"="#F29B80")) + 
  xlim(0,40) + ylim(0,40) + geom_hline(aes(yintercept=4), colour="#990000", linetype="dashed")
