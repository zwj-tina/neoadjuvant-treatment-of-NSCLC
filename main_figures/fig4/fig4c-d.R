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
head(sample.info)
unique(sample.info$pathological_response)
sample.info <- sample.info[sample.info$pathological_response %in% c("pPR","pCR","MPR","nPR"),]
sample.info$pathological_response <- factor(sample.info$pathological_response,
                                            levels = c("pPR","pCR","MPR","nPR"))
response <- list(
  "pCR" = "MPR",
  "MPR" = "MPR",
  "pPR" = "pPR",
  "nPR" = "nPR"
)
sample.info$pathology <- unlist(sapply(as.vector(sample.info$pathological_response), function(x) response[[x]]))
head(sample.info)
rownames(sample.info) <- sample.info$sampleID
sample.info <- sample.info[rownames(Tex_relevant),]
head(sample.info)

Tex_relevant$pathological_response <- sample.info$pathology
unique(Tex_relevant$pathological_response)
Tex_relevant$group <- paste0("group",Tex_relevant$group)
head(Tex_relevant)

Tex_relevant$type <- paste0(Tex_relevant$group,"_",Tex_relevant$pathological_response)
unique(Tex_relevant$type)

Tex_relevant <- Tex_relevant[Tex_relevant$type %in% c("group3_pPR","group2_MPR",
                                                      "group1_MPR","group4_nPR",
                                                      "group2_nPR","group3_MPR",
                                                      "group5_nPR","group4_pPR",
                                                      "group3_nPR","group1_nPR",
                                                      "group5_MPR"),]

head(Tex_relevant)
Tex_relevant$type <- factor(as.vector(Tex_relevant$type), levels = c("group1_MPR","group1_nPR",
                                                                     "group2_MPR","group2_nPR",
                                                                     "group3_MPR","group3_pPR","group3_nPR",
                                                                     "group4_pPR","group4_nPR",
                                                                     "group5_MPR","group5_nPR"))

compaired <- list(c("group1_MPR","group1_nPR"),
                  c("group2_MPR","group2_nPR"),
                  c("group3_MPR","group3_pPR"),
                  c("group3_pPR","group3_nPR"),
                  c("group4_pPR","group4_nPR"),
                  c("group5_MPR","group5_nPR"))

ggboxplot(Tex_relevant, x = "type", y = "number",
          color = "type",add="jitter",
          x.text.angle=90) + labs(x='group', y= 'clonotype number of Tex-relevent cells') + 
  theme(legend.position="none") + 
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)

Tex_relevant$pathological_response <- factor(Tex_relevant$pathological_response, levels = c("MPR","pPR","nPR"))
Tex_relevant$group <- factor(Tex_relevant$group, levels = c("group1","group2",
                                                            "group3","group4",
                                                            "group5"))
ggboxplot(Tex_relevant, x = "group", y = "number",
          color = "pathological_response",add = "jitter",
          x.text.angle=0,size = 0.5,add.params=list(size=0.6)) + 
  scale_color_manual(values=c("nPR" = "#EDB752","pPR" = "#FFC1B3","MPR" = "#AFD483")) + 
  labs(y= 'number of Tex-relevant clonotypes')


#########################################################################################################
Treg <- read.csv("expanded_CD4Treg_clonotype_number.csv")
Treg <- Treg[,-1]
colnames(Treg) <- c("sampleID","number")
rownames(Treg) <- Treg$sampleID
head(Treg)
dim(Treg)

LUAD.group <- read.csv("NMF_LUAD_group_5.csv")
LUAD.group <- LUAD.group[,-1]
rownames(LUAD.group) <- LUAD.group$sampleID
head(LUAD.group)
dim(LUAD.group)

samples <- intersect(Treg$sampleID, LUAD.group$sampleID)
length(samples)

Treg <- Treg[samples,]
LUAD.group <- LUAD.group[samples,]

Treg$group <- LUAD.group$group
head(Treg)

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
sample.info$pathology <- unlist(sapply(as.vector(sample.info$pathological_response), function(x) response[[x]]))
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
Treg <- Treg[Treg$type %in% c("group3_pPR","group2_MPR",
                                                      "group1_MPR","group4_nPR",
                                                      "group2_nPR","group3_MPR",
                                                      "group5_nPR","group4_pPR",
                                                      "group3_nPR","group1_nPR",
                                                      "group5_MPR"),]

head(Treg)
Treg$type <- factor(as.vector(Treg$type), levels = c("group1_MPR","group1_nPR",
                                                                     "group2_MPR","group2_nPR",
                                                                     "group3_MPR","group3_pPR","group3_nPR",
                                                                     "group4_pPR","group4_nPR",
                                                                     "group5_MPR","group5_nPR"))


compaired <- list(c("group3_MPR", "group3_nPR"),
                  c("group3_pPR", "group3_nPR"),
                  c("group3_nPR", "group4_nPR"),
                  c("group4_pPR", "group4_nPR"))

ggboxplot(Treg, x = "type", y = "number",
          color = "type",add="jitter",
          x.text.angle=90) + labs(x='group', y= 'clonotype number of activated Treg clonotypes') + 
  theme(legend.position="none") + 
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)



Treg$pathological_response <- factor(Treg$pathological_response, levels = c("MPR","pPR","nPR"))
Treg$group <- factor(Treg$group, levels = c("group1","group2",
                                            "group3","group4",
                                            "group5"))
ggboxplot(Treg, x = "group", y = "number",
          color = "pathological_response",add = "jitter",
          x.text.angle=0,size = 0.5,add.params=list(size=0.6)) + 
  scale_color_manual(values=c("nPR" = "#EDB752","pPR" = "#FFC1B3","MPR" = "#AFD483")) + 
  labs(y= 'clonotype number of activated Treg clonotypes')

