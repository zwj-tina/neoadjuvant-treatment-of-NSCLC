LUAD.group <- read.csv("NMF_LUAD_group_5.csv")
LUAD.group <- LUAD.group[,-1]
rownames(LUAD.group) <- LUAD.group$sampleID
head(LUAD.group)
dim(LUAD.group)

library(readxl)
sample.info <- as.data.frame(read_excel("sample.xlsx"))
sample.info <- sample.info[,c("sampleID","pathological_response","pathological_response_rate")]
head(sample.info)
rownames(sample.info) <- sample.info$sampleID
sample.info <- sample.info[LUAD.group$sampleID,]
head(sample.info)

LUAD.group$pathological_response <- sample.info$pathological_response
LUAD.group$pathological_response_rate <- sample.info$pathological_response_rate
head(LUAD.group)
LUAD.group <- LUAD.group[LUAD.group$pathological_response %in% c("pCR","MPR","pPR","nPR"),]
LUAD.group$pathological_response_rate <- as.numeric(LUAD.group$pathological_response_rate)
LUAD.group$group <- paste0("group",LUAD.group$group)
LUAD.group$group <- factor(LUAD.group$group, levels = c("group1","group2","group3","group4","group5"))
LUAD.group$type <- paste0(LUAD.group$group,"_",LUAD.group$pathological_response)

Treg <- read.csv("expanded_CD4Treg_clonotype_number.csv")
Treg <- Treg[,-1]
colnames(Treg) <- c("sampleID","number")
rownames(Treg) <- Treg$sampleID

samples <- intersect(LUAD.group$sampleID,Treg$sampleID)
length(samples)
#56

LUAD.group <- LUAD.group[samples,]
Treg <- Treg[samples,]

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

LUAD.group$Treg_clonotype_number <- Treg$number
LUAD.group$Treg_level <- Treg$Treg_level
head(LUAD.group)


Tex_relevant <- read.csv("CD8Tex_relevant_clonotype_number.csv")
Tex_relevant <- Tex_relevant[,-1]
colnames(Tex_relevant) <- c("sampleID","number")
rownames(Tex_relevant) <- Tex_relevant$sampleID
head(Tex_relevant)
Tex_relevant <- Tex_relevant[LUAD.group$sampleID,]

LUAD.group$Tex_clonotype_number <- Tex_relevant$number
head(LUAD.group)
LUAD.group$tolerance <- paste0(LUAD.group$pathological_response,LUAD.group$Treg_level)
LUAD.group <- LUAD.group[LUAD.group$pathological_response %in% c("pPR","nPR"),]
unique(LUAD.group$tolerance)
LUAD.group$tolerance <- factor(LUAD.group$tolerance, levels = c("pPRhigh","pPRlow",
                                                                "nPRhigh","nPRlow"))

#pPRhigh  pPRlow nPRhigh  nPRlow 
#6       8      11      10 
compaired <- list(c("pPRhigh","pPRlow"),
                  c("nPRhigh","nPRlow"))
ggboxplot(LUAD.group, x = "tolerance", y = "Tex_clonotype_number",
          color = "tolerance",add="jitter",
          x.text.angle=0) + labs(x='group', y= 'number of Tex-relevant clonotypes') + 
  theme(legend.position="none")  + 
  scale_color_manual(values=c("#D9BFAE","#8CB4A3","#F2C57C","#5E9EA0"))+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test) 


compaired <- list(c("pPRhigh","pPRlow"),
                  c("nPRhigh","nPRlow"))
ggboxplot(LUAD.group, x = "tolerance", y = "Treg_clonotype_number",
          color = "tolerance",add="jitter",
          x.text.angle=0) + labs(x='group', y= 'number of activated Treg clonotypes') + 
  theme(legend.position="none")  + 
  scale_color_manual(values=c("#D9BFAE","#8CB4A3","#F2C57C","#5E9EA0"))+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test) 
