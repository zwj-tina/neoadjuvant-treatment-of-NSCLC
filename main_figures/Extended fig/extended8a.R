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

compaired <- list(c("group3_pPR","group4_pPR"))
ggboxplot(LUAD.group[LUAD.group$group %in% c("group3","group4"),], x = "type", y = "pathological_response_rate",
          color = "pathological_response",add="jitter",add.params=list(size=0.6),
          x.text.angle=0) + labs(x='group', y= 'pathological response rate') + 
  theme(legend.position="none")  + 
  scale_color_manual(values=c("nPR" = "#EDB752","pPR" = "#FFC1B3","MPR" = "#AFD483","pCR" = "#469832")) + 
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)

LUAD.group$pathological_response <- factor(LUAD.group$pathological_response,levels = c("pCR","MPR","pPR","nPR"))
ggboxplot(LUAD.group[LUAD.group$group %in% c("group3","group4"),], x = "group", y = "pathological_response_rate",
          color = "pathological_response",add="jitter",add.params=list(size=0.6),
          x.text.angle=0) + labs(x='group', y= 'pathological response rate') + 
  scale_color_manual(values=c("nPR" = "#EDB752","pPR" = "#FFC1B3","MPR" = "#AFD483","pCR" = "#469832")) + 
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
