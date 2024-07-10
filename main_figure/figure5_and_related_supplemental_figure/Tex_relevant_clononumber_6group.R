Tex_relevant <- read.csv("CD8Tex_relevant_clonotype_number_over2.csv")
dim(Tex_relevant)
Tex_relevant <- Tex_relevant[,-1]
colnames(Tex_relevant) <- c("sampleID","number")
rownames(Tex_relevant) <- Tex_relevant$sampleID
head(Tex_relevant)


cluster.info <- read.csv("NMF_all_group_5.csv")
dim(cluster.info)
cluster.info <- cluster.info[,-1]
rownames(cluster.info) <- cluster.info$sampleID
head(cluster.info)

#include samples with TCR data and in our clustering data
samples <- intersect(Tex_relevant$sampleID, cluster.info$sampleID)
length(samples)

cluster.info <- cluster.info[samples,]
Tex_relevant <- Tex_relevant[samples,]

Tex_relevant$group <- paste0("group",cluster.info$group)
Tex_relevant$group <- factor(Tex_relevant$group,levels = c("group1","group2","group3","group4","group5"))
head(Tex_relevant)



library(readxl)
sample.info <- as.data.frame(read_excel("neoadjuvant/data/other/sample.xlsx"))
sample.info <- sample.info[,c("sampleID","pathological_response","cancer_type")]
pathology <- c()
for(each in sample.info$pathological_response){
  if(each %in% c("MPR","pCR")){
    pathology <- c(pathology, "MPR")
  }else{
    pathology <- c(pathology, "non-MPR")
  }
}
sample.info$pathology <- pathology
rownames(sample.info) <- sample.info$sampleID
sample.info <- sample.info[sample.info$sampleID %in% c(Tex_relevant$sampleID),]
dim(sample.info)
head(sample.info)

Tex_relevant <- merge(Tex_relevant,sample.info,by = "sampleID")
head(Tex_relevant)

new_group <- c()
for(each in 1:length(Tex_relevant$sampleID)){
  print(Tex_relevant$group[[each]])
  if(Tex_relevant$group[each] %in% c("group3")){
    new_group <- c(new_group, paste0(Tex_relevant$group[[each]],"_",Tex_relevant$cancer_type[[each]]))
  }else{
    new_group <- c(new_group,as.vector(Tex_relevant$group)[[each]])
  }
}
Tex_relevant$new_group <- new_group
head(Tex_relevant)
Tex_relevant$new_group <- factor(Tex_relevant$new_group,levels = c("group1",
                                                               "group2",
                                                               "group3_LUSC",
                                                               "group3_LUAD",
                                                               "group4",
                                                               "group5"))
ggboxplot(Tex_relevant, x = "new_group", y = "number",
          color = "pathology",add="jitter",add.params=list(size=0.5),
          x.text.angle=90) + labs(x='group', y= 'clonotype number of Tex-relevent cells') +
  scale_color_manual(values=c("MPR" = "#2868A6", "non-MPR" = "#B1161C")) + 
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)

Tex_relevant$cluster <- paste0(Tex_relevant$new_group,"_",Tex_relevant$pathology)
Tex_relevant$cluster <- factor(Tex_relevant$cluster,levels = c("group1_MPR","group1_non-MPR",
                                                       "group2_MPR","group2_non-MPR",
                                                       "group3_LUSC_MPR","group3_LUSC_non-MPR",
                                                       "group3_LUAD_MPR","group3_LUAD_non-MPR",
                                                       "group4_MPR","group4_non-MPR",
                                                       "group5_MPR","group5_non-MPR"))

compaired <- list(c("group2_MPR","group2_non-MPR"),
                  c("group3_LUSC_MPR","group3_LUSC_non-MPR"),
                  c("group3_LUAD_MPR","group3_LUAD_non-MPR"),
                  c("group5_MPR","group5_non-MPR"),
                  c("group4_non-MPR","group3_LUSC_non-MPR"),
                  c("group4_non-MPR","group3_LUAD_non-MPR"),
                  c("group3_LUSC_non-MPR","group3_LUAD_non-MPR"),
                  c("group4_non-MPR","group3_LUSC_MPR"),
                  c("group4_non-MPR","group3_LUAD_MPR"))

ggboxplot(Tex_relevant, x = "cluster", y = "number",
          color = "cluster",add="jitter",add.params=list(size=0.5),
          x.text.angle=90) + labs(x='group', y= 'clonotype number of Tex-relevent cells') +
  theme(legend.position="none") +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
