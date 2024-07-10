Treg_clone <- read.csv("expanded_CD4Treg_clonotype_number_over2.csv")
dim(Treg_clone)
Treg_clone <- Treg_clone[,-1]
colnames(Treg_clone) <- c("sampleID","number")
rownames(Treg_clone) <- Treg_clone$sampleID
head(Treg_clone)


cluster.info <- read.csv("NMF_all_group_5.csv")
dim(cluster.info)
cluster.info <- cluster.info[,-1]
rownames(cluster.info) <- cluster.info$sampleID
head(cluster.info)

#include samples with TCR data and in our clustering data
samples <- intersect(Treg_clone$sampleID, cluster.info$sampleID)
length(samples)

cluster.info <- cluster.info[samples,]
Treg_clone <- Treg_clone[samples,]

Treg_clone$group <- paste0("group",cluster.info$group)
Treg_clone$group <- factor(Treg_clone$group,levels = c("group1","group2","group3","group4","group5"))
head(Treg_clone)



library(readxl)
sample.info <- as.data.frame(read_excel("sample.xlsx"))
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
sample.info <- sample.info[sample.info$sampleID %in% c(Treg_clone$sampleID),]
dim(sample.info)
head(sample.info)

Treg_clone <- merge(Treg_clone,sample.info,by = "sampleID")
head(Treg_clone)

new_group <- c()
for(each in 1:length(Treg_clone$sampleID)){
  print(Treg_clone$group[[each]])
  if(Treg_clone$group[each] %in% c("group3")){
    new_group <- c(new_group, paste0(Treg_clone$group[[each]],"_",Treg_clone$cancer_type[[each]]))
  }else{
    new_group <- c(new_group,as.vector(Treg_clone$group)[[each]])
  }
}
Treg_clone$new_group <- new_group
head(Treg_clone)

head(Treg_clone)

Treg_clone$new_group <- factor(Treg_clone$new_group,levels = c("group1",
                                                                   "group2",
                                                                   "group3_LUSC",
                                                                   "group3_LUAD",
                                                                   "group4",
                                                                   "group5"))
ggboxplot(Treg_clone, x = "new_group", y = "number",
          color = "pathology",add="jitter",add.params=list(size=0.5),
          x.text.angle=90) + labs(x='group', y= 'clonotype number of Tex-relevent cells') +
  scale_color_manual(values=c("MPR" = "#2868A6", "non-MPR" = "#B1161C")) + 
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)

Treg_clone$cluster <- paste0(Treg_clone$new_group,"_",Treg_clone$pathology)
Treg_clone$cluster <- factor(Treg_clone$cluster,levels = c("group1_MPR","group1_non-MPR",
                                                               "group2_MPR","group2_non-MPR",
                                                               "group3_LUSC_MPR","group3_LUSC_non-MPR",
                                                               "group3_LUAD_MPR","group3_LUAD_non-MPR",
                                                               "group4_MPR","group4_non-MPR",
                                                               "group5_MPR","group5_non-MPR"))

compaired <- list(c("group2_MPR","group2_non-MPR"),
                  c("group3_LUSC_MPR","group3_LUSC_non-MPR"),
                  c("group3_LUAD_MPR","group3_LUAD_non-MPR"),
                  c("group5_MPR","group5_non-MPR"),
                  c("group3_LUSC_non-MPR","group3_LUAD_non-MPR"),
                  c("group4_non-MPR","group3_LUSC_non-MPR"),
                  c("group4_non-MPR","group3_LUAD_non-MPR"),
                  c("group4_non-MPR","group3_LUSC_MPR"),
                  c("group4_non-MPR","group3_LUAD_MPR"))

ggboxplot(Treg_clone, x = "cluster", y = "number",
          color = "cluster",add="jitter",add.params=list(size=0.5),
          x.text.angle=90) + labs(x='group', y= 'clonotype number of Tex-relevent cells') +
  theme(legend.position="none") +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)
