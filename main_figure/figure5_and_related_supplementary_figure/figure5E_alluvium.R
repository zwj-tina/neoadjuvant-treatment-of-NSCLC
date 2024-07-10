library(ggalluvial)
cluster.info <- read.csv("/home/zhangwj/data_yi/neoadjuvant/revision2/data/NMF_all_group_5.csv")
cluster.info <- cluster.info[,-1]
head(cluster.info)

sample.info <- as.data.frame(read_excel("/home/zhangwj/data_yi/neoadjuvant/data/other/sample.xlsx"))
head(sample.info)
sample.info <- sample.info[sample.info$sampleID %in% cluster.info$sampleID,]
sample.info <- sample.info[,c("sampleID","pathological_response","cancer_type")]

pathological_response_level <- c()
for(each in sample.info$pathological_response){
  if(each %in% c("MPR","pCR")){
    pathological_response_level <- c(pathological_response_level, "MPR")
  }else{
    pathological_response_level <- c(pathological_response_level, "non-MPR")
  }
}
sample.info$pathological_response_level <- pathological_response_level
rownames(sample.info) <- sample.info$sampleID
head(sample.info)
sample.info <- sample.info[cluster.info$sampleID,]
sample.info$group <- cluster.info$group
sample.info$group <- paste0("group",sample.info$group)
sample.info$sub.group <- paste0(sample.info$group,"_",sample.info$pathological_response_level)
sample.info$sub.group <- factor(sample.info$sub.group,
                                levels = c("group1_MPR","group1_non-MPR",
                                           "group2_MPR","group2_non-MPR",
                                           "group3_MPR","group3_non-MPR",
                                           "group4_MPR","group4_non-MPR",
                                           "group5_MPR","group5_non-MPR"))
head(sample.info)


Treg_clone <- read.csv("expanded_CD4Treg_clonotype_number_over2.csv")
dim(Treg_clone)
Treg_clone <- Treg_clone[,-1]
colnames(Treg_clone) <- c("sampleID","number")
rownames(Treg_clone) <- Treg_clone$sampleID
head(Treg_clone)

Treg_clone <- Treg_clone[Treg_clone$sampleID %in% cluster.info$sampleID,]
dim(Treg_clone)
Treg_level <- c()
for(each in Treg_clone$number){
  print(each)
  if(each >= 9){
    Treg_level <- c(Treg_level,"high")
  }else{
    Treg_level <- c(Treg_level,"low")
  }
}
Treg_clone$Treg_level <- Treg_level
colnames(Treg_clone) <- c("sampleID","Treg_number","Treg_level")
head(Treg_clone)
dim(Treg_clone)

head(sample.info)

sample.info <- merge(sample.info, Treg_clone,by = "sampleID")
head(sample.info)

cluster <- c()
for(i in 1:length(sample.info$sampleID)){
  if(sample.info$pathological_response_level[[i]] %in% c("MPR")){
    cluster <- c(cluster, "MPR")
  }else{
    cluster <- c(cluster, paste0(sample.info$pathological_response_level[[i]],"_",sample.info$Treg_level[[i]]))
  }
}
unique(cluster)

sample.info$cluster <- cluster
sample.info$cluster <- factor(sample.info$cluster, levels = c("MPR","non-MPR_high","non-MPR_low"))
head(sample.info)

LUSC <- sample.info[sample.info$cancer_type %in% c("LUSC"),]

mm <- as.data.frame(table(sample.info$sub.group,sample.info$cluster))
colnames(mm) <- c("sub.group","tolerance","number")
mm$tolerance <- factor(mm$tolerance,levels = c("MPR","non-MPR_high",
                                               "non-MPR_low"))
head(mm)

ggplot(data = mm,
       aes(axis1 = sub.group,   # First variable on the X-axis
           axis2 = tolerance,   # Third variable on the X-axis
           y = number)) +
  geom_alluvium(aes(fill = tolerance,order = tolerance)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void() +
  scale_fill_manual(values=c("MPR"="#D9BFAE",
                             "non-MPR_high"="#8CB4A3",
                             "non-MPR_low"="#7998AD"))
