###################################################################################################
#############################################################################################################
#number of each clone type Tex-relevant
#############################################################################################################
TCR.data <- read.csv("TCR_with_new_group.csv")
head(TCR.data)

cluster.info <- read.csv("NMF_LUSC_group_5.csv")
cluster.info <- cluster.info[,-1]
head(cluster.info)
dim(cluster.info)

#include samples with TCR data and in our clustering data
samples <- intersect(TCR.data$sampleID, cluster.info$sampleID)
length(samples)

cluster.info <- cluster.info[cluster.info$sampleID %in% samples,]
TCR.data <- TCR.data[TCR.data$sampleID %in% samples,]
length(unique(TCR.data$sampleID))

#only include cells with shared clones with expanded terminal Tex
TCR.data <- TCR.data[TCR.data$T_new_name %in% c("CD8Texp","expanded_terminal_Tex"),]
dim(TCR.data)

library(readxl)
sample.info <- as.data.frame(read_excel("sample.xlsx"))
response <- list(
  "pCR" = "MPR",
  "MPR" = "MPR",
  "pPR" = "non-MPR",
  "nPR" = "non-MPR",
  "non-MPR" = "non-MPR"
)

sample.info$pathology <- sapply(as.vector(sample.info$pathological_response), function(x) response[[x]])
sample.info

count = 0
for(each in unique(TCR.data$sampleID)){
  P <- TCR.data[TCR.data$sampleID == each,c("clonetype","clonotype_number")]
  P_count <- P %>% distinct(clonetype, .keep_all = TRUE)
  colnames(P_count) <- c("clonotype","clonotype_number")
  P_count$sampleID <- each
  P_count <- P_count[order(P_count$clonotype_number,decreasing = TRUE),]
  if(nrow(P_count) < 3){
    print(nrow(P_count))
    clonotype <- as.vector(P_count$clonotype)
    clonotype_number <- as.vector(P_count$clonotype_number)
    sampleID <- as.vector(P_count$sampleID)
    for(i in 1:(3-nrow(P_count))){
      print(i)
      clonotype <- c(clonotype,paste0(each,"_new_clone_",i))
      clonotype_number <- c(clonotype_number,0)
      sampleID <- c(sampleID, each)
      P_count <- as.data.frame(cbind(clonotype,clonotype_number,sampleID))
    }
  }else{
    P_count <- P_count[1:3,]
  }
  print(P_count)
  if(count == 0){
    df <- P_count
    count <- count + 1
  }else{
    df <- rbind(df, P_count)
  }
}

head(df)

for(each in cluster.info$sampleID){
  if(! each %in% unique(df$sampleID)){
    clonotype <- c(paste0(each,"_new_clone_1"),paste0(each,"_new_clone_2"),paste0(each,"_new_clone_3"))
    clonotype_number <- c(0,0,0)
    sampleID <- c(each,each,each)
    P_count <- as.data.frame(cbind(clonotype,clonotype_number,sampleID))
    print(P_count)
    df <- rbind(df, P_count)
  }
}

length(unique(df$sampleID))


df <- merge(df, cluster.info, by = "sampleID", all.x = TRUE)
head(df)
df$group <- paste0("group",df$group)
head(df)

df <- merge(df, sample.info, by = "sampleID", all.x = TRUE)
head(df)

df$type <- paste0(df$group,"_",df$pathology)
head(df)

df$type <- factor(as.vector(df$type), levels = c("group1_MPR","group1_non-MPR",
                                                     "group2_MPR","group2_non-MPR",
                                                     "group3_MPR","group3_non-MPR",
                                                     "group4_MPR","group4_non-MPR",
                                                     "group5_MPR","group5_non-MPR"))

compaired <- list(c("group3_non-MPR", "group3_MPR"))
df$clonotype_number <- as.numeric(df$clonotype_number)
ggboxplot(df, x = "type", y = "clonotype_number",
          color = "type",add="jitter",
          x.text.angle=90) + labs(x='group', y= 'number of cells in top3 Tex-relevant clonotypes') + 
  theme(legend.position="none") + 
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)

df$pathology <- as.factor(unlist(df$pathology))
df$group <- factor(df$group,levels = c("group1","group2","group3","group4","group5"))
ggboxplot(df, x = "group", y = "clonotype_number",
          color = "pathology",add = "jitter",
          x.text.angle=0,size = 0.5,add.params=list(size=0.6)) + 
  scale_color_manual(values=c("MPR" = "#2868A6", "non-MPR" = "#B1161C")) + 
  labs(y= 'number of cells in top3 Tex-relevant clonotypes')
