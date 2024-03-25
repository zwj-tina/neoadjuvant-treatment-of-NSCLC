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
head(sample.info)

#only include top3 clonopes and fill 0 if the sample don't have 3 clonotypes
count = 0
for(each in unique(TCR.data$sampleID)){
  P <- TCR.data[TCR.data$sampleID == each,c("clonotype","clonotype_number")]
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
df$pathology <- as.factor(unlist(df$pathology))

Treg <- read.csv("expanded_CD4Treg_clonotype_number.csv")
Treg <- Treg[,-1]
colnames(Treg) <- c("sampleID","number")
rownames(Treg) <- Treg$sampleID
head(Treg)
dim(Treg)
Treg <- Treg[Treg$sampleID %in% df$sampleID,]
dim(Treg)
Treg_level <- c()
for(each in Treg$number){
  if(each >= 4){
    Treg_level <- c(Treg_level,"high")
  }else{
    Treg_level <- c(Treg_level,"low")
  }
}
Treg$Treg_level <- Treg_level
head(Treg)


df <- merge(df, Treg, by = "sampleID", all.x = TRUE)
head(df)

df$type <- paste0(df$pathology,"_",df$Treg_level)
head(df)

df <- df[df$pathology %in% c("non-MPR"),]
head(df)

table(df$pathology)

df$type <- factor(df$type, levels = c("non-MPR_low","non-MPR_high"))
df$clonotype_number <- as.numeric(df$clonotype_number)

compaired <- list(c("non-MPR_low","non-MPR_high"))
ggboxplot(df, x = "type", y = "clonotype_number",
          color = "type",add="jitter",add.params=list(size=0.6),
          x.text.angle=0) + labs(x='group', y= 'clone number of each clonotype in Tex-relevant cells') + 
  theme(legend.position="none")  + 
  scale_color_manual(values=c("#D9BFAE","#8CB4A3")) + 
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test) 

