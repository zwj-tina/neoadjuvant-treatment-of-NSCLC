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
head(TCR.data)
a <- table(TCR.data$clonetype,TCR.data$T_new_name)
a <- a / rowSums(a)
a <- as.data.frame(a)
colnames(a) <- c("clonotype","cell.type","Freq")
a <- a[a$cell.type %in% c("expanded_terminal_Tex"),]
head(a)

df <- TCR.data[,c("sampleID","clonetype","clonetype_number")]
colnames(df) <- c("sampleID","clonotype","clonotype_number")
df <- df %>% distinct(clonotype, .keep_all = TRUE)
dim(df)

a <- merge(a,df,by = "clonotype")
head(a)


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

sample.info <- sample.info[,c("sampleID","pathology")]
head(sample.info)

a <- merge(a, sample.info,by = "sampleID")
head(a)

a <- merge(a,cluster.info, by = "sampleID")
head(a)
a$group <- paste0("group",a$group)
head(a)


Treg <- read.csv("expanded_CD4Treg_clonotype_number.csv")
Treg <- Treg[,-1]
colnames(Treg) <- c("sampleID","number")
rownames(Treg) <- Treg$sampleID
head(Treg)
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
dim(Treg)

a <- merge(a,Treg,by = "sampleID")
a

cluster <- a[a$group %in% c("group1"),]
cluster$clonotype_number <- as.numeric(cluster$clonotype_number)
cluster$Freq <- as.numeric(cluster$Freq)
cluster$pathology <- as.factor(unlist(cluster$pathology))

ggscatter(cluster, 
          x = "clonotype_number",
          y = "Freq", color = "pathology",shape = "Treg_level")+
  scale_color_manual(values=c("MPR" = "#2868A6", "non-MPR" = "#B1161C"))+ ylim(0,1)+xlim(0,400)
