library(ggpmisc)
library(ggsci)
library(paletteer)
library(readxl)
#############################################################################################################
#number of Tex-relevant clonotype of LUSC
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

#only remain clonetype
df <- TCR.data[,c("sampleID","clonetype")]
df = df[!duplicated(df$clonetype),]

count.clonetype <- as.data.frame(table(df$sampleID))
colnames(count.clonetype) <- c("sampleID","number")
head(count.clonetype)

# fill 0 for samples with no Tex-relevant clones
sampleID <- as.vector(count.clonetype$sampleID)
number <- as.vector(count.clonetype$number)
for(each in cluster.info$sampleID){
  if(!each %in% count.clonetype$sampleID){
    print(each)
    sampleID <- c(sampleID,each)
    number <- c(number,0)
  }
}

df_number <- as.data.frame(cbind(sampleID, number))
df_number$number <- as.numeric(df_number$number)
head(df_number)

df_number <- merge(df_number, cluster.info, by = "sampleID", all.x = TRUE)
df_number$group <- paste0("group",df_number$group)
head(df_number)

sample.info <- as.data.frame(read_excel("sample.xlsx"))
sample.info <- sample.info[,-c(2,3,4,5,6,7,8,9,10,11,12,13,14,18,19,20,21)]
df_number <- merge(df_number,sample.info, by = "sampleID", all.x = TRUE)
head(df_number)

response <- list(
  "pCR" = "pCR",
  "MPR" = "MPR",
  "pPR" = "non-MPR",
  "nPR" = "non-MPR",
  "non-MPR" = "non-MPR"
)
df_number$pathology <- sapply(as.vector(df_number$pathological_response), function(x) response[[x]])
head(df_number)

#############################################################################################################
#number of Treg clonotype of LUSC
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

#only include expanded CCR8 Treg cells
TCR.data <- TCR.data[TCR.data$T_new_name %in% c("expanded_CCR8_Treg"),]
dim(TCR.data)
df <- TCR.data[,c("sampleID","clonetype")]
df = df[!duplicated(df$clonetype),]

count.clonetype <- as.data.frame(table(df$sampleID))
colnames(count.clonetype) <- c("sampleID","number")
head(count.clonetype)

# fill 0 for samples with no activated Treg clones
sampleID <- as.vector(count.clonetype$sampleID)
number <- as.vector(count.clonetype$number)
for(each in cluster.info$sampleID){
  if(!each %in% count.clonetype$sampleID){
    print(each)
    sampleID <- c(sampleID,each)
    number <- c(number,0)
  }
}

df_number <- as.data.frame(cbind(sampleID, number))
head(df_number)
df_number$number <- as.numeric(df_number$number)
head(df_number)

df_number <- merge(df_number, cluster.info, by = "sampleID", all.x = TRUE)
df_number$group <- paste0("group",df_number$group)
head(df_number)

library(readxl)
sample.info <- as.data.frame(read_excel("sample.xlsx"))
df_number <- merge(df_number,sample.info, by = "sampleID", all.x = TRUE)
head(df_number)

response <- list(
  "pCR" = "MPR",
  "MPR" = "MPR",
  "pPR" = "non-MPR",
  "nPR" = "non-MPR",
  "non-MPR" = "non-MPR"
)

df_number$pathology <- sapply(as.vector(df_number$pathological_response), function(x) response[[x]])
head(df_number)

df_number$type <- paste0(df_number$group,"_",df_number$pathology)
head(df_number)
