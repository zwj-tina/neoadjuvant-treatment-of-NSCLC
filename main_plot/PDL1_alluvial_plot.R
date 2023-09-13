library(readxl)
library(ggplot2)
sample.info <- as.data.frame(read_excel("../main_data/sample.xlsx"))
head(sample.info)

cluster.info <- read.csv("../main_data/NMF_LUSC_group_5.csv")
cluster.info <- cluster.info[,-1]
head(cluster.info)

sample.info <- merge(sample.info, cluster.info, by = "sampleID", all.x = TRUE)
head(sample.info)
sample.info <- sample.info[sample.info$group %in% c("1","2","3","4","5"),]
head(sample.info)

sample.info$group <- paste0("group",sample.info$group)

pathological_response <- c()
for(each in sample.info$pathological_response){
  if(each %in% c("MPR","pCR")){
    pathological_response <- c(pathological_response, "MPR")
  }else{
    pathological_response <- c(pathological_response, "non-MPR")
  }
}

sample.info$pathological_response <- pathological_response
sample.info$cluster <- paste0(sample.info$group,"_",sample.info$pathological_response)

#PD1
PDL1_TPS_group <- c()
for(each in sample.info$PDL1_TPS){
  if (each %in% c("<1%","0")){
    PDL1_TPS_group <- c(PDL1_TPS_group,"<1%") 
  }else if (each %in% c("0.01","0.02","0.05","0.08","0.03","0.3","0.2","0.35","0.4")){
    PDL1_TPS_group <- c(PDL1_TPS_group,"1%-49%") 
  }else if (each %in% c("0.5","0.7","0.6","0.9","0.8","0.85","1","0.55","0.75","0.65")){
    PDL1_TPS_group <- c(PDL1_TPS_group,">=50%") 
  }else{
    PDL1_TPS_group <- c(PDL1_TPS_group,"Not tested") 
  }
}
sample.info$PDL1_TPS_group <- PDL1_TPS_group

sample.info <- sample.info[sample.info$PDL1_TPS_group != "Not tested",]
sample.info$PDL1_TPS_group <- factor(as.vector(sample.info$PDL1_TPS_group), levels = c("<1%","1%-49%",">=50%"))
a <- table(sample.info$cluster,sample.info$PDL1_TPS_group)
df <- as.data.frame(a)
head(df)
colnames(df) <- c("cluster","PDL1_TPS","number")
head(df)
df$cluster <- factor(df$cluster,
                     levels = c("group1_MPR","group1_non-MPR",
                                "group2_MPR","group2_non-MPR",
                                "group3_MPR","group3_non-MPR",
                                "group4_MPR","group4_non-MPR",
                                "group5_MPR","group5_non-MPR"))

head(df)

df$group <- sapply(as.vector(df$cluster), function(x) strsplit(x,"_")[[1]][1])
head(df)
df$pathological_response <- sapply(as.vector(df$cluster), function(x) strsplit(x,"_")[[1]][2])
head(df)
df <- df[!df$number %in% c(0),]
df
ggplot(data = df,
       aes(axis1 = PDL1_TPS,   # First variable on the X-axis
           axis2 = group,   # Third variable on the X-axis
           y = number)) +
  geom_alluvium(aes(fill = pathological_response,order = pathological_response)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void() +
  scale_fill_manual(values=c("MPR" = "#ffcf5e", "non-MPR" = "#4E6EE8"))


