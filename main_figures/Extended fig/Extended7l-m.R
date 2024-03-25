library(readxl)
library(ggalluvial)
sample.info <- as.data.frame(read_excel("sample.xlsx"))
head(sample.info)

cluster.info <- read.csv("NMF_LUAD_group_5.csv")
cluster.info <- cluster.info[,-1]
head(cluster.info)

sample.info <- merge(sample.info, cluster.info, by = "sampleID", all.x = TRUE)
head(sample.info)
sample.info <- sample.info[sample.info$group %in% c("1","2","3","4","5"),]
head(sample.info)

sample.info$group <- paste0("group",sample.info$group)


sample.info <- sample.info[sample.info$pathological_response %in% c("pCR","MPR","pPR","nPR"),]
response <- list(
  "pCR" = "MPR",
  "MPR" = "MPR",
  "pPR" = "pPR",
  "nPR" = "nPR"
)
sample.info$pathological_response <- unlist(sapply(as.vector(sample.info$pathological_response), function(x) response[[x]]))
sample.info$pathological_response <- factor(sample.info$pathological_response, levels = c("MPR","pPR","nPR"))

sample.info$cluster <- paste0(sample.info$group,"_",sample.info$pathological_response)
head(sample.info)


a <- table(sample.info$pathological_response, sample.info$group)
a <- as.data.frame(a / rowSums(a))
colnames(a) <- c("response","group","Freq")
head(a)
ggbarplot(a, x="response", y="Freq", fill = "group",
          x.text.angle=0) + theme(legend.position = "right") + 
  scale_fill_manual(values=c("group1"="#E84C35","group2"="#4FBAD6",
                             "group3"="#00A289","group4"="#3C5487",
                             "group5"="#F29B80")) 

a <- as.data.frame(table(sample.info$group,sample.info$pathological_response))
colnames(a) <- c("group","response","number")
head(a)
ggbarplot(a, x="group", y="number", fill = "response",
          x.text.angle=0) + theme(legend.position = "right") + 
  scale_fill_manual(values=c("nPR" = "#EDB752","pPR" = "#FFC1B3","MPR" = "#AFD483")) 


#PD1
PDL1_group <- c()
for(each in sample.info$PDL1_TPS){
  if (each %in% c("<1%","0")){
    PDL1_group <- c(PDL1_group,"<1%") 
  }else if (each %in% c("0.02","0.05","0.08","0.03","0.01","0.3","0.4","0.2","0.35")){
    PDL1_group <- c(PDL1_group,"1%-49%") 
  }else if (each %in% c("0.5","0.7","0.6","0.9","0.55","0.8","0.85","1","0.75","0.65")){
    PDL1_group <- c(PDL1_group,">=50%") 
  }else{
    PDL1_group <- c(PDL1_group,"Not tested") 
  }
}
sample.info$PDL1_group <- PDL1_group

sample.info <- sample.info[sample.info$PDL1_group != "Not tested",]
sample.info$PDL1_group <- factor(as.vector(sample.info$PDL1_group), levels = c("<1%","1%-49%",">=50%"))

a <- table(sample.info$cluster,sample.info$PDL1_group)
df <- as.data.frame(a)
head(df)
colnames(df) <- c("cluster","PDL1_TPS","number")
head(df)

df$cluster <- factor(df$cluster, levels = c(
  "group1_MPR","group1_nPR",
  "group2_MPR","group2_nPR",
  "group3_MPR","group3_pPR","group3_nPR",
  "group4_pPR","group4_nPR",
  "group5_MPR","group5_nPR"
))


df$group <- sapply(as.vector(df$cluster), function(x) strsplit(x,"_")[[1]][1])
head(df)
df$pathological_response <- sapply(as.vector(df$cluster), function(x) strsplit(x,"_")[[1]][2])
head(df)
df <- df[!df$number %in% c(0),]
df
df$pathological_response <- factor(as.vector(df$pathological_response), levels = c("MPR","pPR","nPR"))


ggplot(data = df,
       aes(axis1 = PDL1_TPS,   # First variable on the X-axis
           axis2 = group,   # Third variable on the X-axis
           y = number)) +
  geom_alluvium(aes(fill = pathological_response,order = pathological_response)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void() +
  scale_fill_manual(values=c("nPR" = "#EDB752","pPR" = "#FFC1B3","MPR" = "#AFD483"))


a <- table(sample.info$cluster,sample.info$PDL1_group)
df <- as.data.frame(a)
head(df)
colnames(df) <- c("cluster","PDL1_TPS","number")
head(df)
df <- df[!df$number %in% c(0),]
ggbarplot(df, x="cluster", y="number", fill = "PDL1_TPS",
          x.text.angle=90) + theme(legend.position = "right") + 
  scale_fill_manual(values=c("<1%"="#B3D9D9","1%-49%"="#4F9D9D",
                             ">=50%"="#3D7878"))

group <- c()
response <- c()
for(each in df$cluster){
  group <- c(group, str_split(each, "_")[[1]][1])
  response <- c(response, str_split(each, "_")[[1]][2])
}
df$group <- group
df$response <- response
head(df)
df$response <- factor(df$response, levels = c("MPR","pPR","nPR"))
ggbarplot(df, x="PDL1_TPS", y="number", fill = "group",
          x.text.angle=90,facet.by = "response", ncol = 4) + theme(legend.position = "right") + 
  scale_fill_manual(values=c("group1"="#E84C35","group2"="#4FBAD6",
                             "group3"="#00A289","group4"="#3C5487",
                             "group5"="#F29B80"))
