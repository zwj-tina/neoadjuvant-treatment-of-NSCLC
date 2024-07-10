Tex_relevant <- read.csv("CD8Tex_relevant_clonotype_number_over2.csv")
Tex_relevant <- Tex_relevant[,-1]
colnames(Tex_relevant) <- c("sampleID","number")
rownames(Tex_relevant) <- Tex_relevant$sampleID

all.group <- read.csv("/home/zhangwj/data_yi/neoadjuvant/revision2/data/NMF_all_group_5.csv")
all.group <- all.group[,-1]
rownames(all.group) <- all.group$sampleID

samples <- intersect(Tex_relevant$sampleID, all.group$sampleID)
length(samples)

Tex_relevant <- Tex_relevant[samples,]
all.group <- all.group[samples,]

Tex_relevant$group <- all.group$group
head(Tex_relevant)

library(readxl)
sample.info <- as.data.frame(read_excel("sample.xlsx"))
sample.info <- sample.info[,c("sampleID","pathological_response","cancer_type")]
head(sample.info)
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
sample.info <- sample.info[rownames(Tex_relevant),]

Tex_relevant <- merge(Tex_relevant,sample.info,by = "sampleID")
head(Tex_relevant)
Tex_relevant$group <- paste0("group",Tex_relevant$group)
colnames(Tex_relevant) <- c("sampleID","Tex_number","group","pathological_response_detail","cancer_type","pathological_response")

Treg <- read.csv("expanded_CD4Treg_clonotype_number_over2.csv")
Treg <- Treg[,-1]
colnames(Treg) <- c("sampleID","number")
rownames(Treg) <- Treg$sampleID
head(Treg)
dim(Treg)
Treg <- Treg[Tex_relevant$sampleID,]
head(Treg)

Tex_relevant$Treg_number <- Treg$number
head(Tex_relevant)

ggscatter(Tex_relevant[Tex_relevant$pathological_response %in% c("non-MPR"),], 
          x = "Tex_number",
          y = "Treg_number", color = "group",shape = "cancer_type")+
  geom_hline(yintercept = 9,linetype = "dashed") + 
  scale_color_manual(values=c("group1"="#E84C35","group2"="#4FBAD6",
                              "group3"="#00A289","group4"="#3C5487",
                              "group5"="#F29B80")) + geom_abline(intercept = 0, slope = 1,linetype = "dashed") +
  xlim(0,100) +ylim(0,100)


ggscatter(Tex_relevant[Tex_relevant$pathological_response %in% c("MPR"),], 
          x = "Tex_number",
          y = "Treg_number", color = "group",shape = "cancer_type")+
  geom_hline(yintercept = 9,linetype = "dashed") + 
  scale_color_manual(values=c("group1"="#E84C35","group2"="#4FBAD6",
                              "group3"="#00A289","group4"="#3C5487",
                              "group5"="#F29B80")) + geom_abline(intercept = 0, slope = 1,linetype = "dashed") +
  xlim(0,100) +ylim(0,100)
