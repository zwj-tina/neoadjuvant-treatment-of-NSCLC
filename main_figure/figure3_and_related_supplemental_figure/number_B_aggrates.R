library(readxl)
#只保留了T1的数目
B.info <- as.data.frame(read_excel("B_cell_aggregates.xlsx"))
head(B.info)

cluster.info <- read.csv("NMF_all_group_5.csv")
cluster.info <- cluster.info[,-1]
head(cluster.info)

df <- merge(cluster.info, B.info, by = "sampleID", all.x = TRUE)
head(df)
df <- df[df$group %in% c("1","2","3","4","5"),]
head(df)
df[is.na(df)] <- "unknown"
head(df)
df <- df[df$number_of_B_cell_aggregates != "unknown",]
dim(df)
df$number_of_B_cell_aggregates <- as.numeric(df$number_of_B_cell_aggregates)
df$number_of_B_cell_aggregates <- 4 * df$number_of_B_cell_aggregates
df$group <- paste0("group", df$group)
df$group <- factor(as.vector(df$group), levels = c("group1","group2","group3",
                                                   "group4","group5"))

compaired <- list(c("group1", "group2"),c("group2", "group3"),
                  c("group2", "group4"),c("group2", "group5"))

ggboxplot(df, x = "group", y = "number_of_B_cell_aggregates",
          color = "group",add="jitter",add.params=list(size=0.7),
          x.text.angle=0) + labs(x='group', y= 'number_of_B_cell_aggregates / cm2') + 
  theme(legend.position="none") +
  scale_color_manual(values=c("group1"="#E84C35","group2"="#4FBAD6",
                              "group3"="#00A289","group4"="#3C5487",
                              "group5"="#F29B80"),guide = "none") +  
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test) 


new_group <- c()
for(each in df$group){
  if(each %in% c("group1","group3","group4","group5")){
    new_group <- c(new_group, "other_group")
  }else{
    new_group <- c(new_group, "group2")
  }
}

df$new_group <- new_group
df$new_group <- factor(df$new_group, levels = c("group2","other_group"))

compaired <- list(c("group2","other_group"))
ggboxplot(df, x = "new_group", y = "number_of_B_cell_aggregates",
          color = "new_group",add="jitter",add.params=list(size=0.7),
          x.text.angle=0) + labs(x='group', y= 'number_of_B_cell_aggregates / cm2') + 
  theme(legend.position="none") +
  scale_color_manual(values=c("other_group"="#8491B4CC","group2"="#4FBAD6"),guide = "none") +  
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test) 


