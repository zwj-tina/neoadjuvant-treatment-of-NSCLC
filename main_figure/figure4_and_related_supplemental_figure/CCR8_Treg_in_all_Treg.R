#############################################################################################################
info <- read.csv("all_sub_cell_type.csv")
head(info)
length(unique(info$sampleID))

info <- info[info$sub_cell_type %in% c("CD4T_Treg_CCR8","CD4T_Treg_FOXP3","CD4T_Treg_MKI67"),]

cluster.info <- read.csv("NMF_all_group_5.csv")
cluster.info <- cluster.info[,-1]
head(cluster.info)
dim(cluster.info)


info <- info[info$sampleID %in% cluster.info$sampleID,]
length(unique(info$sampleID))


df <- table(info$sampleID,info$sub_cell_type)
ratio <- as.data.frame(df / rowSums(df))
head(ratio)
colnames(ratio) <- c("sampleID","cell.type","Freq")
head(ratio)

ratio <- merge(ratio,cluster.info)
head(ratio)
ratio$group <- paste0("group",ratio$group)
ratio$group <- factor(ratio$group, levels = c("group1","group2","group3","group4","group5"))
head(ratio)

CCR8 <- ratio[ratio$cell.type %in% c("CD4T_Treg_CCR8"),]

compaired <- list(c("group3","group4"))

ggboxplot(CCR8[CCR8$group %in% c("group3","group4"),], x = "group", y = "Freq",
          color = "group",add="jitter",add.params=list(size=0.5),
          x.text.angle=0) + labs(x='group', y= 'CCR8+ Treg in all Tregs') +
  theme(legend.position="none") +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test) + 
  scale_color_manual(values=c("group3"="#00A289","group4"="#3C5487"))
                                                                                                                              "group3"="#00A289","group4"="#3C5487",
