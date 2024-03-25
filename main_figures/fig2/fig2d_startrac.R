library(Startrac)
library(tictoc)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(circlize)
library(ggpmisc)
library(ggsci)

####################################################################
TCR.data <- read.csv("/home/zwj/raid1/neoadjuvant/data/process_data/T_with_TCR_obs_V2.csv")
head(TCR.data)

#only include CD8 T cells
TCR.data <- TCR.data[TCR.data$sub_cell_type %in% c("CD8T_Tem_GZMK+GZMH+","CD8T_Tm_IL7R",
                                                   "CD8T_Tex_CXCL13","CD8T_Tem_GZMK+NR4A1+",
                                                   "CD8T_NK-like_FGFBP2","CD8T_Trm_ZNF683",
                                                   "CD8T_ISG15","CD8T_terminal_Tex_LAYN",
                                                   "CD8T_MAIT_KLRB1","CD8T_prf_MKI67"),]
head(TCR.data)

cluster.info <- read.csv("NMF_LUSC_group_5.csv")
cluster.info <- cluster.info[,-1]
head(cluster.info)

TCR.data <- merge(TCR.data, cluster.info, by = "sampleID", all.x = TRUE)
head(TCR.data)
TCR.data <- TCR.data[TCR.data$group %in% c("2","3","4","5"),]
head(TCR.data)

in.dat <- TCR.data[,c("sampleID","cellID","clonetype","expansion",
                      "sub_cell_type")]
head(in.dat)

colnames(in.dat) <- c("patient","Cell_Name","clone.id","clone.status",
                      "majorCluster")
head(in.dat)
in.dat$loc = "T"

head(in.dat)

out <- Startrac.run(in.dat, proj="NSCLC",verbose=F)

expan <- out@cluster.data
head(expan)

expan$majorCluster <- factor(expan$majorCluster, levels = c("CD8T_ISG15",
                                                            "CD8T_MAIT_KLRB1",
                                                            "CD8T_prf_MKI67",
                                                            "CD8T_Tem_GZMK+NR4A1+",
                                                            "CD8T_Tm_IL7R",
                                                            "CD8T_Tem_GZMK+GZMH+",
                                                            "CD8T_NK-like_FGFBP2",
                                                            "CD8T_terminal_Tex_LAYN",
                                                            "CD8T_Trm_ZNF683",
                                                            "CD8T_Tex_CXCL13"))

head(expan)

mycol <- paletteer_d("ggthemes::Classic_20", n=20)
compaired <- list(c("CD8T_NK-like_FGFBP2", "CD8T_Tex_CXCL13"))

ggboxplot(expan, x = "majorCluster", y = "expa",
          color = "majorCluster",add="jitter",
          x.text.angle=45) + labs(x='cell.type', y= 'STARTRAC-expa') + 
  theme(legend.position="none") + 
  scale_color_npg() + scale_fill_npg() + 
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)

