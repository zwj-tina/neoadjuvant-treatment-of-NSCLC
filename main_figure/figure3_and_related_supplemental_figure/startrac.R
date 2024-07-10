library(Startrac)
library(ggpubr)
library(ggplot2)
library(circlize)
library(ggpmisc)
library(ggsci)

####################################################################
TCR.data <- read.csv("T_with_TCR.csv")
# TCR.data must include the clone type, expansion state and the cluster
head(TCR.data)
dim(TCR.data)

#only include CD8 T cells
TCR.data <- TCR.data[TCR.data$sub_cell_type %in% c("your sub cell types"),]
head(TCR.data)

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
expan$majorCluster <- factor(expan$majorCluster, levels = c("your sub cell types"))

compaired <- list(c(),c())
ggboxplot(expan, x = "majorCluster", y = "expa",
          color = "majorCluster",add="jitter",add.params=list(size=0.5),
          x.text.angle=45) + labs(x='cell.type', y= 'expa') + 
  theme(legend.position="none")  +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)

