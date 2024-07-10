########################################################################################
#load library
#########################################################################################
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(readr))
suppressMessages(library(harmony))
suppressMessages(library(ggpubr))
suppressMessages(library(ggpmisc))
suppressMessages(library(ggrepel))
suppressMessages(library(readxl))
##############################################################################################
#B cell
##########################################################################
scRNA <- readRDS("Treg_count.rds")
dim(scRNA)
Idents(scRNA) <- "sub_cell_type"

scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

scRNA <- FindVariableFeatures(scRNA, selection.method = "vst",nfeatures = 1000)
length(scRNA@assays$RNA@var.features)

all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))

ElbowPlot(scRNA)

Idents(scRNA) <- "sub_cell_type"
cluster1.markers <- FindMarkers(scRNA, ident.1 = "CD4T_Treg_CCR8",min.pct = 0.5,
                                test.use = "wilcox_limma",slot = "counts")
head(cluster1.markers, n = 5)

data.markers <- cluster1.markers

data.markers$symbol <- rownames(data.markers)
data.markers$logP <- -log10(data.markers$p_val_adj + 1e-100)
dim(data.markers)
data.markers$Group = "not-significant"
data.markers$Group[which((data.markers$p_val_adj < 0.05) & (data.markers$avg_log2FC > 1.4))] = "CD4T_Treg_CCR8"
data.markers$Group[which((data.markers$p_val_adj < 0.05) & (data.markers$avg_log2FC < -0.9))] = "CD4T_Treg_FOXP3"
table(data.markers$Group)

data.markers$label = ""
#对差异基因的p值进行从小到大的排序
data.markers <- data.markers[order(data.markers$avg_log2FC, decreasing = TRUE),]
#高表达基因中选取p_val_adj最小的10个
up.genes <- head(data.markers$symbol[which(data.markers$Group == "CD4T_Treg_CCR8")], 13)
#低表达基因中选取p_val_adj最小的10个
down.genes <- tail(data.markers$symbol[which(data.markers$Group == "CD4T_Treg_FOXP3")], 7)
#将up.genes和down.genes合并并加入到Label
data.top10.genes <- c(as.character(up.genes), as.character(down.genes))
data.markers$label[match(data.top10.genes, data.markers$symbol)] <- data.top10.genes
ggscatter(data.markers, x = "avg_log2FC", y = "logP", color = "Group",
          palette = c("#CC0000","#2f5688","#BBBBBB"),
          size = 1, font.label = 18,
          repel = T, xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)") + 
  geom_text_repel(size=3,point.padding = NA,label = data.markers$label, max.overlaps = 1000)
