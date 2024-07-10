#downloaded data from GSE169246
TNBC <- readRDS("/home/zhangwj/data_yi/neoadjuvant/data/other_data/zhangyy/zyy_TNBC.rds")

dim(TNBC)
Idents(TNBC) <- "cellType_in_paper"
DimPlot(TNBC, reduction = "umap", label = FALSE,pt.size = 0.1) +NoLegend()

FeaturePlot(TNBC, features = c("GBP1"), cols = c("lightgrey" ,"#FD3131"),pt.size = 0.1) 
FeaturePlot(TNBC, features = c("NKG7"), cols = c("lightgrey" ,"#FD3131"))

NKT.part <- subset(TNBC, cellType_in_paper %in% c("t_Tn-LEF1","t_ILC1-IL32","t_CD8_Tem-GZMK",
                                                  "t_CD4-CXCL13","t_ILC1-GZMK",
                                                  "t_CD8_MAIT-KLRB1",
                                                  "t_CD4_Treg-FOXP3",
                                                  "t_CD8_Trm-ZNF683","t_CD8_Teff-GNLY","t_CD4_Tcm-LMNA",   
                                                  "t_ILC3-AREG",     
                                                  "t_CD8-CXCL13","t_ILC1-IFNG",      
                                                  "t_ILC1-FGFBP2","t_Tact-IFI6",     
                                                  "t_ILC1-ZNF683",   
                                                  "t_ILC1-CD160","t_ILC3-IL7R","t_ILC1-CX3CR1",    
                                                  "t_ILC1-SELL","t_CD4_Tact-XIST","t_ILC2-SPON2",    
                                                  "t_ILC1-CNOT2","t_Tprf-MKI67",   
                                                  "t_ILC1-VCAM1"))
dim(NKT.part)
NKT.part <- NormalizeData(NKT.part, normalization.method = "LogNormalize", scale.factor = 10000)

NKT.part <- FindVariableFeatures(NKT.part, selection.method = "vst",nfeatures = 1000)
#delet IgG/H/L
NKT.part@assays$RNA@var.features <- NKT.part@assays$RNA@var.features[-which(NKT.part@assays$RNA@var.features %in% grep("^IG[KHL]",NKT.part@assays$RNA@var.features,value=T))]
NKT.part@assays$RNA@var.features <- NKT.part@assays$RNA@var.features[-which(NKT.part@assays$RNA@var.features %in% grep("^MT",NKT.part@assays$RNA@var.features,value=T))]
#NKT.part@assays$RNA@var.features <- NKT.part@assays$RNA@var.features[-which(NKT.part@assays$RNA@var.features %in% grep("^RP[LS]",NKT.part@assays$RNA@var.features,value=T))]
length(NKT.part@assays$RNA@var.features)

all.genes <- rownames(NKT.part)
NKT.part <- ScaleData(NKT.part, features = all.genes)
NKT.part <- RunPCA(NKT.part, features = VariableFeatures(object = NKT.part))

ElbowPlot(NKT.part)

#NKT.part <- subset(NKT.part, sample %in% c("medial 2","distal 2",
#                                           "distal 1a","proximal 3",
#                                           "distal 3"))
#remove batch effect
NKT.part <- RunHarmony(NKT.part, c("sampleID"))

NKT.part <- NKT.part %>% 
  FindNeighbors(reduction = "harmony", dims = 1:15) %>% 
  FindClusters(resolution = 0.5) %>% 
  RunUMAP(reduction = "harmony", dims = 1:15) %>% 
  identity()

NKT.part <- FindClusters(NKT.part,resolution = 1)

DimPlot(NKT.part, reduction = "umap", label = TRUE,pt.size = 1) 
NKT.part$cellType_in_paper

FeaturePlot(NKT.part, features = c(""),pt.size = 1, cols = c("lightgrey" ,"#FD3131")) 

NKT.markers <- FindAllMarkers(NKT.part, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
NKT.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> NKT.top5

FeaturePlot(NKT.part, features = c("FCGR3A"), cols = c("lightgrey" ,"#FD3131")) 
FeaturePlot(NKT.part, features = c("FGFBP2"), cols = c("lightgrey" ,"#FD3131"))

DimPlot(NKT.part, reduction = "umap", label = TRUE,pt.size = 0.1) + NoLegend()

Idents(NKT.part) <- "seurat_clusters"
new.cluster.ids <- c("other","other","other","other","other","other",
                     "other","other","other","NK_CD16hi_FGFBP2","other",
                     "other","other","other","other","other","NK_CD16hi_FGFBP2",
                     "other")
names(new.cluster.ids) <- levels(NKT.part)
NKT.part <- RenameIdents(NKT.part, new.cluster.ids)
DimPlot(NKT.part, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() +
  scale_color_manual(values = c("#8ECFC9","#FA7F6F"))
NKT.part$new.cell.type <- Idents(NKT.part)

a <- table(NKT.part$sampleID,NKT.part$new.cell.type)
a <- as.data.frame(a / rowSums(a))
colnames(a) <- c("sampleID","cell.type","Freq")
head(a)

obs <- NKT.part@meta.data
obs <- obs[,c("sampleID","patientID","tissue","treatment_status",
                    "ICB_treatment","treatment","treatment_response")]
obs <- obs %>% distinct(sampleID, .keep_all = TRUE)
head(obs)

df <- merge(a, obs, by = "sampleID", all.x = TRUE)
head(df)
df$treatment_response <- factor(df$treatment_response, levels = c("PD","SD","PR"))

df <- df[df$cell.type == "NK_CD16hi_FGFBP2",]
df <- df[df$treatment %in% c("Chemo"),]
df <- df[df$treatment_status %in% c("Post-treatment"),]

treatment.response <- c()
for(each in df$treatment_response){
  if(each %in% c("PD","SD")){
    treatment.response <- c(treatment.response,"PD/SD")
  }else{
    treatment.response <- c(treatment.response,"PR")
  }
}
df$treatment.response <- treatment.response
compaired <- list(c("PD/SD","PR"))

ggboxplot(df, x = "treatment.response", y = "Freq",
          color = "treatment.response",add = "jitter",
          x.text.angle=0,size = 0.5,pt.size = 1, facet.by = "treatment") +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test) +
  scale_color_manual(values = c("PD/SD" = "#88AB8E","PR"="#E97777")) + theme(legend.position="none")
