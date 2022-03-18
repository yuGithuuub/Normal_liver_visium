# STutility
library("STutility")
library(harmony)
library(Seurat)
prepare_df <- read.csv("./data_prepare.csv") #
prepare_df
infoTable <- prepare_df
se <- InputFromTable(infotable = infoTable, 
                     min.gene.count = 100, 
                     min.gene.spots = 5,
                     min.spot.count = 500,
                     platform =  "Visium",verbose = T)

ST.FeaturePlot(se, features = c("nFeature_RNA"), cols = c("lightgray", "mistyrose", "red", "darkred", "black"), ncol = 4, pt.size = 1,label.by = "typeF",pt.alpha = 0.5)
se@meta.data
se@active.ident <- as.factor(se$orig)

se$orig
se <- SCTransform(se, method = "glmGamPoi", vars.to.regress ="orig")
save(se,file = "../object/visium_all_sctransform.rdata")

se <- RunPCA(se, assay = "SCT", verbose = FALSE)
se <- RunHarmony(se,reduction = "pca",dims.use = c(1:30),group.by.vars = "orig",assay.use = "SCT")
se <- RunUMAP(se, reduction = "harmony", dims = 1:30)
se <- FindNeighbors(se, reduction = "harmony", dims = 1:30)
se <- FindClusters(se, verbose = FALSE,resolution = 0.5)

p1 <- DimPlot(se, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(se, label = TRUE, label.size = 3)
p1 + p2


DimPlot(se, reduction = "umap", label = TRUE)
DimPlot(se, reduction = "umap", label = TRUE, split.by = "orig")

MARKER_se <- FindAllMarkers(se,only.pos = T)

se_L5 <- SubsetSTData(se,expression = orig %in% c("L5"))

tablel5_umap <- data.frame(barcode = substr(colnames(se_L5),start = 1,stop = 18),UMAP_1 = se_L5@reductions$umap@cell.embeddings[,"UMAP_1"],UMAP_2= se_L5@reductions$umap@cell.embeddings[,"UMAP_2"])
tablel5_cluster <- data.frame(barcode = substr(colnames(se_L5),start = 1,stop = 18),category = se_L5@active.ident)

rownames(tablel5_umap) <- NULL
rownames(tablel5_cluster) <- NULL

write.csv(tablel5_umap,file = "l5_umap.csv")
write.csv(tablel5_cluster,file = "l5_category.csv")

se_L18 <- SubsetSTData(se,expression = orig %in% c("L18"))

tablel18_umap <- data.frame(barcode = substr(colnames(se_L18),start = 1,stop = 18),UMAP_1 = se_L18@reductions$umap@cell.embeddings[,"UMAP_1"],UMAP_2= se_L18@reductions$umap@cell.embeddings[,"UMAP_2"])
tablel18_cluster <- data.frame(barcode = substr(colnames(se_L18),start = 1,stop = 18),category = se_L18@active.ident)

rownames(tablel18_umap) <- NULL
rownames(tablel18_cluster) <- NULL

write.csv(tablel18_umap,file = "l18_umap.csv")
write.csv(tablel18_cluster,file = "l18_category.csv")

#################Futher analysis combined the image and transcriptome data were inported to Loupe Brower Software.#########

