####
# load packages ----
library(Seurat)
library(future)
library(ggplot2)
library(tidyverse)
options(future.globals.maxSize = 8000 * 1024^2)


####
# load data: re-created and filtered seurat object ----
seu <- readRDS("RDSfiles/seu_002_filtered.RDS")


####
# subset  ----

# fovs of each biopsy sample (taken manually from the whole slide images)
fovs1_1 <- 1:23;fovs1_2 <- 24:45;fovs1_3 <- 46:59;fovs1_4 <- 60:74;fovs1_5 <- 75:94;fovs1_6 <- 95:105;
fovs2_1 <- 1:18;fovs2_2 <- 19:32;fovs2_3 <- 33:40;fovs2_4 <- 41:53;fovs2_5 <- 54:63;fovs2_6 <- 64:78;

# just use sample fovs2_1
seu <- subset(seu, subset = Run_Tissue_name == "Hayakawa-Slide1")
seu <- subset(seu, subset = fov %in% fovs1_1)    # sample "A1", pre-tx, 18F, D/C, MES2
seucp <- seu


####
# reference mapping following Seurat vignette of mapping ----
# https://satijalab.org/seurat/articles/integration_mapping

# label transfer
seu_ref <- readRDS("data/reference_data/seu_100_combine_npc30.RDS")
Idents(seu_ref) <- "cellgroup"
DimPlot(seu_ref)
seu <- seucp
seu <- NormalizeData(seu)
anchors <- FindTransferAnchors(reference = seu_ref, query = seu, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = seu_ref$cellgroup, dims = 1:30)
seu <- AddMetaData(seu, metadata = predictions)
table(seu$predicted.id)
VlnPlot(seu, c("MS4A1", "EPCAM", "C1QC", "VIM", "CD3D"), group.by = "predicted.id")
ImageDimPlot(seu, group.by = "predicted.id", size = 2, axes = TRUE)

# # UMAP projection
# 
# # to avoid error, first re-run RunUMAP() with "return.model=TRUE"
# # https://github.com/satijalab/seurat/issues/3615
# seu_ref <- RunUMAP(seu_ref, reduction = "harmony", dims = 1:30, verbose = FALSE, return.model = TRUE)
# 
# # do UMAP projection
# seu <- TransferData(anchorset = anchors, reference = seu_ref, query = seu, refdata = list(cellgroup = "cellgroup"))
# seu <- IntegrateEmbeddings(anchorset = anchors, reference = seu_ref, query = seu, new.reduction.name = "ref.pca")
# seu <- ProjectUMAP(query = seu, query.reduction = "ref.pca", reference = seu_ref, reference.reduction = "pca", reduction.model = "umap")
# 
# # the below is the alternative wrapper function for the above three lines
# seu <- MapQuery(anchorset = anchors, reference = seu_ref, query = seu, 
#                 refdata = list(cellgroup = "cellgroup"), reference.reduction = "pca", reduction.model = "umap")
# 
# DimPlot(seu, reduction = "ref.umap", group.by = "predicted.cellgroup")


####
# Image based plots ----
crop <- Crop(seu[["slide1"]], y = c(8.5, 9.0), x = c(3.8, 4.3))    # x and y are flipped in the above plot
seu[["crop"]] <- crop
DefaultBoundary(seu[["crop"]]) <- "segmentation"
ImageDimPlot(seu, fov = "crop", group.by = "predicted.id", cols = "polychrome", axes = TRUE)
ImageDimPlot(seu, fov = "crop", molecules = c("OLFM4", "LCN2", "EPCAM", "C1QA"), mols.size = 0.5, axes = TRUE, nmols = 10000, cols = "polychrome")
ImageFeaturePlot(seu, fov = "crop", features = "OLFM4", size = 2, axes = TRUE, max.cutoff = "q90")


####
# niche analysis ----
seu <- BuildNicheAssay(object = seu, fov = "slide1", group.by = "predicted.id", niches.k = 3, neighbors.k = 30)
celltype.plot <- ImageDimPlot(seu, group.by = "predicted.id", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(seu, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches") + 
  scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D"))
celltype.plot | niche.plot

table(seu$predicted.id, seu$niches)
