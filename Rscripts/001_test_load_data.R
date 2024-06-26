####
# load packages ----
library(Seurat)
library(future)
library(ggplot2)
library(tidyverse)
options(future.globals.maxSize = 8000 * 1024^2)


####
# load data ----
seu <- readRDS("data_4/seuratObject.RDS")


####
# workaround for cell order error of the original seurat object from AtoMx ----
cts <- seu[["RNA"]]$counts
meta <- seu[[]]
img1 <- seu@images$Hayakawa.Slide1
img2 <- seu@images$Hayakawa.Slide2
seu <- CreateSeuratObject(counts = cts, meta.data = meta)    # keep only the raw counts and the meta data
seu[["slide1"]] <- img1
seu[["slide2"]] <- img2

####
# filter out low quality cells  ----
seu <- subset(seu, subset = qcCellsFlagged == FALSE)
saveRDS(seu, file = "RDSfiles/seu_002_filtered.RDS")

####
# subset  ----

# fovs of each biopsy sample (taken manually from the whole slide images)
fovs1_1 <- 1:23;fovs1_2 <- 24:45;fovs1_3 <- 46:59;fovs1_4 <- 60:74;fovs1_5 <- 75:94;fovs1_6 <- 95:105;
fovs2_1 <- 1:18;fovs2_2 <- 19:32;fovs2_3 <- 33:40;fovs2_4 <- 41:53;fovs2_5 <- 54:63;fovs2_6 <- 64:78;

# just use sample fovs2_1
seu <- subset(seu, subset = Run_Tissue_name == "Hayakawa-Slide2")
seu <- subset(seu, subset = fov %in% fovs2_1)
seucp <- seu


####
# clustering  ----
npcs = 30
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
seu <- RunPCA(seu, npcs = npcs)
seu <- FindNeighbors(seu, dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.6, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:npcs, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
FeaturePlot(seu, features = "Max.CD45", cols = c("lightgrey","darkred"), label = FALSE, repel = TRUE) + NoAxes()
FeaturePlot(seu, features = "PTPRC", cols = c("lightgrey","darkred"), label = FALSE, repel = TRUE) + NoAxes()
FeaturePlot(seu, features = "Max.PanCK", cols = c("lightgrey","darkred"), label = FALSE, repel = TRUE) + NoAxes()
FeaturePlot(seu, features = "EPCAM", cols = c("lightgrey","darkred"), label = FALSE, repel = TRUE) + NoAxes()
# sample quality poor? the cells do not cluster well... 


####
# Image based plots ----
ImageDimPlot(seu, molecules = "OLFM4")
ImageFeaturePlot(seu, features = "OLFM4", axes = TRUE) + coord_fixed()    # coord_fixed do nothing
crop <- Crop(seu[["slide2"]], y = c(9.2, 9.8), x = c(0.5, 1))    # x and y are flipped in the above plot
seu[["crop"]] <- crop
ImageFeaturePlot(seu, fov = "crop", features = "OLFM4", size = 2, axes = TRUE, max.cutoff = "q90")
ImageDimPlot(seu, fov = "crop", molecules = "OLFM4", size = 2, axes = TRUE)
DefaultBoundary(seu[["crop"]]) <- "segmentation"
Idents(seu) <- "orig.ident"
ImageDimPlot(seu, fov = "crop", molecules = c("OLFM4", "LGR5"), axes = TRUE, nmols = 10000, cols = "polychrome")
ImageFeaturePlot(seu, fov = "crop", features = "Max.CD45", size = 2, axes = TRUE, max.cutoff = "q90")
