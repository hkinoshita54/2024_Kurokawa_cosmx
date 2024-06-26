####
# load packages ----
library(Seurat)
library(future)
library(ggplot2)
library(tidyverse)
library(pheatmap)
options(future.globals.maxSize = 8000 * 1024^2)


####
# load data: re-created and filtered seurat object ----
seu <- readRDS("RDSfiles/seu_030_labele_transfered.RDS")
seu_all <- seu

####
# subset  ----

# fovs of each biopsy sample (taken manually from the whole slide images)
fovs1_1 <- 1:23;fovs1_2 <- 24:45;fovs1_3 <- 46:59;fovs1_4 <- 60:74;fovs1_5 <- 75:94;fovs1_6 <- 95:105;
fovs2_1 <- 1:18;fovs2_2 <- 19:32;fovs2_3 <- 33:40;fovs2_4 <- 41:53;fovs2_5 <- 54:63;fovs2_6 <- 64:78;

# just use sample fovs2_1
seu <- seu_all
seu <- subset(seu, subset = Run_Tissue_name == "Hayakawa-Slide1")
seu <- subset(seu, subset = fov %in% fovs1_1)    # sample "A1", pre-tx, 18F, D/C, MES2


####
# Image based plots ----
ImageDimPlot(seu, group.by = "celltype", size = 1, cols = "polychrome", axes = TRUE)
seu <- BuildNicheAssay(object = seu, fov = "slide1", group.by = "celltype", niches.k = 10, neighbors.k = 10)
ImageDimPlot(seu, group.by = "niches", size = 1, cols = "polychrome", axes = TRUE)
mat <- table(seu$celltype, seu$niches)
pheatmap(scale(mat))
