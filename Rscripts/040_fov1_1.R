####
# load packages ----
library(Seurat)
library(future)
library(ggplot2)
library(tidyverse)
options(future.globals.maxSize = 8000 * 1024^2)
library(Giotto)


####
# load data: re-created and filtered seurat object ----
seu <- readRDS("RDSfiles/seu_031_labele_transfered_2.RDS")
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
ImageDimPlot(seu, group.by = "celltype_cr", size = 1, cols = "polychrome", axes = TRUE)
crop <- Crop(seu[["slide1"]], y = c(8.5, 9.0), x = c(3.8, 4.3))    # x and y are flipped in the above plot
seu[["crop"]] <- crop
DefaultBoundary(seu[["crop"]]) <- "segmentation"
ImageDimPlot(seu, fov = "crop", group.by = "celltype_cr", cols = "polychrome", axes = TRUE)
# ImageDimPlot(seu, fov = "crop", molecules = c("OLFM4", "LCN2", "EPCAM", "C1QA"), mols.size = 0.5, axes = TRUE, nmols = 10000, cols = "polychrome")
# ImageFeaturePlot(seu, fov = "crop", features = "Max.PanCK", size = 2, axes = TRUE, max.cutoff = "q90")



# to plot selected types of cells with other types as "other"

# all <- levels(seu$celltype)
# sel <- c("TA", "Rib-C", "LAM-C", "PLGC-C", "Colono", "BEST4-C", "Sec-Pro", "Goblet", "Mature-G", "Paneth", "EE", "Tuft")
# rest <- all[! all %in% sel]
# dict <- setNames(rep("other", length(rest)), rest)
# seu$plot <- str_replace_all(seu$celltype, pattern = dict)
# seu$plot <- factor(seu$plot, levels = c("other", sel))

# alternatively
sel <- c("TA", "Colono", "Sec-Pro", "Goblet", "Paneth", "EE", "Tuft")
levels <- levels(seu$celltype_cr)
new_lev <- levels
new_lev[! new_lev %in% sel] <- "other"
seu$plot <- plyr::mapvalues(seu$celltype_cr, from = levels, to = new_lev)
seu$plot <- factor(seu$plot, levels = c("other", sel))
ImageDimPlot(seu, fov = "crop", group.by = "plot", cols = "polychrome", axes = TRUE)
ImageDimPlot(seu, fov = "crop", group.by = "cellgroup", cols = "polychrome", axes = TRUE)

