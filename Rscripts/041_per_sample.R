####
# load packages ----
library(Seurat)
library(future)
library(ggplot2)
library(tidyverse)
options(future.globals.maxSize = 8000 * 1024^2)
library(speckle)
# library(RColorBrewer)
# library(viridis)


####
# load data: re-created and filtered seurat object ----
seu <- readRDS("RDSfiles/seu_031_labele_transfered_2.RDS")
seu_all <- seu

####
# subset  ----

# fovs of each biopsy sample (taken manually from the whole slide images)
fovs1_1 <- 1:23;fovs1_2 <- 24:45;fovs1_3 <- 46:59;fovs1_4 <- 60:74;fovs1_5 <- 75:94;fovs1_6 <- 95:105;
fovs2_1 <- 1:18;fovs2_2 <- 19:32;fovs2_3 <- 33:40;fovs2_4 <- 41:53;fovs2_5 <- 54:63;fovs2_6 <- 64:78;
fov1_list <- list(fovs1_1, fovs1_2, fovs1_3, fovs1_4, fovs1_5, fovs1_6)
fov2_list <- list(fovs2_1, fovs2_2, fovs2_3, fovs2_4, fovs2_5, fovs2_6)

# put sample_id in slide1
seu1 <- subset(seu, subset = Run_Tissue_name == "Hayakawa-Slide1")
seu1$sample_id <- 0
for (i in 1:6){seu1$sample_id[seu1$fov %in% fov1_list[[i]]] <- i}
table(seu1$sample_id)

# put sample_id in slide2
seu2 <- subset(seu, subset = Run_Tissue_name == "Hayakawa-Slide2")
seu2$sample_id <- 0
for (i in 1:6){seu2$sample_id[seu2$fov %in% fov2_list[[i]]] <- 6+i}
table(seu2$sample_id)

# merge
seu <- merge(seu1, seu2)
load("RDSfiles/meta_celltype_cr.RData")
seu$cellgroup <- factor(seu$cellgroup, levels = c("epithelial", "stromal", "Bcell", "Tcell", "myeloid"))
seu$celltype_cr <- factor(seu$celltype_cr, levels = c(epi_levels, str_levels, Bcell_levels, Tcell_levels, mye_levels))
ImageDimPlot(seu, fov = "slide1", group.by = "sample_id", cols = "polychrome", axes = TRUE)
ImageDimPlot(seu, fov = "slide2", group.by = "sample_id", cols = "polychrome", axes = TRUE)

# add exp_group data
seu$time <- "pre"
seu$time[seu$sample_id %in% c(4, 5, 6, 10, 11, 12)] <- "post"
seu$time <- factor(seu$time, levels = c("pre", "post"))

seu$drug <- "tol"
seu$drug[seu$sample_id %in% c(7, 8, 9, 10, 11, 12)] <- "int"
seu$drug <- factor(seu$drug, levels = c("tol", "int"))

seu$exp_group <- paste(seu$drug, seu$time, sep = "_")
seu$exp_group <- factor(seu$exp_group, levels = c("tol_pre", "tol_post", "int_pre", "int_post"))


seu <- JoinLayers(seu)
saveRDS(seu, file = "RDSfiles/seu_031_labele_transfered_2.RDS")


