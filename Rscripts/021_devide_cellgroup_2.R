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
# divide the data into 5 cell groups ----
# use the data by Garrido-Trigo et al. (nat commun 2023)
# label transfer following Seurat vignette of mapping
# https://satijalab.org/seurat/articles/integration_mapping

# load reference
seu_ref <- readRDS("data/reference_data/seu_101_combine_npc30.RDS")
Idents(seu_ref) <- "cellgroup"
DimPlot(seu_ref)

# label transfer
seu <- NormalizeData(seu)
anchors <- FindTransferAnchors(reference = seu_ref, query = seu, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = seu_ref$cellgroup, dims = 1:30)
seu$cellgroup <- predictions$predicted.id
seu$cellgroup <- factor(seu$cellgroup, levels = c("epithelial", "stromal", "Bcell", "Tcell", "myeloid"))
ImageDimPlot(seu, fov = "slide1", group.by = "cellgroup", size = 1, axes = TRUE)
ImageDimPlot(seu, fov = "slide2", group.by = "cellgroup", size = 1, axes = TRUE)

# get cell names in each group
epi_names <- colnames(seu[, seu$cellgroup == "epithelial"])
str_names <- colnames(seu[, seu$cellgroup == "stromal"])
Bcell_names <- colnames(seu[, seu$cellgroup == "Bcell"])
Tcell_names <- colnames(seu[, seu$cellgroup == "Tcell"])
mye_names <- colnames(seu[, seu$cellgroup == "myeloid"])

save(epi_names, str_names, Bcell_names, Tcell_names, mye_names, file = "RDSfiles/cellnames_2.RData")
saveRDS(seu, file = "RDSfiles/seu_021_cellgroup_2.RDS")
