####
# load packages ----
library(Seurat)
library(future)
library(ggplot2)
library(tidyverse)
options(future.globals.maxSize = 8000 * 1024^2)


####
# load data: re-created and filtered seurat object ----
seu <- readRDS("RDSfiles/seu_021_cellgroup_2.RDS")
seu_all <- seu
load("RDSfiles/cellnames_2.RData")


####
# label transfer, epithelial ----
seu <- subset(seu_all, cells = epi_names)

# load reference
seu_ref <- readRDS("data/reference_data/seu_050_epithelial_npc30.RDS")
Idents(seu_ref) <- "celltype"
DimPlot(seu_ref)
epi_levels <- levels(seu_ref$celltype)

# label transfer
# seu <- NormalizeData(seu)
anchors <- FindTransferAnchors(reference = seu_ref, query = seu, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = seu_ref$celltype, dims = 1:30)
seu$celltype <- predictions$predicted.id
seu$celltype <- factor(seu$celltype, levels = epi_levels)
# ImageDimPlot(seu, fov = "slide1", group.by = "celltype", size = 1, axes = TRUE)
# ImageDimPlot(seu, fov = "slide2", group.by = "celltype", size = 1, axes = TRUE)

seu_epi <- seu


####
# label transfer, stromal ----
seu <- subset(seu_all, cells = str_names)

# load reference
seu_ref <- readRDS("data/reference_data/seu_060_stromal_npc20.RDS")
Idents(seu_ref) <- "celltype"
DimPlot(seu_ref)
str_levels <- levels(seu_ref$celltype)

# label transfer
# seu <- NormalizeData(seu)
anchors <- FindTransferAnchors(reference = seu_ref, query = seu, dims = 1:20, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = seu_ref$celltype, dims = 1:20)
seu$celltype <- predictions$predicted.id
seu$celltype <- factor(seu$celltype, levels = str_levels)
# ImageDimPlot(seu, fov = "slide1", group.by = "celltype", size = 1, axes = TRUE)
# ImageDimPlot(seu, fov = "slide2", group.by = "celltype", size = 1, axes = TRUE)

seu_str <- seu


###
# label transfer, Bcell ----
seu <- subset(seu_all, cells = Bcell_names)

# load reference
seu_ref <- readRDS("data/reference_data/seu_070_Bcell_npc30.RDS")
Idents(seu_ref) <- "celltype"
DimPlot(seu_ref)
Bcell_levels <- levels(seu_ref$celltype)

# label transfer
# seu <- NormalizeData(seu)
anchors <- FindTransferAnchors(reference = seu_ref, query = seu, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = seu_ref$celltype, dims = 1:30)
seu$celltype <- predictions$predicted.id
seu$celltype <- factor(seu$celltype, levels = Bcell_levels)
# ImageDimPlot(seu, fov = "slide1", group.by = "celltype", size = 1, axes = TRUE)
# ImageDimPlot(seu, fov = "slide2", group.by = "celltype", size = 1, axes = TRUE)

seu_Bcell <- seu


###
# label transfer, Tcell ----
seu <- subset(seu_all, cells = Tcell_names)

# load reference
seu_ref <- readRDS("data/reference_data/seu_080_Tcell_npc30_removed31_2.RDS")
Idents(seu_ref) <- "celltype"
DimPlot(seu_ref)
Tcell_levels <- levels(seu_ref$celltype)

# label transfer
# seu <- NormalizeData(seu)
anchors <- FindTransferAnchors(reference = seu_ref, query = seu, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = seu_ref$celltype, dims = 1:30)
seu$celltype <- predictions$predicted.id
seu$celltype <- factor(seu$celltype, levels = Tcell_levels)
# ImageDimPlot(seu, fov = "slide1", group.by = "celltype", size = 1, axes = TRUE)
# ImageDimPlot(seu, fov = "slide2", group.by = "celltype", size = 1, axes = TRUE)

seu_Tcell <- seu


###
# label transfer, myeloid ----
seu <- subset(seu_all, cells = mye_names)

# load reference
seu_ref <- readRDS("data/reference_data/seu_090_myeloid_npc20.RDS")
Idents(seu_ref) <- "celltype"
DimPlot(seu_ref)
mye_levels <- levels(seu_ref$celltype)

# label transfer
# seu <- NormalizeData(seu)
anchors <- FindTransferAnchors(reference = seu_ref, query = seu, dims = 1:20, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = seu_ref$celltype, dims = 1:20)
seu$celltype <- predictions$predicted.id
seu$celltype <- factor(seu$celltype, levels = mye_levels)
# ImageDimPlot(seu, fov = "slide1", group.by = "celltype", size = 1, axes = TRUE)
# ImageDimPlot(seu, fov = "slide2", group.by = "celltype", size = 1, axes = TRUE)

seu_mye <- seu


####
# merge and save the metadata (celltype) ----
seu <- merge(x = seu_epi, y = list(seu_str, seu_Bcell, seu_Tcell, seu_mye))
celltype <- select(seu[[]], cell_id, celltype)
celltype$celltype <- factor(celltype$celltype, levels = c(epi_levels, str_levels, Bcell_levels, Tcell_levels, mye_levels))
saveRDS(celltype, file = "RDSfiles/celltype.RDS")

# add to the previous seu obj
seu <- readRDS("RDSfiles/seu_031_labele_transfered_2.RDS")
seu[[]] <- inner_join(seu[[]], celltype, by = "cell_id")
saveRDS(seu, file = "RDSfiles/seu_031_labele_transfered_2.RDS")

