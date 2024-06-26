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
Idents(seu_ref) <- "celltype_cr"
DimPlot(seu_ref)
epi_levels <- levels(seu_ref$celltype_cr)

# label transfer
# seu <- NormalizeData(seu)
anchors <- FindTransferAnchors(reference = seu_ref, query = seu, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = seu_ref$celltype_cr, dims = 1:30)
seu$celltype_cr <- predictions$predicted.id
seu$celltype_cr <- factor(seu$celltype_cr, levels = epi_levels)
ImageDimPlot(seu, fov = "slide1", group.by = "celltype_cr", size = 1, axes = TRUE)
ImageDimPlot(seu, fov = "slide2", group.by = "celltype_cr", size = 1, axes = TRUE)

seu_epi <- seu


####
# label transfer, stromal ----
seu <- subset(seu_all, cells = str_names)

# load reference
seu_ref <- readRDS("data/reference_data/seu_060_stromal_npc20.RDS")
Idents(seu_ref) <- "celltype_cr"
DimPlot(seu_ref)
str_levels <- levels(seu_ref$celltype_cr)

# label transfer
# seu <- NormalizeData(seu)
anchors <- FindTransferAnchors(reference = seu_ref, query = seu, dims = 1:20, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = seu_ref$celltype_cr, dims = 1:20)
seu$celltype_cr <- predictions$predicted.id
seu$celltype_cr <- factor(seu$celltype_cr, levels = str_levels)
ImageDimPlot(seu, fov = "slide1", group.by = "celltype_cr", size = 1, axes = TRUE)
ImageDimPlot(seu, fov = "slide2", group.by = "celltype_cr", size = 1, axes = TRUE)

seu_str <- seu


###
# label transfer, Bcell ----
seu <- subset(seu_all, cells = Bcell_names)

# load reference
seu_ref <- readRDS("data/reference_data/seu_070_Bcell_npc30.RDS")
Idents(seu_ref) <- "celltype_cr"
DimPlot(seu_ref)
Bcell_levels <- levels(seu_ref$celltype_cr)

# label transfer
# seu <- NormalizeData(seu)
anchors <- FindTransferAnchors(reference = seu_ref, query = seu, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = seu_ref$celltype_cr, dims = 1:30)
seu$celltype_cr <- predictions$predicted.id
seu$celltype_cr <- factor(seu$celltype_cr, levels = Bcell_levels)
ImageDimPlot(seu, fov = "slide1", group.by = "celltype_cr", size = 1, axes = TRUE)
ImageDimPlot(seu, fov = "slide2", group.by = "celltype_cr", size = 1, axes = TRUE)

seu_Bcell <- seu


###
# label transfer, Tcell ----
seu <- subset(seu_all, cells = Tcell_names)

# load reference
seu_ref <- readRDS("data/reference_data/seu_080_Tcell_npc30_removed31_2.RDS")
Idents(seu_ref) <- "celltype_cr"
DimPlot(seu_ref)
Tcell_levels <- levels(seu_ref$celltype_cr)

# label transfer
# seu <- NormalizeData(seu)
anchors <- FindTransferAnchors(reference = seu_ref, query = seu, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = seu_ref$celltype_cr, dims = 1:30)
seu$celltype_cr <- predictions$predicted.id
seu$celltype_cr <- factor(seu$celltype_cr, levels = Tcell_levels)
ImageDimPlot(seu, fov = "slide1", group.by = "celltype_cr", size = 1, axes = TRUE)
ImageDimPlot(seu, fov = "slide2", group.by = "celltype_cr", size = 1, axes = TRUE)

seu_Tcell <- seu


###
# label transfer, myeloid ----
seu <- subset(seu_all, cells = mye_names)

# load reference
seu_ref <- readRDS("data/reference_data/seu_090_myeloid_npc20.RDS")
Idents(seu_ref) <- "celltype_cr"
DimPlot(seu_ref)
mye_levels <- levels(seu_ref$celltype_cr)

# label transfer
# seu <- NormalizeData(seu)
anchors <- FindTransferAnchors(reference = seu_ref, query = seu, dims = 1:20, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = seu_ref$celltype_cr, dims = 1:20)
seu$celltype_cr <- predictions$predicted.id
seu$celltype_cr <- factor(seu$celltype_cr, levels = mye_levels)
ImageDimPlot(seu, fov = "slide1", group.by = "celltype_cr", size = 1, axes = TRUE)
ImageDimPlot(seu, fov = "slide2", group.by = "celltype_cr", size = 1, axes = TRUE)

seu_mye <- seu


####
# merge and save the meta data ----
seu <- merge(x = seu_epi, y = list(seu_str, seu_Bcell, seu_Tcell, seu_mye))
seu <- JoinLayers(seu)
seu$celltype_cr <- factor(seu$celltype_cr, levels = c(epi_levels, str_levels, Bcell_levels, Tcell_levels, mye_levels))
seu$cellgroup <- factor(seu$cellgroup, levels = c("epithelial", "stromal", "Bcell", "Tcell", "myeloid"))
meta <- seu[[]]
save(meta, file = "RDSfiles/meta_celltype_cr.RData")


####
# re-create seurat object from the original seurat object ----

load(file = "RDSfiles/meta_celltype_cr.RData")
load(file = "RDSfiles/imgs_cts.RData")
seu <- CreateSeuratObject(counts = cts, meta.data = meta)
seu[["slide1"]] <- img1
seu[["slide2"]] <- img2
seu <- subset(seu, subset = qcCellsFlagged == FALSE)
ImageDimPlot(seu, fov = "slide1", group.by = "celltype_cr", size = 1, axes = TRUE, cols = "polychrome")
ImageDimPlot(seu, fov = "slide2", group.by = "celltype_cr", size = 1, axes = TRUE, cols = "polychrome")
saveRDS(seu, file = "RDSfiles/seu_031_labele_transfered_2.RDS")
save(Bcell_levels, epi_levels, mye_levels, str_levels, Tcell_levels, file = "RDSfiles/levels_celltype_cr.RData")


