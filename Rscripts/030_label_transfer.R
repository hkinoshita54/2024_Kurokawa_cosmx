####
# load packages ----
library(Seurat)
library(future)
library(ggplot2)
library(tidyverse)
options(future.globals.maxSize = 8000 * 1024^2)


####
# load data: re-created and filtered seurat object ----
seu <- readRDS("RDSfiles/seu_020_cellgroup.RDS")
seu_all <- seu
load("RDSfiles/cellnames.RData")


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
ImageDimPlot(seu, fov = "slide1", group.by = "celltype", size = 1, axes = TRUE)
ImageDimPlot(seu, fov = "slide2", group.by = "celltype", size = 1, axes = TRUE)

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
ImageDimPlot(seu, fov = "slide1", group.by = "celltype", size = 1, axes = TRUE)
ImageDimPlot(seu, fov = "slide2", group.by = "celltype", size = 1, axes = TRUE)

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
ImageDimPlot(seu, fov = "slide1", group.by = "celltype", size = 1, axes = TRUE)
ImageDimPlot(seu, fov = "slide2", group.by = "celltype", size = 1, axes = TRUE)

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
ImageDimPlot(seu, fov = "slide1", group.by = "celltype", size = 1, axes = TRUE)
ImageDimPlot(seu, fov = "slide2", group.by = "celltype", size = 1, axes = TRUE)

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
ImageDimPlot(seu, fov = "slide1", group.by = "celltype", size = 1, axes = TRUE)
ImageDimPlot(seu, fov = "slide2", group.by = "celltype", size = 1, axes = TRUE)

seu_mye <- seu


####
# merge and save the metadata (celltype) ----
seu <- merge(x = seu_epi, y = list(seu_str, seu_Bcell, seu_Tcell, seu_mye))

# the following did not work (when some cells are removed, the original image can not be added)
# seu[[c("slide1", "slide2", "slide1.2", "slide2.2", "slide1.3", "slide2.3", "slide1.4", "slide2.4", "slide1.5", "slide2.5")]] <- NULL
# seu <- JoinLayers(seu)
# seu$celltype <- factor(seu$celltype, levels = c(epi_levels, str_levels, Bcell_levels, Tcell_levels, mye_levels))
# seu$cellgroup <- factor(seu$cellgroup, levels = c("epithelial", "stromal", "Bcell", "Tcell", "myeloid"))
# load(file = "RDSfiles/imgs.RData")
# seu[["slide1"]] <- img1
# seu[["slide2"]] <- img2
# ImageDimPlot(seu, fov = "slide1", group.by = "celltype", size = 1, axes = TRUE, cols = "polychrome")
# ImageDimPlot(seu, fov = "slide2", group.by = "celltype", size = 1, axes = TRUE, cols = "polychrome")
# saveRDS(seu, file = "RDSfiles/seu_030_labele_transfered.RDS")


####
# reduce number of cell types following the original paper ----
# dict = c(
#   "LAM-C" = "Colono", "PLGC-C" = "Colono", "Mature-G" = "Goblet",
#   "Rib-F" = "IER-F", "MT-F" = "S1-F", "Act.-EC" = "EC", "LEC" = "EC",
#   "IgGL-P" = "IgG-P", "IgAL-P" = "IgA-P", 
#   "CD8-CTL" = "CD8", "CD8-TRM" = "CD8", "CD8-FGFBP2" = "CD8", 
#   "CD4-naive" = "CD4", "CD4-ANXA1" = "CD4",
#   "DN-TNF" = "DN", "DN-EOMES" = "DN",
#   "Rib-M0" = "M0", "DC-CD1C" = "DC", "DC-CCL22" = "DC"
# )
# seu$celltype_2 <- str_replace_all(seu$celltype, pattern = dict)
# levels <- levels(seu$celltype)
# levels_2 <- str_replace_all(levels, pattern = dict) %>% unique
# seu$celltype_2 <- factor(seu$celltype_2, levels = levels_2)
# ImageDimPlot(seu, fov = "slide1", group.by = "celltype_2", size = 1, axes = TRUE, cols = "polychrome")
# saveRDS(seu, file = "RDSfiles/seu_030_labele_transfered.RDS")

####
# reduce number of cell types following the original paper ----
# dict_2 = c(
#   "Rib-C" = "Colono", 
#   "S2a-F" = "S2-F", "S2b-F" = "S2-F", "IER-F" = "S1-F",
#   "IgA-HSP-P" = "IgA-P", "IER-P" = "IgA-P", "IGLL5-P" = "IgA-P",
#   "Naive-B" = "Bcell", "GC-B" = "Bcell", "Mem.-B" = "Bcell",
#   "CD4-CCL20" = "CD4", "Rib-T" = "Tcell-NOS", "MT-T" = "Tcell-NOS",
#   "IDA-Mac" = "M1", "Infl.-Mono" = "M1"
# )
# seu$celltype_3 <- str_replace_all(seu$celltype_2, pattern = dict_2)
# levels_3 <- str_replace_all(levels_2, pattern = dict_2) %>% unique
# seu$celltype_3 <- factor(seu$celltype_3, levels = levels_3)
# ImageDimPlot(seu, fov = "slide1", group.by = "celltype_3", size = 1, axes = TRUE, cols = "polychrome")
# saveRDS(seu, file = "RDSfiles/seu_030_labele_transfered.RDS")

# correct level "EC" to "BEC"
# levels <- levels(seu$celltype)
# levels[levels == "EC"] <- "BEC"
# levels(seu$celltype) <- levels
# saveRDS(seu, file = "RDSfiles/seu_030_labele_transfered.RDS")
