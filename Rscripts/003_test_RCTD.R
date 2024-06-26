####
# load packages ----
library(Seurat)
library(future)
library(ggplot2)
library(tidyverse)
library(spacexr)
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


####
# reference mapping by RCTD following Seurat vignette of mapping ----
# https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#human-lung-nanostring-cosmx-spatial-molecular-imager

query.counts <- GetAssayData(seu, assay = "RNA", layer = "counts")[, Cells(seu[["slide1"]])]
coords <- GetTissueCoordinates(seu[["slide1"]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# load reference data
seu_ref <- readRDS("data/reference_data/seu_100_combine_npc30.RDS")
Idents(seu_ref) <- "cellgroup"
DimPlot(seu_ref)
seu_ref <- JoinLayers(seu_ref)
counts <- GetAssayData(seu_ref, assay = "RNA", layer = "counts")
cluster <- as.factor(seu_ref$cellgroup)
names(cluster) <- colnames(seu_ref)
nUMI <- seu_ref$nCount_RNA
names(nUMI) <- colnames(seu_ref)
nUMI <- colSums(counts)
reference <- Reference(counts, cluster, nUMI)

# run RCTD with many cores
RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

annotations.df <- RCTD@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
seu$predicted.cellgroup <- annotations
keep.cells <- Cells(seu)[!is.na(seu$predicted.cellgroup)]
seu <- subset(seu, cells = keep.cells)

ImageDimPlot(seu, fov = "slide1", group.by = "predicted.cellgroup", axes = TRUE, size = 2)
# it does not look good because many cells are omitted

