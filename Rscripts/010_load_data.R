###
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
save(img1, img2, cts, file = "RDSfiles/imgs_cts.RData")
saveRDS(seu, file = "RDSfiles/seu_002_filtered.RDS")
