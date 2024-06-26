####
# load packages ----
library(reticulate)
use_condaenv("sc-env")
library(tidyverse)
library(Seurat)
library(future)
# library(scCustomize)
options(future.globals.maxSize = 8000 * 1024^2)
library(Matrix)


####
# load data ----
seu <- readRDS("RDSfiles/seu_031_labele_transfered_2.RDS")
seu_all <- seu
dir.create("adata")


# slide1 ----
seu <- subset(seu_all, subset = Run_Tissue_name == "Hayakawa-Slide1")

# as.anndata(x = seu, file_path = "./adata", file_name = "s1.h5ad") > did not work...

# manually prepare cts, meta, coords, etc. 
# https://smorabit.github.io/tutorials/8_velocyto/

# save metadata table
seu$barcode <- colnames(seu)
# seu$UMAP_1 <- seu@reductions$umap@cell.embeddings[,1]
# seu$UMAP_2 <- seu@reductions$umap@cell.embeddings[,2]
write.csv(seu@meta.data, file='adata/seu1_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
counts_matrix <- LayerData(seu, assay = 'RNA', layer = 'counts')
writeMM(counts_matrix, file = 'adata/seu1_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
# write.csv(seu@reductions$pca@cell.embeddings, file = 'adata/seu_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene' = rownames(counts_matrix)), file = 'adata/seu1_gene_names.csv',
  quote = F, row.names = F, col.names = F
)

# write spatial coordinates
coords <- GetTissueCoordinates(seu[["slide1"]][["centroids"]])
write.csv(coords, file='adata/seu1_coords.csv', quote=F, row.names=F)


# slide2 ----
seu <- subset(seu_all, subset = Run_Tissue_name == "Hayakawa-Slide2")

# as.anndata(x = seu, file_path = "./adata", file_name = "s1.h5ad") > did not work...

# manually prepare cts, meta, coords, etc. 
# https://smorabit.github.io/tutorials/8_velocyto/

# save metadata table
seu$barcode <- colnames(seu)
# seu$UMAP_1 <- seu@reductions$umap@cell.embeddings[,1]
# seu$UMAP_2 <- seu@reductions$umap@cell.embeddings[,2]
write.csv(seu@meta.data, file='adata/seu2_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
counts_matrix <- LayerData(seu, assay = 'RNA', layer = 'counts')
writeMM(counts_matrix, file = 'adata/seu2_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
# write.csv(seu@reductions$pca@cell.embeddings, file = 'adata/seu_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene' = rownames(counts_matrix)), file = 'adata/seu2_gene_names.csv',
  quote = F, row.names = F, col.names = F
)

# write spatial coordinates
coords <- GetTissueCoordinates(seu[["slide2"]][["centroids"]])
write.csv(coords, file='adata/seu2_coords.csv', quote=F, row.names=F)


# levels of celltype ----
load("RDSfiles/levels_celltype_cr.RData")
levels = c(epi_levels, str_levels, Bcell_levels, Tcell_levels, mye_levels)
write.table(
  data.frame('levels' = levels), file = 'adata/levels.csv',
  quote = F, row.names = F, col.names = F
)
