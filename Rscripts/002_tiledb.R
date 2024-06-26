####
# load packages ----

library(tiledb)
library(tiledbsc)
library(ggplot2)
library(Seurat)
library(RColorBrewer)



tiledbURI <- "data/82bb400f-e4e4-49bc-b804-44b880b11b06_TileDB/"
tiledb_scdataset <- tiledbsc::SOMACollection$new(uri = tiledbURI, verbose = FALSE)
names(tiledb_scdataset$somas)
names(tiledb_scdataset$somas$RNA$members)

# raw counts
counts <- tiledb_scdataset$somas$RNA$X$members$counts$to_matrix(batch_mode = TRUE)
dim(counts)
counts[1:4,1:4]

#normalized counts (first sample)
norm <- tiledb_scdataset$somas$RNA_normalized_461882b3.3c71.46c6.88b6.3df74f27ec01_1$X$members$data$to_matrix(batch_mode = TRUE)
dim(norm)
norm[1:4,1:4]

# metadata
metadata <- tiledb_scdataset$somas$RNA$obs$to_dataframe()
dim(metadata)
metadata[1:4,1:10]

cellCoords <- tiledb_scdataset$somas$RNA$obs$to_dataframe(
  attrs = c("x_FOV_px", "y_FOV_px", "x_slide_mm", "y_slide_mm", "slide_ID_numeric", "Run_Tissue_name", "fov"))
head(cellCoords)

ggplot(cellCoords, aes(x=x_slide_mm, y=y_slide_mm))+
  geom_point(alpha = 0.05, size = 0.01)+
  facet_wrap(~Run_Tissue_name)+
  coord_equal()+
  labs(title = "Cell coordinates in XY space")

transcriptCoords <- tiledb::tiledb_array(
  tiledb_scdataset$somas$RNA$obsm$members$transcriptCoords$uri,
  return_as="data.frame")[]
head(transcriptCoords)
tail(transcriptCoords)


slide <- 2
fov <- 35

slideName <- unique(cellCoords$Run_Tissue_name[cellCoords$slide_ID_numeric == slide])

fovCoords <- cellCoords[cellCoords$slide_ID_numeric == slide & cellCoords$fov == fov,]
fovTranscriptCoords <- transcriptCoords[transcriptCoords$slideID == slide & transcriptCoords$fov == fov,]
targetCounts <- table(fovTranscriptCoords$target)
targets <- names(targetCounts[which(targetCounts >= 2500 & targetCounts <= 3000)])
fovTranscriptCoords <- fovTranscriptCoords[fovTranscriptCoords$target %in% targets,]

ggplot(fovCoords, aes(x=x_FOV_px, y=y_FOV_px))+
  geom_point(alpha = 0.6, size = 0.1, color = "black")+
  geom_point(data = fovTranscriptCoords, mapping = aes(x=x_FOV_px, y=y_FOV_px, color = target), size = 0.3, alpha = 0.2)+
  theme_bw()+
  coord_equal()+
  guides(colour = guide_legend(override.aes = list(size=1, alpha=1)))+
  labs(color = "RNA Target", title = paste0("RNA Transcripts in\n", slideName, "\nFOV", fov))

# analysis results
pca <- tiledb_scdataset$somas$RNA$obsm$members$dimreduction_approximatepca_50a7107d.448a.4936.9f36.583d48cdf064_1$to_matrix()
pca[1:4,1:4]

neighbs <- tiledb_scdataset$somas$RNA$obsp$members$graph_nn_c0739609.c124.4f54.8e60.46cd777bfad2_1_nn$to_seurat_graph()
neighbs[1:4,1:25]

# to seurat
RNA_seurat <- tiledb_scdataset$to_seurat(somas = c("RNA"), batch_mode = TRUE)
RNA_seurat
RNA_seurat@meta.data <- metadata




####
# below did not work w/o predetermined annotation ----


# colorColumn = "cellType"
umap_TileDB <- as.data.frame(tiledb_scdataset$somas$RNA$obsm$members$dimreduction_approximateumap_d487376e.4e36.4e20.b28b.3500f97d433f_1$to_matrix())
# umap_TileDB$colorBy <- tiledb_scdataset$somas$RNA$obs$to_dataframe(attrs = colorColumn)[[1]]

colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
colorSchemes <- c("PuOr", "Dark2", "Set2", "BrBG")
colrs <- colrs[colorSchemes,]
col_vec <- unique(unlist(mapply(brewer.pal, colrs$maxcolors, colorSchemes)))

col_vec <- col_vec[-grep(col_vec, pattern = "#F|#E|#D")]

plotting <- function(data, title, Xcol, Ycol, Xname, Yname, color, 
                     size = 0.02, alpha = 0.05){
  gp <- ggplot(data, aes(x = .data[[Xcol]], y = .data[[Ycol]], 
                         color = .data[[color]]))+
    geom_point(size = size, alpha = alpha)+
    coord_equal()+
    labs(title = title, color = colorColumn,
         x = Xname, y = Yname)+
    scale_color_manual(values = col_vec)+
    guides(colour = guide_legend(override.aes = list(size=1,
                                                     alpha=1)))
  
  return(gp)
}

umapPlot <- function(data, title, colorBy = "colorBy", ...){
  return(plotting(data, title, Xcol = "APPROXIMATEUMAP_1", 
                  Ycol = "APPROXIMATEUMAP_2", Xname = "UMAP1", 
                  Yname = "UMAP2", color = colorBy, ...))
}

pcaPlot <- function(data, title, colorBy = "colorBy", ...){
  return(plotting(data, title, Xcol = "PCA_1", 
                  Ycol = "PCA_2", Xname = "PCA_1", 
                  Yname = "PCA_2", color = colorBy, ...))
}

xySlidePlot <- function(data, title, colorBy = "colorBy", ...){
  return(plotting(data, title, Xcol = "x_slide_mm", 
                  Ycol = "y_slide_mm", Xname = "x_slide_mm", 
                  Yname = "y_slide_mm", color = colorBy, ...))
}




# umapPlot(umap_TileDB, "TileDB - R")

