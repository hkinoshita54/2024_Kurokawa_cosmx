# install packages following 
# https://nanostring.com/wp-content/uploads/2023/01/LiverPublicDataRelease.html

if (!require("tiledb", quietly = TRUE)) remotes::install_github("TileDB-Inc/TileDB-R", force = TRUE, ref = "0.17.0")
if (!require("tiledbsc", quietly = TRUE)) remotes::install_github("tiledb-inc/tiledbsc", force = TRUE, ref = "8157b7d54398b1f957832f37fff0b173d355530e")
options(timeout = 600000000) ### set this to avoid timeout error
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

devtools::install_github("drieslab/Giotto@suite")
library(Giotto)
installGiottoEnvironment()
devtools::install_github("drieslab/GiottoData")

devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("jokergoo/circlize")
install.packages('NMF')
BiocManager::install("BiocNeighbors")
devtools::install_github("jinworks/CellChat")


# install umap-learn python package by using reticulate
# reticulate::install_python(version = '3.12.2')    # may not be necessaary
library(reticulate)
conda_remove("r-reticulate")
conda_create(envname = "r-reticulate", packages = "umap-learn", python_version = "3.12.2")

# run the following BEFORE loading packages
reticulate::use_condaenv("/Users/gut/Library/r-miniconda/envs/r-reticulate", required=T)
reticulate::py_module_available("umap")


# scCoustomize
# https://samuel-marsh.github.io/scCustomize/articles/Installation.html
install.packages("scCustomize")
BiocManager::install(c("ComplexHeatmap", "dittoSeq", "DropletUtils", "Nebulosa"))
install.packages(c("ggpubr", "hdf5r", "rliger"))

