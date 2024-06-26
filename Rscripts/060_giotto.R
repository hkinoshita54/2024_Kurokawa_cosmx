####
# load packages ----
library(Seurat)
library(future)
library(ggplot2)
library(tidyverse)
options(future.globals.maxSize = 8000 * 1024^2)
library(Giotto)



####
# Giotto ----

pal10 = c("#66C5CC","#F6CF71","#F89C74","#DCB0F2","#87C55F","#9EB9F3","#FE88B1","#C9DB74","#8BE0A4","#B3B3B3")
viv10 = c("#E58606","#5D69B1","#52BCA3","#99C945","#CC61B0","#24796C","#DAA51B","#2F8AC4","#764E9F","#A5AA99")
results_folder = 'results_giotto'
my_python_path = NULL
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  return_plot = FALSE,
                                  python_path = my_python_path)


## provide path to nanostring folder
data_path = '/Users/gut/WORKSPACE/2024_Kurokawa_cosmx/data_4/flatFiles/HayakawaSlide1/'

## create giotto cosmx object > still not working...
fov_join = createGiottoCosMxObject(cosmx_dir = data_path,
                                   data_to_use = 'subcellular', # only subcellular
                                   FOVs = 1:105,
                                   instructions = instrs)

showGiottoFeatInfo(fov_join)
showGiottoSpatialInfo(fov_join)






tx_coord_all = data.table::fread(paste0(data_path, 'HayakawaSlide1_tx_file.csv.gz'))


# correct "FOV" in fov filse to "fov" to be read in squidpy ----
fov = data.table::fread(paste0(data_path, 'HayakawaSlide1_fov_positions_file.csv.gz'))
colnames(fov)[colnames(fov) == "FOV"] <- "fov"
write_csv(fov, file = paste0(data_path, 'HayakawaSlide1_fov2_positions_file.csv.gz'))
fov = data.table::fread(paste0(data_path, 'HayakawaSlide1_fov2_positions_file.csv.gz'))
