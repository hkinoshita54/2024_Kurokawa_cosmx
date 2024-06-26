####
# load packages ----
library(Seurat)
library(future)
library(ggplot2)
library(tidyverse)
options(future.globals.maxSize = 8000 * 1024^2)
library(speckle)


####
# set conditions ----
description = "050_DA_speckle"
path = paste0("plot/", description)
dir.create(path)



####
# load data: re-created and filtered seurat object ----
seu <- readRDS("RDSfiles/seu_031_labele_transfered_2.RDS")
seu_all <- seu
load("RDSfiles/levels_celltype_cr.RData")



####
# differential abundance testing by speckle ----
propres <- propeller(clusters = seu$cellgroup, sample = seu$sample_id, group = seu$exp_group)
write.table(propres, file = "results/speckle/cellgroup_exp_group.txt", sep ="\t", col.names = T,row.names = F)
props <- getTransformedProps(clusters = seu$cellgroup, sample = seu$sample_id)
transprop <- as.data.frame(props$TransformedProps) %>% pivot_wider(names_from = sample, values_from = Freq) %>% t()
write.table(transprop, file = "results/speckle/cellgroup.tsv", sep = "\t")
ggplot(as.data.frame(props$Proportions), aes(x=sample, y=Freq, fill=clusters)) + 
  geom_bar(position="stack", stat="identity") +
  # scale_fill_brewer(palette = "Dark2") +
  theme_classic()
ggsave("cellgroup_prop.png", path = path, width = 4, height = 4, units = "in", dpi = 150)


####
# epithelial ----
seu <- subset(seu_all, subset = cellgroup == "epithelial")
seu$celltype_cr <- factor(seu$celltype_cr, levels = epi_levels)
propres <- propeller(clusters = seu$celltype_cr, sample = seu$sample_id, group = seu$exp_group)
write.table(propres, file = "results/speckle/epi_exp_group.txt", sep ="\t", col.names = T,row.names = F)
props <- getTransformedProps(clusters = seu$celltype_cr, sample = seu$sample_id)
transprop <- as.data.frame(props$TransformedProps) %>% pivot_wider(names_from = sample, values_from = Freq) %>% t()
write.table(transprop, file = "results/speckle/transprop_epi.tsv", sep = "\t")
ggplot(as.data.frame(props$Proportions), aes(x=sample, y=Freq, fill=clusters)) + 
  geom_bar(position="stack", stat="identity") +
  # scale_fill_brewer(palette = "Dark2") +
  theme_classic()
ggsave("epi_prop.png", path = path, width = 4, height = 4, units = "in", dpi = 150)


####
# stromal ----
seu <- subset(seu_all, subset = cellgroup == "stromal")
seu$celltype_cr <- factor(seu$celltype_cr, levels = str_levels)
propres <- propeller(clusters = seu$celltype_cr, sample = seu$sample_id, group = seu$exp_group)
write.table(propres, file = "results/speckle/str_exp_group.txt", sep ="\t", col.names = T,row.names = F)
props <- getTransformedProps(clusters = seu$celltype_cr, sample = seu$sample_id)
transprop <- as.data.frame(props$TransformedProps) %>% pivot_wider(names_from = sample, values_from = Freq) %>% t()
write.table(transprop, file = "results/speckle/transprop_str.tsv", sep = "\t")
ggplot(as.data.frame(props$Proportions), aes(x=sample, y=Freq, fill=clusters)) + 
  geom_bar(position="stack", stat="identity") +
  # scale_fill_brewer(palette = "Dark2") +
  theme_classic()
ggsave("str_prop.png", path = path, width = 4, height = 4, units = "in", dpi = 150)


####
# Bcell ----
seu <- subset(seu_all, subset = cellgroup == "Bcell")
seu$celltype_cr <- factor(seu$celltype_cr, levels = Bcell_levels)
propres <- propeller(clusters = seu$celltype_cr, sample = seu$sample_id, group = seu$exp_group)
write.table(propres, file = "results/speckle/Bcell_exp_group.txt", sep ="\t", col.names = T,row.names = F)
props <- getTransformedProps(clusters = seu$celltype_cr, sample = seu$sample_id)
transprop <- as.data.frame(props$TransformedProps) %>% pivot_wider(names_from = sample, values_from = Freq) %>% t()
write.table(transprop, file = "results/speckle/transprop_Bcell.tsv", sep = "\t")
ggplot(as.data.frame(props$Proportions), aes(x=sample, y=Freq, fill=clusters)) + 
  geom_bar(position="stack", stat="identity") +
  # scale_fill_brewer(palette = "Dark2") +
  theme_classic()
ggsave("Bcell_prop.png", path = path, width = 4, height = 4, units = "in", dpi = 150)


####
# Tcell ----
seu <- subset(seu_all, subset = cellgroup == "Tcell")
seu$celltype_cr <- factor(seu$celltype_cr, levels = Tcell_levels)
propres <- propeller(clusters = seu$celltype_cr, sample = seu$sample_id, group = seu$exp_group)
write.table(propres, file = "results/speckle/Tcell_exp_group.txt", sep ="\t", col.names = T,row.names = F)
props <- getTransformedProps(clusters = seu$celltype_cr, sample = seu$sample_id)
transprop <- as.data.frame(props$TransformedProps) %>% pivot_wider(names_from = sample, values_from = Freq) %>% t()
write.table(transprop, file = "results/speckle/transprop_Tcell.tsv", sep = "\t")
ggplot(as.data.frame(props$Proportions), aes(x=sample, y=Freq, fill=clusters)) + 
  geom_bar(position="stack", stat="identity") +
  # scale_fill_brewer(palette = "Dark2") +
  theme_classic()
ggsave("Tcell_prop.png", path = path, width = 4, height = 4, units = "in", dpi = 150)


####
# myeloid ----
seu <- subset(seu_all, subset = cellgroup == "myeloid")
seu$celltype_cr <- factor(seu$celltype_cr, levels = mye_levels)
propres <- propeller(clusters = seu$celltype_cr, sample = seu$sample_id, group = seu$exp_group)
write.table(propres, file = "results/speckle/mye_exp_group.txt", sep ="\t", col.names = T,row.names = F)
props <- getTransformedProps(clusters = seu$celltype_cr, sample = seu$sample_id)
transprop <- as.data.frame(props$TransformedProps) %>% pivot_wider(names_from = sample, values_from = Freq) %>% t()
write.table(transprop, file = "results/speckle/transprop_mye.tsv", sep = "\t")
ggplot(as.data.frame(props$Proportions), aes(x=sample, y=Freq, fill=clusters)) + 
  geom_bar(position="stack", stat="identity") +
  # scale_fill_brewer(palette = "Dark2") +
  theme_classic()
ggsave("mye_prop.png", path = path, width = 4, height = 4, units = "in", dpi = 150)