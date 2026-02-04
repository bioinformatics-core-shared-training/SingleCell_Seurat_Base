library(Seurat)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
library(tidyverse)
set.seed = 123 
options(future.globals.maxSize = 3.0 * 1e9)

seurat_object_500 <- readRDS("../RObjects/Filtered.500.rds")
seurat_object_full <- readRDS("../RObjects/Filtered.full.rds")

seurat_object_500 <- SCTransform(seurat_object_500, 
                                 vars.to.regress = "percent.mt", verbose = FALSE)
seurat_object_full <- SCTransform(seurat_object_full,
                                 vars.to.regress = "percent.mt", verbose = FALSE)

saveRDS(seurat_object_500, file = "../RObjects/SCT.500.rds")
saveRDS(seurat_object_full, file = "../RObjects/SCT.full.rds")
