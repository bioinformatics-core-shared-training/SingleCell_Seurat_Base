library(Seurat)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
library(tidyverse)
set.seed = 123 

seurat_object_500 <- readRDS("../RObjects/DI.500.rds")

seurat_object_500 <- FindNeighbors(seurat_object_500, 
                                reduction = "harmony", 
                                k.param = 17,
                                dims = 1:15)
seurat_object_500 <- FindClusters(seurat_object_500,,
                                resolution = 1,
                                algorithm = 4,
                                random.seed = 123,
                                cluster.name = "leiden_cluster")
saveRDS(seurat_object_500, file = "../RObjects/Clustered.500.rds")