library(Seurat)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
library(tidyverse)
set.seed = 123 
options(future.globals.maxSize = 3.0 * 1e9)

seurat_object_500 <- readRDS("../RObjects/Filtered.500.rds")
seurat_object_full <- readRDS("../RObjects/Filtered.full.rds")

##### For plot image
seurat_object_full_plot <- NormalizeData(seurat_object_full)
seurat_object_full_plot <- FindVariableFeatures(seurat_object_full_plot, 
                                      selection.method = "vst", 
                                      nfeatures = 3000)
VariableFeaturePlot(seurat_object_full_plot)

##### Normalisation #######

seurat_object_500 <- SCTransform(seurat_object_500, 
                                 vars.to.regress = "percent.mt", verbose = FALSE)
seurat_object_full <- SCTransform(seurat_object_full,
                                 vars.to.regress = "percent.mt", verbose = FALSE)

saveRDS(seurat_object_500, file = "../RObjects/SCT.500.rds")
saveRDS(seurat_object_full, file = "../RObjects/SCT.full.rds")

####### Dimensionality reduction ########

seurat_object_500 <- RunPCA(seurat_object_500, 
                        features = VariableFeatures(object = seurat_object_500))
seurat_object_full <- RunPCA(seurat_object_full,
                        features = VariableFeatures(object = seurat_object_full))

seurat_object_500 <- RunTSNE(seurat_object_500, 
                         reduction = "pca", 
                         dims = 1:15)
seurat_object_full <- RunTSNE(seurat_object_full,
                         reduction = "pca",
                         dims = 1:15)

seurat_object_500 <- RunUMAP(seurat_object_500, 
                         reduction = "pca", 
                         dims = 1:15)
seurat_object_full <- RunUMAP(seurat_object_full,
                         reduction = "pca",
                         dims = 1:15)

saveRDS(seurat_object_500, file = "../RObjects/DimRed.500.rds")
saveRDS(seurat_object_full, file = "../RObjects/DimRed.full.rds")

###### Intergration ######
