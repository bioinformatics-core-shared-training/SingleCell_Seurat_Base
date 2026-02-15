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

###### Integration ######

DIobj.500 <- readRDS("../RObjects/Filtered.500.rds")
DIobj.500[["RNA"]] <- split(DIobj.500[["RNA"]],
                                f = DIobj.500$SampleGroup)
DIobj.500 <- SCTransform(
  DIobj.500,
  assay = "RNA",
  vars.to.regress = "percent.mt",
  verbose = FALSE
)
DIobj.500 <- RunPCA(DIobj.500, 
                        features = VariableFeatures(object = DIobj.500))

DIobj.500 <- IntegrateLayers(object = DIobj.500, 
                                        method = HarmonyIntegration,
                                        orig.reduction = "pca", 
                                        new.reduction = "harmony",
                                        theta = 0.1)

DIobj.500 <- RunUMAP(DIobj.500,
                                reduction = "harmony",
                                dims = 1:15)
DIobj.500 <- JoinLayers(DIobj.500, assay = "RNA")

saveRDS(DIobj.500, file = "../RObjects/DI.500.rds")
