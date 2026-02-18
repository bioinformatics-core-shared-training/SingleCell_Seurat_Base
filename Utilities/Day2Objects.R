library(Seurat)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
library(tidyverse)
set.seed = 123 
options(future.globals.maxSize = 3.0 * 1e9)

seurat_object_500 <- readRDS("../RObjects/Filtered.500.rds")
seurat_object_full <- readRDS("RObjects/Filtered.full.rds")

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
DIobj.full <- readRDS("RObjects/Filtered.full.rds")

temp_metadata <- DIobj.full@meta.data %>%
  rownames_to_column("Cell") %>%
  mutate(SampleGroup = str_remove(orig.ident, "-.*")) %>% 
  mutate(SampleName = orig.ident) %>%
  column_to_rownames("Cell")
DIobj.full@meta.data <- temp_metadata

DIobj.500[["RNA"]] <- split(DIobj.500[["RNA"]],
                                f = DIobj.500$SampleGroup)
DIobj.full[["RNA"]] <- split(DIobj.full[["RNA"]],
                                f = DIobj.full$SampleGroup)
DIobj.500 <- SCTransform(
  DIobj.500,
  assay = "RNA",
  vars.to.regress = "percent.mt",
  verbose = FALSE
)
DIobj.full <- SCTransform(
  DIobj.full,
  assay = "RNA",
  vars.to.regress = "percent.mt",
  verbose = FALSE
)
DIobj.500 <- RunPCA(DIobj.500, 
                        features = VariableFeatures(object = DIobj.500))
DIobj.full <- RunPCA(DIobj.full,
                        features = VariableFeatures(object = DIobj.full))
DIobj.full <- RunUMAP(DIobj.full, reduction = "pca", dims = 1:20)
DIobj.full <- FindNeighbors(DIobj.full, 
                                      reduction = "pca",
                            k.param = 30,
                                      dims = 1:20)
DIobj.full <- FindClusters(DIobj.full, 
                                     cluster.name = "uncorrected_clusters")
DimPlot(DIobj.full, reduction = "umap", group.by = "orig.ident")

#####
harmony_object_theta1 <- IntegrateLayers(object = DIobj.full, 
                                        method = HarmonyIntegration,
                                        orig.reduction = "pca", 
                                        new.reduction = "harmony_theta_1",
                                        theta = 0.05)

harmony_object_theta1 <- RunUMAP(harmony_object_theta1,
                   reduction = "harmony_theta_1",
                                dims = 1:20)
harmony_object_theta1 <- FindNeighbors(harmony_object_theta1, reduction = "harmony_theta_1",
                        k.param = 30,
                                dims = 1:20)
harmony_object_theta1 <- FindClusters(harmony_object_theta1, 
                               cluster.name = "harmony_clusters1_k30")
theta1_plot <- DimPlot(harmony_object_theta1, 
                        reduction = "umap", 
                        group.by = "orig.ident") + 
  ggtitle("Harmony integrated data")
theta1_clusters_plot <- DimPlot(harmony_object_theta1, reduction = "umap", 
                                 group.by ="harmony_clusters1_k30") + 
  ggtitle("Harmony data clusters")

theta1_plot + theta1_clusters_plot

FeaturePlot(harmony_object_theta1, features = c("CD79A", "CST3", "CD3D", "HBA1"), 
            reduction = "umap")

library(bluster)
ari.T <- pairwiseRand(DIobj.full$uncorrected_clusters,
                      harmony_object_theta1$harmony_clusters1_k30,
                      mode = "index")

ari.T

harmony_object_theta1 <- JoinLayers(harmony_object_theta1, assay = "RNA")
saveRDS(harmony_object_theta1, "RObjects/DIobj.full.rds")

#####

DIobj.500 <- IntegrateLayers(object = DIobj.500, 
                                        method = HarmonyIntegration,
                                        orig.reduction = "pca", 
                                        new.reduction = "harmony",
                                        theta = 0.1)
DIobj.full <- IntegrateLayers(object = DIobj.full,
                                        method = HarmonyIntegration,
                                        orig.reduction = "pca",
                                        new.reduction = "harmony",
                                        theta = 0.1)

DIobj.500 <- RunUMAP(DIobj.500,
                                reduction = "harmony",
                                dims = 1:15)
DIobj.500 <- JoinLayers(DIobj.500, assay = "RNA")

saveRDS(DIobj.500, file = "../RObjects/DI.500.rds")
