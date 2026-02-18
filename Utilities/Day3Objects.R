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

seurat_object_500 <- PrepSCTFindMarkers(seurat_object_500, assay = "SCT")

temp_metadata <- seurat_object_500@meta.data %>%
  mutate(Idents = case_when(
    leiden_cluster == 1 ~ "B (c1)",
    leiden_cluster == 2 ~ "T (c2)",
    leiden_cluster == 3 ~ "B (c3)",
    leiden_cluster == 4 ~ "B (c4)",
    leiden_cluster == 5 ~ "T (c5)",
    leiden_cluster == 6 ~ "B (c6)",
    leiden_cluster == 7 ~ "B (c7)",
    leiden_cluster == 8 ~ "B (c8)",
    leiden_cluster == 9 ~ "CD20+ B (c9)",
    leiden_cluster == 10 ~ "T (c10)",
    leiden_cluster == 11 ~ "NK T (c11)",
    leiden_cluster == 12 ~ "B (c12)",
    leiden_cluster == 13 ~ "Monocytes (c13)",
    leiden_cluster == 14 ~ "Erythrocytes (c14)",
    leiden_cluster == 15 ~ "Erythrocytes (c15)",
    leiden_cluster == 16 ~ "Unknown (c16)"
  ))

seurat_object_500 <- AddMetaData(seurat_object_500, metadata = temp_metadata$Idents, col.name = "Idents")

Idents(seurat_object_500) <- "Idents"
saveRDS(seurat_object_500, file = "../RObjects/Annotated.500.rds")

