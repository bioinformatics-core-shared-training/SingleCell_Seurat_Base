library(Seurat)
library(SingleCellExperiment)
library(tidyverse)

# Load the Seurat object
seurat_object <- readRDS("RObjects/Annotated.full.ETV6.PBMMC.rds")


# Convert the Seurat object to a SingleCellExperiment object
sce_object <- as.SingleCellExperiment(seurat_object)

sce_object

head(colData(sce_object))

library(scater)
plotReducedDim(sce_object, "UMAP", colour_by = "leiden_cluster")
DimPlot(seurat_object, reduction = "umap", group.by = "leiden_cluster")

colLabels(sce_object) <- Idents(seurat_object)
