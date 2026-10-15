############################################################################
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
############################################################################
# Load the Seurat object
seurat_object <- readRDS("RObjects/Annotated.full.ETV6.PBMMC.rds")

# explore seurat object
seurat_object

seurat_object[[]] |> 
  head()
table(seurat_object$SampleName)
table(seurat_object$SampleGroup)
############################################################################
# Convert the Seurat object to a SingleCellExperiment object
sce_object <- as.SingleCellExperiment(seurat_object)

# Explore sce object
sce_object

head(colData(sce_object))
reducedDim(sce_object, 'PCA')
###########################################################################
library(scater)

plotReducedDim(sce_object, "UMAP", colour_by = "leiden_cluster")

DimPlot(seurat_object, reduction = "umap", group.by = "leiden_cluster")

colData(sce_object) |> 
  head()

colLabels(sce_object) <- Idents(seurat_object)

colData(sce_object) |> 
  head()
###########################################################################