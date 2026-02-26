## ----libraries, message=FALSE, warning=FALSE--------------------------------
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)


## ----load data--------------------------------------------------------------
# Load the Seurat object
seurat_object <- readRDS("RObjects/Annotated.full.ETV6.PBMMC.rds")



## ----convert to sce, warning=FALSE------------------------------------------
# Convert the Seurat object to a SingleCellExperiment object
sce_object <- as.SingleCellExperiment(seurat_object)


## ----sce--------------------------------------------------------------------
sce_object


## ----check colData----------------------------------------------------------
head(colData(sce_object))


## ----plotumap---------------------------------------------------------------
library(scater)
plotReducedDim(sce_object, "UMAP", colour_by = "leiden_cluster")
DimPlot(seurat_object, reduction = "umap", group.by = "leiden_cluster")

