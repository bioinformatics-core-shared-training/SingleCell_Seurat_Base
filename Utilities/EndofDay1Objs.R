library(Seurat)
library(sctransform)
library(ggplot2)
library(tidyverse)
set.seed = 123 
setwd("SingleCell_Seurat_Base/Utilities")

sampleinfo <- read_tsv("../Data/sample_sheet.tsv")

samples <- sampleinfo$Sample
list_of_files <- str_c("../Data/CellRanger_Outputs/", 
                       samples, 
                       "/outs/filtered_feature_bc_matrix")
names(list_of_files) <- sampleinfo$SampleName

expression_matrix <- Read10X(data.dir = list_of_files)
seurat_object = CreateSeuratObject(counts = expression_matrix, min.cells = 1)

orig <- seurat_object

seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

metadata <- seurat_object@meta.data

keep <- c()
for (sam in sampleinfo$SampleName) {
  cells <- metadata %>%
    filter(orig.ident == sam) %>%
    filter(nFeature_RNA > (median(nFeature_RNA) - 2 * mad(nFeature_RNA)),
           nCount_RNA > (median(nCount_RNA) - 2 * mad(nCount_RNA)),
           percent.mt < (median(percent.mt) + 2 * mad(percent.mt))) %>%
    rownames_to_column("Cell") %>%
    pull(Cell)
  keep <- c(keep,cells)
}

new_metadata <- metadata %>%
  rownames_to_column("Cell") %>%
  mutate(Keep = ifelse(Cell %in% keep, "Keep", "Remove")) %>%
  column_to_rownames("Cell")

seurat_object@meta.data <- new_metadata

VlnPlot(seurat_object, 
        features = c("nFeature_RNA"), 
        cols = rep(c("white"), each = 11),
        layer = "counts", 
        pt.size=0) +
  geom_point(mapping = aes(color = seurat_object@meta.data$Keep), size = 0.5) + theme(legend.position = 'none')

VlnPlot(seurat_object, 
        features = c("nCount_RNA"), 
        cols = rep(c("white"), each = 11),
        layer = "counts", 
        pt.size=0) +
  geom_point(mapping = aes(color = seurat_object@meta.data$Keep), size = 0.5) + theme(legend.position = 'none')

VlnPlot(seurat_object, 
        features = c("percent.mt"), 
        cols = rep(c("white"), each = 11),
        layer = "counts", 
        pt.size=0) +
  geom_point(mapping = aes(color = seurat_object@meta.data$Keep), size = 0.5) + theme(legend.position = 'none')

seurat_object_filtered <- subset(seurat_object, cells = keep)

saveRDS(seurat_object_filtered, file = "../RObjects/Filtered.full.rds")

downsampled_seurat <- subset(seurat_object_filtered, downsample = 500)

saveRDS(downsampled_seurat, file = "../RObjects/Filtered.500.rds")
