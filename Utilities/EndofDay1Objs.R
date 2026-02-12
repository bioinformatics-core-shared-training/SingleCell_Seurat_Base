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

temp_metadata <- seurat_object@meta.data %>%
  rownames_to_column("Cell") %>%
  mutate(SampleGroup = str_remove(orig.ident, "-.*")) %>% 
  mutate(SampleName = orig.ident) %>%
  column_to_rownames("Cell")
seurat_object@meta.data <- temp_metadata

orig <- seurat_object

seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

metadata <- seurat_object@meta.data

sample_name_keep_cells <- metadata %>%
  rownames_to_column("Cell") %>%
  group_by(SampleName) %>%
  filter(nFeature_RNA > (median(nFeature_RNA) - 2 * mad(nFeature_RNA)),
         nCount_RNA > (median(nCount_RNA) - 2 * mad(nCount_RNA)),
         percent.mt < (median(percent.mt) + 2 * mad(percent.mt))) %>%
  ungroup() %>%
  pull(Cell)

nrow(metadata) - length(sample_name_keep_cells)

filtered_seurat_object <- subset(seurat_object, 
                                 cells = sample_name_keep_cells)

saveRDS(filtered_seurat_object, file = "../RObjects/Filtered.full.rds")

downsampled_seurat <- subset(filtered_seurat_object, downsample = 500)

downsampled_seurat <- subset(downsampled_seurat, 
                             features = rownames(downsampled_seurat)[Matrix::rowSums(downsampled_seurat[["RNA"]]$counts) > 0])

saveRDS(downsampled_seurat, file = "../RObjects/Filtered.500.rds")





