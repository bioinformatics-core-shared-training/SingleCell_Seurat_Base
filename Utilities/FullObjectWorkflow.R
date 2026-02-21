library(Seurat)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
library(tidyverse)
set.seed = 123 
options(future.globals.maxSize = 3.0 * 1e9)

# QC

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

# Normalisation

seurat_object[["RNA"]] <- split(seurat_object[["RNA"]],
                                f = seurat_object$SampleGroup)

seurat_object <- SCTransform(
  seurat_object,
  assay = "RNA",
  vars.to.regress = "percent.mt",
  verbose = FALSE
)
  
seurat_object <- RunPCA(seurat_object, 
                        features = VariableFeatures(object = seurat_object))

saveRDS(seurat_object, file = "../RObjects/SCT.full.rds")

# Integration 

harmony_object <- IntegrateLayers(object = seurat_object, 
                                         method = HarmonyIntegration,
                                         orig.reduction = "pca", 
                                         new.reduction = "harmony",
                                         theta = 0.05)
  
harmony_object <- RunUMAP(harmony_object, 
                          reduction = "harmony", 
                          dims = 1:20)
harmony_object <- JoinLayers(harmony_object, assay = "RNA")

saveRDS(harmony_object, file = "../RObjects/DI.full.rds")

# Clustering
harmony_object <- readRDS("../RObjects/DI.full.rds")
harmony_object <- FindNeighbors(harmony_object,
                                reduction = "harmony",
                                k.param = 48,
                                dims = 1:20)
harmony_object <- FindClusters(harmony_object,
                              resolution = 0.8,
                              algorithm = 4,
                              random.seed = 123,
                              cluster.name = "leiden_cluster")
saveRDS(harmony_object, file = "../RObjects/Clustered.full.rds")

# Annotation

seurat_object <- PrepSCTFindMarkers(harmony_object, assay = "SCT")

known_genes <- c(
  "HBA1", # erythrocytes
  "CST3", # monocytes
  "CD3E", # T cells
  "NKG7", # NK T cells
  "CD79A",  # B cells
  "MS4A1" # CD20 B cells
)

VlnPlot(seurat_object,
        features = known_genes,
        group.by = "leiden_cluster",
        pt.size = 0) + NoLegend() + scale_x_discrete(guide =
                                                       guide_axis(angle = 45))
FeaturePlot(seurat_object,
                            features = known_genes,
                            order = TRUE,
                            label = TRUE,
                            repel = TRUE) + NoLegend()
markers_all <- FindAllMarkers(seurat_object,
                              assay = "SCT",
                              group.by = "leiden_cluster")
saveRDS(markers_all, file = "../RObjects/markers_all.rds")
top_markers_all <- lapply(unique(markers_all$cluster), function(x) {
  markers_all %>%
    filter(cluster == x) %>%
    filter(pct.1 > 0.70, abs(avg_log2FC) > 1) %>%
    top_n(n = 15, wt = abs(avg_log2FC)) %>%
    rownames()
})

# examining this list reveals several known markers of immune cells
top_markers_all

temp_metadata <- seurat_object@meta.data %>%
  mutate(Idents = case_when(
    leiden_cluster == 1 ~ "B (c1)",
    leiden_cluster == 2 ~ "T (c2)",
    leiden_cluster == 3 ~ "B (c3)",
    leiden_cluster == 4 ~ "B (c4)",
    leiden_cluster == 5 ~ "B (c5)",
    leiden_cluster == 6 ~ "B (c6)",
    leiden_cluster == 7 ~ "B (c7)",
    leiden_cluster == 8 ~ "Erythrocytes (c8)",
    leiden_cluster == 9 ~ "Unknown (c9)",
    leiden_cluster == 10 ~ "CD20+ B (c10)",
    leiden_cluster == 11 ~ "NK T (c11)",
    leiden_cluster == 12 ~ "B (c12)",
    leiden_cluster == 13 ~ "T (c13)",
    leiden_cluster == 14 ~ "Monocytes (c14)",
    leiden_cluster == 15 ~ "Erythrocytes (c15)",
    leiden_cluster == 16 ~ "B (c16)",
    leiden_cluster == 17 ~ "Erythrocytes (c17)",
    leiden_cluster == 18 ~ "T (c18)",
    leiden_cluster == 19 ~ "T (c19)",
    leiden_cluster == 20 ~ "T (c20)",
    leiden_cluster == 21 ~ "Monocytes (c21)",
    leiden_cluster == 22 ~ "B (c22)",
    leiden_cluster == 23 ~ "B (c23)"
    
  ))

seurat_object <- AddMetaData(seurat_object, metadata = temp_metadata$Idents, col.name = "Idents")

Idents(seurat_object) <- "Idents"

saveRDS(seurat_object, file = "../RObjects/Annotated.full.rds")

### divide to ETV6RUNX1 vs PBMMC and PRET vs HDD
#seurat_object <- readRDS("../RObjects/Annotated.full.rds")
seurat_object_EP <- subset(seurat_object, 
                        subset = SampleGroup == "ETV6RUNX1" | SampleGroup == "PBMMC")
# drop factor levels to avoid issues with downstream analyses
seurat_object_EP@meta.data$orig.ident <- droplevels(seurat_object_EP@meta.data$orig.ident)
levels(seurat_object_EP@meta.data$orig.ident)
seurat_object_EP@meta.data$SampleName <- droplevels(seurat_object_EP@meta.data$SampleName)
levels(seurat_object_EP@meta.data$SampleName)

saveRDS(seurat_object_EP, "../RObjects/Annotated.full.ETV6.PBMMC.rds")
library(SingleCellExperiment)
sce_object_EP <- as.SingleCellExperiment(seurat_object_EP)
colLabels(sce_object_EP) <- Idents(seurat_object_EP)
saveRDS(sce_object_EP, "../RObjects/Annotated.full.ETV6.PBMMC.sce.rds")


seurat_object_HP <- subset(seurat_object,
                           subset = SampleGroup == "PRET" | SampleGroup == "HHD")
# drop factor levels to avoid issues with downstream analyses

seurat_object_HP@meta.data$orig.ident <- droplevels(seurat_object_HP@meta.data$orig.ident)
levels(seurat_object_HP@meta.data$orig.ident)
seurat_object_HP@meta.data$SampleName <- droplevels(seurat_object_HP@meta.data$SampleName)
levels(seurat_object_HP@meta.data$SampleName)

saveRDS(seurat_object_HP, "../RObjects/Annotated.full.PRET.HHD.rds")
sce_object_HP <- as.SingleCellExperiment(seurat_object_HP)
colLabels(sce_object_HP) <- Idents(seurat_object_HP)
saveRDS(sce_object_HP, "../RObjects/Annotated.full.PRET.HHD.sce.rds")