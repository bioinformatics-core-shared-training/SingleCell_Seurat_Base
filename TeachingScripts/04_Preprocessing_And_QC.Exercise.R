sampleinfo <- read_tsv("Data/sample_sheet.tsv")



temp_meta.data <- seurat_object@meta.data %>%
  rownames_to_column("Cell") %>%
  separate(col = "orig.ident", into = c("SampleGroup"), sep = "-", remove = FALSE) %>%
  mutate(SampleName = orig.ident) %>%
  column_to_rownames("Cell")

seurat_object@meta.data <- temp_meta.data
head(seurat_object@meta.data)

seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

# colours <- c("cyan3","cyan3","cyan3","cyan3","darkgoldenrod1","darkgoldenrod1","blue","blue","blue","lightgreen","lightgreen")
# 
# VlnPlot(seurat_object,
#         features = c("nCount_RNA"),
#         cols = colours,
#         layer = "counts",
#         pt.size=0) +
#   ggtitle("Total count of UMIs")







new_metadata <- metadata %>%
  rownames_to_column("Cell") %>%
  mutate(Keep = ifelse(Cell %in% all_keep_cells, "Keep", "Remove")) %>%
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



filtered_seurat_object <- subset(seurat_object, cells = all_keep_cells)

sessionInfo()
