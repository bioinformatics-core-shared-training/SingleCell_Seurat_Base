# load the libraries
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(tidyverse)
library(patchwork)


#### Load data ####

# load in the data
seurat_object <- readRDS("RObjects/Clustered.500.rds")

# visualise UMAP with cluster labels
DimPlot(seurat_object,
        reduction = "umap",
        group.by = "leiden_cluster",
        label = TRUE,
        label.size = 5) +
  NoLegend()

# visualise known monocyte marker gene CST3
FeaturePlot(seurat_object,
            reduction = "umap",
            features = "CST3")


#### Identifying cluster marker genes ####

# prepare the data for marker gene detection
seurat_object <- PrepSCTFindMarkers(seurat_object,
                                    assay = "SCT")

# Compare cells in cluster 2 against all other cells
markers_cluster2 <- FindMarkers(seurat_object,
                                ident.1 = 2,
                                ident.2 = NULL,
                                assay = "SCT",
                                group.by = "leiden_cluster")
head(markers_cluster2)

# filter for genes that:
#  - are expressed in at least 70% of cells in cluster 2
#  - and have a log fold change greater than 1
markers_cluster2_filtered <- markers_cluster2 %>%
  filter(pct.1 > 0.70 & abs(avg_log2FC) > 1)

head(markers_cluster2_filtered)

# visualise the expression of one of the marker genes on the UMAP
fplot <- FeaturePlot(seurat_object, 
                     features = "SKAP1",
                     label = TRUE,
                     reduction = "umap")

# As a violin plot
vplot <- VlnPlot(seurat_object, features = "SKAP1",
                 group.by = "leiden_cluster",
                 pt.size = 0) + 
  NoLegend()

# SKAP1 is a TCR gene
# It also appears in cluster 11 next to it
fplot + vplot + plot_layout(ncol = 2)

# find markers for all clusters automatically
markers_all <- FindAllMarkers(seurat_object,
                              assay = "SCT",
                              group.by = "leiden_cluster")

head(markers_all)

# filter the results for cluster 2
# we should get the same results as before
markers_all_cluster2 <- markers_all %>%
  filter(cluster == 2 & pct.1 > 0.70 & abs(avg_log2FC) > 1)

head(markers_all_cluster2)


#### Exercise: cluster 13 markers ####

# Use either the FindMarkers() function or the FindAllMarkers() function
# to identify potential marker genes for cluster 13
# YOUR CODE HERE








#### Heatmap visualisations ####

# top 20 marker genes for cluster 13
top_markers_cluster13 <- markers_cluster13_filtered %>%
  slice_max(n = 20, order_by = abs(avg_log2FC)) %>%
  rownames()

# heatmap of expression for these genes
# we change default colour scale
# and remove colour legend (clusters)
DoHeatmap(seurat_object,
          features = top_markers_cluster13,
          group.by = "leiden_cluster") +
  scale_fill_viridis_c() +
  guides(colour = "none")

# dot plot of expression for these genes
# we change default colour scale
# and rotate x-axis labels for better visibility
DotPlot(seurat_object, 
        features = top_markers_cluster13,
        group.by = "leiden_cluster") +
  scale_colour_viridis_c() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))


#### Manual cell labelling ####

# loop through all marker genes 
# and extract top-ranked gene names for each cluster
top_markers_all <- lapply(unique(markers_all$cluster), function(x) {
  markers_all %>%
    filter(cluster == x) %>%
    filter(pct.1 > 0.70, abs(avg_log2FC) > 1) %>%
    slice_max(n = 15, order_by = abs(avg_log2FC)) %>%
    pull(gene)
})

# we get several known markers of immune cells
top_markers_all

# cell type specific genes
known_genes <- c(
  "HBA1", # erythrocytes
  "CST3", # monocytes
  "CD3E", # T cells
  "NKG7", # NK T cells
  "CD79A", # B cells
  "MS4A1" # CD20 B cells
)

# violin plot
VlnPlot(seurat_object,
        features = known_genes,
        group.by = "leiden_cluster",
        pt.size = 0) + 
  NoLegend() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))

# add new labels to the Idents column in the metadata
seurat_object$Idents <- recode_values(
  seurat_object$leiden_cluster,
  "1" ~ "B (c1)",
  "2" ~ "T (c2)",
  "3" ~ "B (c3)",
  "4" ~ "B (c4)",
  "5" ~ "T (c5)",
  "6" ~ "B (c6)",
  "7" ~ "B (c7)",
  "8" ~ "B (c8)",
  "9" ~ "CD20+ B (c9)",
  "10" ~ "T (c10)",
  "11" ~ "NK T (c11)",
  "12" ~ "B (c12)",
  "13" ~ "Monocytes (c13)",
  "14" ~ "Erythrocytes (c14)",
  "15" ~ "Erythrocytes (c15)",
  "16" ~ "Unknown (c16)"
)

# set the Idents to the new labels
Idents(seurat_object) <- "Idents"

# visualise UMAP with new labels
DimPlot(seurat_object,
        reduction = "umap",
        group.by = "Idents",
        label = TRUE,
        label.size = 5) +
  NoLegend()
