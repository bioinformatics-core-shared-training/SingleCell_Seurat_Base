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

#######

seurat_object <- readRDS("RObjects/DIobj.full.rds")


k_values <- c(39, 40, 41,42, 43,44,45,46,47,48,49,50,51)
for (k in k_values) {
  seurat_object <- FindNeighbors(seurat_object,
                                 reduction = "harmony_theta_1", 
                                 k.param = k, 
                                 dims = 1:20)
  seurat_object <- FindClusters(seurat_object,
                                resolution = 0.8,
                                algorithm = 4,
                                random.seed = 123,
                                cluster.name = paste0("Leiden_k", k, "_res0.8"))
}
library(cluster)
pca_data <- Embeddings(seurat_object, reduction = "pca")
mean_silhouette_widths <- sapply(k_values, function(k) {
  clusters <- seurat_object[[paste0("Leiden_k", k, "_res0.8")]] %>%
    pull()# get cluster assignments for each k
  silhouette_width <- silhouette(as.numeric(clusters), dist(pca_data)) # calculate silhouette width
  mean(silhouette_width[, "sil_width"]) # calculate mean silhouette width
})

silhouette_df <- data.frame(
  k = k_values,
  mean_silhouette_width = mean_silhouette_widths
)
silplot <- silhouette_df %>% 
  ggplot(aes(x = k, y = mean_silhouette_width)) +
  geom_line(lwd=2)
num_clusters <- sapply(k_values, function(k) {
  # get number of clusters for each k
  seurat_object[[paste0("Leiden_k", k, "_res0.8")]] %>%
    pull() %>%
    unique() %>%
    length()
})
clusters_df <- data.frame(
  k = k_values,
  num_clusters = num_clusters
)
cluplot <- clusters_df %>%
  ggplot(aes(x = k, y = num_clusters)) +
  geom_line(lwd=2)

cluplot + silplot

# k=48
seurat_object <- FindNeighbors(seurat_object,
                               reduction = "harmony_theta_1", 
                               k.param = 48, 
                               dims = 1:20)
res_values <- c(0.6,0.7,0.8,0.9)
for (res in res_values) {
  seurat_object <- FindClusters(seurat_object,
                                resolution = res,
                                algorithm = 4,
                                random.seed = 123,
                                cluster.name = paste0("Leiden_k48_res", res))
}
mean_silhouette_widths_res <- sapply(res_values, function(res) {
  clusters <- seurat_object[[paste0("Leiden_k48_res", res)]] %>%
    pull()
  silhouette_width <- silhouette(as.numeric(clusters), dist(pca_data)) 
  mean(silhouette_width[, "sil_width"]) 
})

silhouette_df_res <- data.frame(
  resolution = res_values,
  mean_silhouette_width = mean_silhouette_widths_res
)

silplot_res <- silhouette_df_res %>% 
  ggplot(aes(x = resolution, y = mean_silhouette_width)) +
  geom_line(lwd=2)

silplot_res

num_clusters_res <- sapply(res_values, function(res) {
  # get number of clusters for each resolution
  seurat_object[[paste0("Leiden_k48_res", res)]] %>%
    pull() %>%
    unique() %>%
    length()
})
clusters_df_res <- data.frame(
  res = res_values,
  num_clusters = num_clusters_res
)
cluplot_res <- clusters_df_res %>%
  ggplot(aes(x = res, y = num_clusters)) +
  geom_line(lwd=2)

cluplot_res + silplot_res

# stick with 0.8