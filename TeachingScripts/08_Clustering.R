# Load libraries
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(cluster) # for silhouette function
library(tidyverse)
library(patchwork)

# set default ggplot theme
theme_set(theme_classic())

# Load the data
seurat_object <- readRDS("RObjects/DI.500.rds")

# Confirm the default assay is set to SCT
DefaultAssay(seurat_object)


#### Clustering ####

# construct the neighbour graph, by default it constructs both:
#   - nearest neighbour (NN) graph
#   - shared nearest neighbour (SNN) graph
seurat_object <- FindNeighbors(seurat_object,
                               reduction = "harmony",
                               k.param = 20, # k = 20 is the default
                               dims = 1:15)

# Check the graphs slot
Graphs(seurat_object)

# identify clusters in the nearest neighbour graph (SNN used by default)
seurat_object <- FindClusters(seurat_object,
                              resolution = 0.8, # default is 0.8
                              algorithm = 1, # default is 1 (louvain)
                              cluster.name = "Louvain_k20_res0.8")

# confirm the new column in the meta.data slot
head(seurat_object[[]])

# visualise the clusters on a UMAP plot
DimPlot(seurat_object, 
        reduction = "umap", 
        group.by = "Louvain_k20_res0.8")


#### Exercise: Try different parameter combinations ####
# Experiment with different values of `k` and `resolution`

# Cluster with resolution 0.4 using the previous k value of 20
# YOUR CODE HERE

# Cluster with resolution 1.6 using the previous k value of 20
# YOUR CODE HERE

# Find neighbours with k = 10 and cluster with resolution 0.8
# YOUR CODE HERE

# Find neighbours with k = 40 and cluster with resolution 0.8
# YOUR CODE HERE






#### Exercise: Leiden clustering ####
# Repeat the clustering using the Leiden method (algorithm = 4)

# Find neighbours with k = 20 and cluster with resolution 0.8 using the Leiden method
# YOUR CODE HERE

# Visualise the clusters on a UMAP plot
# YOUR CODE HERE




#### Testing multiple resolutions ####

# FindClusters supports testing multiple resolution values in a single run
seurat_object <- FindClusters(seurat_object,
                              resolution = c(0.4, 0.8, 1.6, 2.0),
                              algorithm = 1)

# Confirm the new columns in the meta.data slot
# named: SCT_snn_res<resolution_value>
head(seurat_object[[]])


#### Testing multiple k ####

# Testing multiple values for k with a for loop
k_values <- c(5, 10, 15, 20, 25)

# For each value of k, run FindNeighbors and FindClusters
for (k in k_values) {
  seurat_object <- FindNeighbors(seurat_object,
                                 reduction = "harmony",
                                 k.param = k,
                                 dims = 1:15)
  
  seurat_object <- FindClusters(seurat_object,
                                resolution = 0.8,
                                algorithm = 1,
                                cluster.name = paste0("Louvain_k", k,
                                                      "_res0.8"))
}
colnames(seurat_object[[]])




#### Assessing cluster behaviour ####

# get Harmony data as we used for clustering
harmony_data <- Embeddings(seurat_object, reduction = "harmony")[, 1:15]

# calculate distance matrix for the harmony data
harmony_dist <- dist(harmony_data)

# get our choice of cluster assignments
# because it's a factor, we convert it to integer
clusters <- seurat_object$Louvain_k20_res0.8
clusters <- as.integer(as.character(clusters))

# calculate silhouette width for each cell
silhouette_width <- silhouette(clusters,
                               harmony_dist) %>%
  # convert to data frame for plotting
  as_tibble()

# look at the top 15 rows
head(silhouette_width, n = 15)

# silhouette width distribution coloured by closest cluster
p1 <- silhouette_width %>%
  # create variable for closest cluster based on silhouette width
  mutate(closest_cluster = ifelse(sil_width > 0, cluster, neighbor)) %>%
  ggplot(aes(x = factor(cluster), 
             y = sil_width, 
             colour = factor(closest_cluster))) +
  ggbeeswarm::geom_quasirandom() +
  labs(x = "Cluster", 
       y = "Silhouette width", 
       colour = "Closest cluster")

# UMAP plot for comparison
p2 <- DimPlot(seurat_object,
              reduction = "umap",
              group.by = "Louvain_k20_res0.8",
              label = TRUE) +
  NoLegend()

# plot them together
p1 + p2


#### Silhouette width for Leiden clustering ####

# In case you didn't complete the previous exercise, make sure to run
# clustering with leiden, k = 20, resolution = 0.8
seurat_object <- FindNeighbors(seurat_object,
                               reduction = "harmony",
                               k.param = 20,
                               dims = 1:15)
seurat_object <- FindClusters(seurat_object,
                              resolution = 0.8,
                              algorithm = 4,
                              cluster.name = "Leiden_k20_res0.8",
                              random.seed = 123) # set seed (warning without)


# get our choice of cluster assignments
clusters_leiden <- seurat_object$Leiden_k20_res0.8
clusters_leiden <- as.integer(as.character(clusters_leiden))

# calculate silhouette width for each cell
silhouette_width_leiden <- silhouette(clusters_leiden,
                                      harmony_dist) %>%
  # convert to data frame for plotting
  as_tibble()

# silhouette width distribution coloured by closest cluster
p1_leiden <- silhouette_width_leiden %>%
  # create variable for closest cluster based on silhouette width
  mutate(closest_cluster = ifelse(sil_width > 0, cluster, neighbor)) %>%
  ggplot(aes(x = factor(cluster), 
             y = sil_width, 
             colour = factor(closest_cluster))) +
  ggbeeswarm::geom_quasirandom() +
  labs(x = "Cluster", 
       y = "Silhouette width", 
       colour = "Closest cluster")

# UMAP plot for comparison
p2_leiden <- DimPlot(seurat_object,
                     reduction = "umap",
                     group.by = "Leiden_k20_res0.8")

# plot them together
p1_leiden + p2_leiden

#### Calculate mean silhouette width for different values of k ####

# function to calculate mean silhouette width for a given clustering
calc_mean_sil <- function(clusters, dist_matrix) {
  silhouette_width <- silhouette(clusters, dist_matrix)
  mean(silhouette_width[, "sil_width"])
}

# function to calculate number of clusters for a given clustering
calc_num_clusters <- function(clusters) {
  length(unique(clusters))
}

# list of k-values we're interested in
k_values <- c(5, 10, 15, 20, 25)

# get the seurat metadata for convenience
seurat_meta <- seurat_object[[]]

# calculate mean silhouette width and number of clusters for each k
# the lapply() function loops through each k value
# and applies the code inside the function 
silhouette_stats <- lapply(k_values, function(k) {
  # column name for the clustering with this k value
  cluster_col <- paste0("Louvain_k", k, "_res0.8")
  
  # get the cluster assignments for this k value
  clusters <- seurat_meta[, cluster_col]
  
  # convert to integer
  clusters <- as.integer(as.character(clusters))
  
  # create a data frame to store the results
  tibble(k = k,
         num_clusters = calc_num_clusters(clusters),
         mean_silhouette_width = calc_mean_sil(clusters, harmony_dist))
}) %>%
  # bind the list into a single data frame
  bind_rows()

# look at the results
silhouette_stats

# Mean silhouette plot
p_mean_sil <- silhouette_stats %>%
  ggplot(aes(x = k, y = mean_silhouette_width)) +
  geom_line(lwd = 2)

# Number of clusters plot
p_num_clusters <- silhouette_stats %>%
  ggplot(aes(x = k, y = num_clusters)) +
  geom_line(lwd = 2)

# plot them together
p_num_clusters + p_mean_sil


#### Exercise: test different parameter combinations ####

# Test a combination of resolution and k values
# For example k = 10, 14, 17, 20, 25
# and resolution = 0.6, 0.8, 1.0, 1.2, 1.4, 1.6

# Test different k values first, clustering with Leiden with resolution 0.8
# YOUR CODE HERE

# Test different clustering resolution values using Leiden with k = 17
# YOUR CODE HERE








#### Finalise the clustering ####

# Add our final choice of clustering to the seurat object
# Create SNN graph with k = 17
seurat_object <- FindNeighbors(seurat_object,
                               reduction = "harmony",
                               k.param = 17,
                               dims = 1:15)

# Cluster cells with resolution = 1 using the Leiden method
seurat_object <- FindClusters(seurat_object,
                              resolution = 1,
                              algorithm = 4,
                              random.seed = 123,
                              cluster.name = "Leiden_k17_res1")

# Visualise the clusters on a UMAP plot
DimPlot(seurat_object, 
        reduction = "umap", 
        group.by = "Leiden_k17_res1")

# Visualise the expression of a marker gene on the UMAP plot and across clusters
fplot <- FeaturePlot(seurat_object, 
                     features = "HBA1", 
                     reduction = "umap")

# visualise the distribution of expression across clusters
vplot <- VlnPlot(seurat_object, features = "HBA1",
                 group.by = "Leiden_k17_res1",
                 pt.size = 0) + 
  NoLegend()

# plot them together
fplot + vplot
