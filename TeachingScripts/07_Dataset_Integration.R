library(Seurat)
library(sctransform)
library(glmGamPoi)
library(tidyverse)
library(patchwork)

############################################################################
# Load the data
seurat_object <- readRDS("RObjects/Filtered.500.rds")
seurat_object@assays$RNA |> dim()

# assess sample groups in the data
unique(seurat_object$SampleGroup)
############################################################################

# Split the data by sample group
seurat_object[["RNA"]] <- split(seurat_object[["RNA"]],
                                f = seurat_object$SampleGroup)

# explore ETV6RUNX1 group counts matrix
seurat_object@assays$RNA$counts.ETV6RUNX1 |> dim()

# explore seurat object after the split
seurat_object
############################################################################

# Group level Re-processing (Normalization and variable gene selection) 
# Seurat will treat each sample group separately
seurat_object <- SCTransform(
  seurat_object,
  assay = "RNA",
  vars.to.regress = "percent.mt",
  verbose = FALSE
)


# Layers in the SCT assay
Layers(seurat_object[["SCT"]])

# Names of the fitted SCT models: expect 4
levels(seurat_object[["SCT"]])
############################################################################

# Run PCA - PCA uses "scale.data" from STC assay
# Scale data = pearson residuals
seurat_object <- RunPCA(seurat_object, 
                         features = VariableFeatures(object = seurat_object))

# explore dimensions of scale.data matrix
seurat_object@assays$SCT$scale.data |> dim()

# Explore PCA
seurat_object@reductions$pca
############################################################################

# Run UMAP on PCA space
seurat_object <- RunUMAP(seurat_object,
                         reduction = "pca",
                         dims = 1:15)

# explore seurat object
seurat_object
############################################################################

# Visualize the uncorrected data
# color by sample
uncorrected_plot <- DimPlot(seurat_object,
                            reduction = "umap",
                            group.by = "SampleName") + 
  ggtitle("Uncorrected data Sample Name")

uncorrected_plot

# Visualize the uncorrected data
# color by SampleGroup
uncorrected_plot2 <- DimPlot(seurat_object,
                            reduction = "umap",
                            group.by = "SampleGroup") + 
  ggtitle("Uncorrected data Sample Group")

uncorrected_plot2

uncorrected_plot + uncorrected_plot2
############################################################################

# Run clustering on PCA space (before integration)

# clustering step 1: generate graph for clustering
seurat_object <- FindNeighbors(seurat_object, reduction = "pca",
                                dims = 1:15)
# explores seurat objectgraphs
seurat_object@graphs

# clustering step 2: generate actual clusters 
seurat_object <- FindClusters(seurat_object, 
                              cluster.name = "uncorrected_clusters")
# explore clusters column in the metadata data frame
seurat_object[[]] |> head()
############################################################################

uncorrected_clusters_plot <- DimPlot(seurat_object, reduction = "umap", 
                                     group.by ="uncorrected_clusters") + ggtitle("Uncorrected data clusters")

# show all three UMAP plots in single window
uncorrected_plot + uncorrected_plot2 + uncorrected_clusters_plot                               

# Sample Group level split: Cluster composition before integration
table(seurat_object$uncorrected_clusters, seurat_object$orig.ident)
############################################################################


############################################################################
# Run Harmony integration
# Harmony uses PCA space
harmony_object <- IntegrateLayers(object = seurat_object, 
                                  method = HarmonyIntegration,
                                  orig.reduction = "pca", 
                                  new.reduction = "harmony")
# Explore harmony space
harmony_object@reductions$harmony@cell.embeddings |> dim()
############################################################################

# Run UMAP on harmony Space (after integration)
harmony_object <- RunUMAP(harmony_object, 
                          reduction = "harmony", 
                          dims = 1:15)
############################################################################

# Run clustering on harmony Space (after integration)
# clustering step 1: generate graph for clustering
harmony_object <- FindNeighbors(harmony_object, reduction = "harmony",
                                dims = 1:15)
# clustering step 2: generate actual clusters 
harmony_object <- FindClusters(harmony_object, 
                              cluster.name = "harmony_clusters")
harmony_object[[]] |> head()
############################################################################
# UMAP plot color by Sample Name
harmony_plot <- DimPlot(harmony_object, 
                        reduction = "umap", 
                        group.by = "SampleName") + 
  ggtitle("Harmony integrated data")

# UMAP plot color by Sample Group
harmony_plot2 <- DimPlot(harmony_object, 
                        reduction = "umap", 
                        group.by = "SampleGroup") + 
  ggtitle("Harmony integrated data Sample Group")

# UMAP plot color by Cluster
harmony_clusters_plot <- DimPlot(harmony_object, reduction = "umap", 
                                     group.by ="harmony_clusters") + ggtitle("Harmony data clusters")

# Plot all thee UMAPs in single window
harmony_plot +  harmony_plot2 + harmony_clusters_plot 

# Sample Group level split: Cluster composition before integration
table(seurat_object$uncorrected_clusters, seurat_object$SampleName) 
# Sample Group level split: Cluster composition after integration
table(harmony_object$harmony_clusters, harmony_object$SampleName)



############################################################################
################ Exercise: 7.2.2 ########################################### 
# 1. Starting with our filtered data, split the object by SampleName.

############################################################################
# 2. Reprocess the data with SCTransform and run dimensionality reduction on the new layers.

############################################################################
# 3. Visualise your uncorrected data.

############################################################################
# 4. Run Harmony integration using the IntegrateLayers function.

############################################################################
# 5. Visualise the results and compare them with the sample group level integration.

############################################################################
############################################################################


# Run Harmony integration with adjusted theta
# theta parameter controls the diversity of data sets within clusters
# Higher theta gives more diverse clusters, lower value has opposite effect
# seurat_object is from Sample Group level split
harmony_object_theta <- IntegrateLayers(object = seurat_object, 
                                  method = HarmonyIntegration,
                                  orig.reduction = "pca", 
                                  new.reduction = "harmony_theta_01",
                                  theta = 0.1)

# Run UMAP on harmony_theta_01 space
harmony_object_theta <- RunUMAP(harmony_object_theta,
                          reduction = "harmony_theta_01", 
                          dims = 1:15)

# UMAP plot color by Sample Name
harmony_object_theta_plot <- DimPlot(harmony_object_theta, 
                                     reduction = "umap", 
                                     group.by = "SampleName") + 
  ggtitle("Harmony integrated data with theta = 0.1")

harmony_object_theta_plot

############################################################################
# cluster the cells by using harmony_theta_01 space  
harmony_object_theta <- FindNeighbors(harmony_object_theta, 
                                       reduction = "harmony_theta_01",
                                dims = 1:15)
harmony_object_theta <- FindClusters(harmony_object_theta,
                              cluster.name = "harmony_theta_01_clusters")


harmony_object_theta_clusters_plot <- DimPlot(harmony_object_theta, 
                                              reduction = "umap", 
                                              group.by ="harmony_theta_01_clusters") + 
  ggtitle("Harmony integrated data clusters with theta = 0.1")

harmony_object_theta_plot + harmony_object_theta_clusters_plot
############################################################################

# Make a barplot of the clusters and samples 
# From Sample Group level split
data.frame(Cluster = seurat_object$uncorrected_clusters, Sample = seurat_object$SampleName) %>%
  ggplot(aes(x = Cluster)) +
    geom_bar(aes(fill = Sample), position = "fill") +
    labs(title = "Uncorrected data") +
  scale_y_continuous(labels = scales::percent)

# Make a barplot of the clusters and samples
data.frame(Cluster = harmony_object_theta$harmony_theta_01_clusters, Sample = harmony_object_theta$SampleName) %>%
  ggplot(aes(x = Cluster)) +
    geom_bar(aes(fill = Sample), position = "fill") +
    labs(title = "Harmony Corrected data") +
  scale_y_continuous(labels = scales::percent)

# CD79A (B cells)
# CST3 (monocytes)
# CD3D (T cells)
# HBA1 (erythrocytes)

FeaturePlot(seurat_object, features = c("CD79A", "CST3", "CD3D", "HBA1"), 
            reduction = "umap") +
  plot_annotation( 'UMAP before correction' )


FeaturePlot(harmony_object_theta, features = c("CD79A", "CST3", "CD3D", "HBA1"), 
            reduction = "umap")+
  plot_annotation( 'UMAP After correction' )


library(bluster)
# Calculate the adjusted Rand index
# Input is a vector of the clusters
# seurat_object from SampleGroup level split
# harmony_object from SampleGroup level split
ari.SG <- pairwiseRand(seurat_object$uncorrected_clusters,
                    harmony_object$harmony_clusters,
                         mode = "index")

ari.SG


# seurat_object_sample from Sample level split
# harmony_object_sample from Sample  level split
ari.SN <- pairwiseRand(seurat_object_sample$uncorrected_clusters,
                    harmony_object_sample$harmony_sample_clusters,
                         mode = "index")

ari.SN


# seurat_object from Sample Group level split
# harmony_object_theta from Sample Group level split with theta = 0.1

ari.T <- pairwiseRand(seurat_object$uncorrected_clusters,
                    harmony_object_theta$harmony_theta_01_clusters,
                         mode = "index")

ari.T

# Rejoin the data
harmony_object_theta <- JoinLayers(harmony_object_theta, assay = "RNA")

sessionInfo()

