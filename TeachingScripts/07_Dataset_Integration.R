## ----libraries----------------------------------------
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(tidyverse)
library(patchwork)


## ----split-----------------------------------------------------------------------
# Load the data
seurat_object <- readRDS("../RObjects/Filtered.500.rds")
# Split the data by sample group
seurat_object[["RNA"]] <- split(seurat_object[["RNA"]],
                                f = seurat_object$SampleGroup)
seurat_object


## ----reprocess-------------------------------------------------------------------
# Re-process, Seurat will treat each sample group separately
seurat_object <- SCTransform(
  seurat_object,
  assay = "RNA",
  vars.to.regress = "percent.mt",
  verbose = FALSE
)
seurat_object <- RunPCA(seurat_object, 
                         features = VariableFeatures(object = seurat_object))
seurat_object <- RunUMAP(seurat_object,
                         reduction = "pca",
                         dims = 1:15)
seurat_object


## ----uncorrected-----------------------------------------------------------------
# Visualise the uncorrected data
uncorrected_plot <- DimPlot(seurat_object,
                            reduction = "umap",
                            group.by = "SampleName") + 
  ggtitle("Uncorrected data")


## ----cluster---------------------------------------------------------------------
# Run clustering
seurat_object <- FindNeighbors(seurat_object, reduction = "pca",
                                dims = 1:15)
seurat_object <- FindClusters(seurat_object, 
                              cluster.name = "uncorrected_clusters")

uncorrected_clusters_plot <- DimPlot(seurat_object, reduction = "umap", 
                                     group.by ="uncorrected_clusters") + ggtitle("Uncorrected data clusters")
 
uncorrected_plot + uncorrected_clusters_plot                               


## ----cluster_unc_tab-------------------------------------------------------------
# Make a table of the clusters and samples
table(seurat_object$uncorrected_clusters, seurat_object$orig.ident)


## ----harmony---------------------------------------------------------------------
# Run Harmony integration
harmony_object <- IntegrateLayers(object = seurat_object, 
                                  method = HarmonyIntegration,
                                  orig.reduction = "pca", 
                                  new.reduction = "harmony")


## ----harmony_plot----------------------------------------------------------------
# Run UMAP
harmony_object <- RunUMAP(harmony_object, 
                          reduction = "harmony", 
                          dims = 1:15)
# Run clustering
harmony_object <- FindNeighbors(harmony_object, reduction = "harmony",
                                dims = 1:15)
harmony_object <- FindClusters(harmony_object, 
                              cluster.name = "harmony_clusters")
# make the plots
harmony_plot <- DimPlot(harmony_object, 
                        reduction = "umap", 
                        group.by = "SampleName") + 
  ggtitle("Harmony integrated data")
harmony_clusters_plot <- DimPlot(harmony_object, reduction = "umap", 
                                     group.by ="harmony_clusters") + ggtitle("Harmony data clusters")
 
harmony_plot + harmony_clusters_plot     




## ----sample_level----------------------------------------------------------------
# Read in our filtered data
seurat_object_sample <- readRDS("../RObjects/Filtered.500.rds")

# Split the data by sample name
seurat_object_sample[["RNA"]] <- split(seurat_object_sample[["RNA"]],
                                f = seurat_object_sample$SampleName)
# Re-process, Seurat will treat each sample separately
seurat_object_sample <- SCTransform(
  seurat_object_sample,
  assay = "RNA",
  vars.to.regress = "percent.mt",
  verbose = FALSE
)
seurat_object_sample <- RunPCA(seurat_object_sample, 
                         features = VariableFeatures(object = seurat_object_sample))
seurat_object_sample <- RunUMAP(seurat_object_sample,
                         reduction = "pca",
                         dims = 1:15)
uncorrected_plot_sample <- DimPlot(seurat_object_sample,
                            reduction = "umap",
                            group.by = "SampleName") + 
  ggtitle("Uncorrected data - Sample")

seurat_object_sample <- FindNeighbors(seurat_object_sample, 
                                      reduction = "pca",
                                dims = 1:15)
seurat_object_sample <- FindClusters(seurat_object_sample, 
                              cluster.name = "uncorrected_clusters")

uncorrected_clusters_sample_plot <- DimPlot(seurat_object, reduction = "umap", 
                                     group.by ="uncorrected_clusters") + ggtitle("Uncorrected data clusters - Sample")

uncorrected_clusters_sample_plot



## ----harmony_sample--------------------------------------------------------------
# Run Harmony integration
harmony_object_sample <- IntegrateLayers(object = seurat_object_sample,
                                  method = HarmonyIntegration,
                                  orig.reduction = "pca",
                                  new.reduction = "harmony_sample")
harmony_object_sample <- RunUMAP(harmony_object_sample,
                          reduction = "harmony_sample",
                          dims = 1:15)
harmony_object_sample_plot <- DimPlot(harmony_object_sample, 
                                     reduction = "umap", 
                                     group.by = "SampleName") + 
  ggtitle("Harmony integrated data - Sample")

uncorrected_plot_sample + harmony_object_sample_plot

harmony_object_sample <- FindNeighbors(harmony_object_sample, 
                                       reduction = "harmony_sample",
                                dims = 1:15)
harmony_object_sample <- FindClusters(harmony_object_sample, 
                              cluster.name = "harmony_sample_clusters")

harmony_object_sample_clusters_plot <- DimPlot(harmony_object_sample, reduction = "umap", 
                                     group.by ="harmony_sample_clusters") + ggtitle("Harmony integrated data clusters - Sample")


 
uncorrected_clusters_sample_plot + harmony_object_sample_clusters_plot



## ----theta-----------------------------------------------------------------------
# Run Harmony integration with adjusted theta
harmony_object_theta <- IntegrateLayers(object = seurat_object, 
                                  method = HarmonyIntegration,
                                  orig.reduction = "pca", 
                                  new.reduction = "harmony_theta_01",
                                  theta = 0.1)

harmony_object_theta <- RunUMAP(harmony_object_theta,
                          reduction = "harmony_theta_01", 
                          dims = 1:15)
harmony_object_theta_plot <- DimPlot(harmony_object_theta, 
                                     reduction = "umap", 
                                     group.by = "SampleName") + 
  ggtitle("Harmony integrated data with theta = 0.1")

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


## ----barplotUnc------------------------------------------------------------------
#Make a barplot of the clusters and samples
data.frame(Cluster = seurat_object$uncorrected_clusters, Sample = seurat_object$SampleName) %>%
  ggplot(aes(x = Cluster)) +
    geom_bar(aes(fill = Sample), position = "fill") +
    labs(title = "Uncorrected data") +
  scale_y_continuous(labels = scales::percent)


## ----barplotHarm-----------------------------------------------------------------
#Make a barplot of the clusters and samples
data.frame(Cluster = harmony_object_theta$harmony_theta_01_clusters, Sample = harmony_object_theta$SampleName) %>%
  ggplot(aes(x = Cluster)) +
    geom_bar(aes(fill = Sample), position = "fill") +
    labs(title = "Harmony Corrected data") +
  scale_y_continuous(labels = scales::percent)


## ----markers---------------------------------------------------------------------
# CD79A (B cells)
# CST3 (monocytes)
# CD3D (T cells)
# HBA1 (erythrocytes)

FeaturePlot(seurat_object, features = c("CD79A", "CST3", "CD3D", "HBA1"), 
            reduction = "umap")


FeaturePlot(harmony_object_theta, features = c("CD79A", "CST3", "CD3D", "HBA1"), 
            reduction = "umap")



## ----ariSG-----------------------------------------------------------------------
library(bluster)
# Calculate the adjusted Rand index
# Input is a vector of the clusters
ari.SG <- pairwiseRand(seurat_object$uncorrected_clusters,
                    harmony_object$harmony_clusters,
                         mode = "index")

ari.SG


## ----ariSN-----------------------------------------------------------------------
ari.SN <- pairwiseRand(seurat_object_sample$uncorrected_clusters,
                    harmony_object_sample$harmony_sample_clusters,
                         mode = "index")

ari.SN


## ----ariTheta--------------------------------------------------------------------
ari.T <- pairwiseRand(seurat_object$uncorrected_clusters,
                    harmony_object_theta$harmony_theta_01_clusters,
                         mode = "index")

ari.T


## ----rejoin----------------------------------------------------------------------
# Rejoin the data
harmony_object_theta <- JoinLayers(harmony_object_theta, assay = "RNA")


## ----session_info----------------------------------------------------------------
sessionInfo()

