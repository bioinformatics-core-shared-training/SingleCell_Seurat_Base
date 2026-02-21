## ----echo=FALSE-----------------------------------------------------------------------
# Setup ----


## ----libraries, message=FALSE, warning=FALSE------------------------------------------
library(Seurat)
library(ggplot2)
library(tidyverse)
library(DESeq2)


## ----load-data------------------------------------------------------------------------
# Load the data
seurat_object <- readRDS("RObjects/Annotated.full.ETV6.PBMMC.rds")
# plot umap
DimPlot(seurat_object,
        reduction = "umap",
        group.by = "Idents",
        label = TRUE,
        label.size = 5) +
  NoLegend()


## ----echo=FALSE-----------------------------------------------------------------------

# Create pseudo-bulk samples -----


## ----table_label_sample---------------------------------------------------------------
# tabulate number of cells per label + sample combination
table(seurat_object$Idents, seurat_object$SampleName)


## ----pseudo-bulk-aggregation----------------------------------------------------------
# Aggregate counts by sample and celltype
pseudo_bulk <- AggregateExpression(seurat_object,
                                  group.by = c("Idents", "SampleName"),
                                  assays = "RNA",
                                  return.seurat = TRUE)
# check the new assay
pseudo_bulk


## ----check-pseudo-bulk----------------------------------------------------------------
head(Cells(pseudo_bulk))


## ----check-pseudo-bulk-meta-----------------------------------------------------------
head(pseudo_bulk@meta.data)


## ----add-sample-group-----------------------------------------------------------------
# add SampleGroup variable to the metadata of the pseudo-bulk object
temp_metadata <- pseudo_bulk@meta.data %>%
  mutate(SampleGroup = str_remove(SampleName, "-.*")) %>%
  mutate(Cluster = Idents)
pseudo_bulk@meta.data <- temp_metadata
# check the new metadata
head(pseudo_bulk@meta.data)


## ----filter-pseudo-bulk---------------------------------------------------------------
# filter out pseudosamples with less than 20 cells
# determine which pseudosamples have at least than 20 cells
ps_keep <- table(seurat_object$Idents, seurat_object$SampleName) %>%
  as.data.frame() %>%
  filter(Freq >= 20) %>%
  mutate(PseudoSample = str_c(Var1, "_", Var2)) %>%
  pull(PseudoSample)
# subset our pseudo-bulk object to keep only the pseudosamples with at least 20 cells
pseudo_bulk <- subset(pseudo_bulk, cells = ps_keep)


## ----filter-genes-pseudo-bulk---------------------------------------------------------
# filter out genes with less than 5 counts across all pseudosamples
# determine which genes have at least 5 counts across all pseudosamples
genes_keep <- rowSums(pseudo_bulk@assays$RNA$counts) >= 5
# subset our pseudo-bulk object to keep only the genes with at least 5 counts across
pseudo_bulk <- subset(pseudo_bulk, features = names(genes_keep)[genes_keep])


## ----deseq2---------------------------------------------------------------------------
Idents(pseudo_bulk) <- pseudo_bulk$SampleGroup

# Run DESeq2 on pseudobulk samples
de <- FindMarkers(pseudo_bulk,
                  ident.1 = "ETV6RUNX1",
                  ident.2 = "PBMMC",
                  test.use = "DESeq2")
head(de)



## ----add-cluster-sample-group---------------------------------------------------------
# add cluster information to the Idents of our pseudo-bulk object
pseudo_bulk$Cluster_SampleGroup <- str_c(pseudo_bulk$Cluster,
                                         pseudo_bulk$SampleGroup, sep = "_")
Idents(pseudo_bulk) <- pseudo_bulk$Cluster_SampleGroup


## ----pca-pseudo-bulk------------------------------------------------------------------
# PCA plot of pseudo-bulk samples
c6 <- subset(pseudo_bulk, subset = Cluster == "B (c6)")
c6_counts <-  GetAssayData(c6, assay = "RNA", layer = "counts")
pca <- prcomp(t(c6_counts))
pca_df <- data.frame(pca$x, SampleGroup = c6$SampleGroup)
ggplot(pca_df, aes(x = PC1, y = PC2, color = SampleGroup)) +
  geom_point() +
  theme_classic()


## ----deseq2-cluster-------------------------------------------------------------------
# Run DESeq2 on pseudobulk samples for cluster 6
de_cluster6 <- FindMarkers(pseudo_bulk,
                           ident.1 = "B (c6)_ETV6RUNX1",
                           ident.2 = "B (c6)_PBMMC",
                           test.use = "DESeq2")
head(de_cluster6)


## ----ma-plot--------------------------------------------------------------------------
# MA plot of DE results for cluster 6
# calculate mean expression for each gene
mean_expression <- rowMeans(c6_counts)
# add mean expression to DE results
de_cluster6$mean_expression <- mean_expression[rownames(de_cluster6)]
# plot MA plot
ggplot(de_cluster6, aes(x = log2(mean_expression), y = avg_log2FC)) +
  geom_point() +
  theme_classic() +
  xlab("Log2 mean expression") +
  ylab("Log2 fold change") +
  ggtitle("MA plot for cluster 6")



## ----pvalue-histogram-----------------------------------------------------------------
# p-value histogram for DE results for cluster 6
hist(de_cluster6$p_val, breaks = 50, 
     main = "P-value histogram for cluster 6", 
     xlab = "P-value")


## ----filter-de-cluster6---------------------------------------------------------------
# filter the DE results for cluster 6
de_cluster6_sig <- de_cluster6 %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) >1)
dim(de_cluster6_sig)


## ----expression-plots-----------------------------------------------------------------
# expression plots for one of the top differentially expressed genes
top_gene <- "MDK"
# counts for the top gene
top_gene_counts <- GetAssayData(pseudo_bulk, 
                                assay = "RNA", 
                                layer = "counts")[top_gene, ]
# plot expression of the top gene across the clusters and sample groups
top_gene_df <- data.frame(Expression = top_gene_counts,
                          Cluster = pseudo_bulk$Cluster,
                          SampleGroup = pseudo_bulk$SampleGroup,
                          SampleName = pseudo_bulk$SampleName)
ggplot(top_gene_df, aes(x = SampleName, y = Expression, color = SampleGroup
)) +
  geom_point() +
  xlab("Cluster") +
  ylab("Expression of MDK") +
  ggtitle("Expression of MDK across clusters and sample groups") +
  facet_wrap(~ Cluster) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## ----umapMDK--------------------------------------------------------------------------
FeaturePlot(seurat_object, features = top_gene, split.by = "SampleGroup")




## ----session-info---------------------------------------------------------------------
sessionInfo()

