## ----libraries----------------------------------------------------------------
# load packages
library(BiocParallel)
library(scran)
library(scater)
library(miloR)
library(tidyverse)
library(patchwork)


## ----loaddata-----------------------------------------------------------------
# Load data
sce <- readRDS("RObjects/Annotated.full.ETV6.PBMMC.sce.rds")




## ----setobject----------------------------------------------------------------
# check the contents of the object
sce

# plot UMAP done on the batch-corrected data
plotReducedDim(sce, dimred = "UMAP", 
               colour_by = "label", 
               text_by = "label")

# convert to a Milo object

milo <- Milo(sce)

milo


## -----------------------------------------------------------------------------

# Build KNN graph ----


## ----build_graph--------------------------------------------------------------
# add KNN graph to Milo object
milo <- buildGraph(milo, 
                   k = 60, 
                   d = 50, 
                   reduced.dim = "HARMONY", 
                   BPPARAM = MulticoreParam(7))




## ----make_neighbourhoods, message=FALSE---------------------------------------
# sample index cells to define neighbourhoods
milo <- makeNhoods(milo, 
                   prop = 0.1, 
                   k = 48, 
                   d = 50, 
                   reduced_dims = "HARMONY")

# check our object again
milo


## ----plot_neighbourhood_sizes, message=FALSE----------------------------------
# distribution of neighbourhood sizes
plotNhoodSizeHist(milo) +
  geom_vline(xintercept = 100, col = "salmon")


## ----count_cells, message=FALSE-----------------------------------------------
# count cells in each neighbourhood
milo <- countCells(milo, 
                   meta.data = colData(milo),
                   samples = "SampleName")

# Milo now has a counts matrix
head(nhoodCounts(milo))


## ----echo=FALSE---------------------------------------------------------------

# Run DA analysis ----


## ----calcultate_neighbourhood_distance, message=FALSE-------------------------
# calculate distances between neighbourhoods - for p-value correction
# This step takes a few minutes
milo <- calcNhoodDistance(milo, d = 50, reduced.dim = "HARMONY")


## ----milo_design, message=FALSE-----------------------------------------------
# define a table for our model design
sample_info <- unique(colData(milo)[,c("SampleName", "SampleGroup")])
rownames(sample_info) <- sample_info$SampleName

sample_info


## ----test_neighbourhoods, message=FALSE---------------------------------------
# run DA test
da_results <- testNhoods(milo, 
                         design = ~ SampleGroup, 
                         design.df = sample_info, 
                         reduced.dim = "HARMONY")

# results are returned as a data.frame
da_results %>%
  arrange(SpatialFDR) %>%
  head()


## ----plot_DA_pvalue_histogram, message=FALSE----------------------------------
# p-value histogram
ggplot(da_results, aes(PValue)) + 
  geom_histogram(bins = 50)


## ----plot_DA_volcano, message=FALSE-------------------------------------------
# volcano plot
# each point in this plot corresponds to a neighbourhood (not a cell)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point(aes(colour = FDR < 0.1)) +
  geom_hline(yintercept = 1) 


## ----build_neighbourhood_graph, message=FALSE---------------------------------
# build neighbourhood graph embedding
milo <- buildNhoodGraph(milo)


## ----visualise_results, message=FALSE-----------------------------------------
# our original UMAP with our previously annotated cell labels
umap_plot <- plotReducedDim(milo, 
                            dimred = "UMAP", 
                            colour_by = "label", 
                            text_by = "label")

# the neighbourhood map adjusted to match UMAP embedding
nh_graph_plot <- plotNhoodGraphDA(milo, 
                                  da_results, 
                                  layout = "UMAP",
                                  alpha = 0.05)

# the two plots together side-by-side
umap_plot + nh_graph_plot +
  plot_layout(guides="collect")




## ----annotate_nhoods----------------------------------------------------------
# annotate our neighbourhood DA results with our cell labels
da_results <- annotateNhoods(milo, da_results, coldata_col = "label")
head(da_results)


## ----histogram_label_fraction-------------------------------------------------
# histogram of fraction of cells in the neighbourhood with the same label
ggplot(da_results, aes(label_fraction)) + 
  geom_histogram(bins = 50)


## ----excluded_mixed_neighbourhoods, message=FALSE-----------------------------
# add "mixed" label to neighbourhoods with less 70% consistency
da_results$label <- ifelse(da_results$label_fraction < 0.7, 
                           "Mixed", 
                           da_results$label)

head(da_results)


## ----beeswarm_plot_by_label, message=FALSE------------------------------------
# distribution of logFC across neighbourhood labels
plotDAbeeswarm(da_results, group.by = "label")


## ----echo=FALSE---------------------------------------------------------------
sessionInfo()

