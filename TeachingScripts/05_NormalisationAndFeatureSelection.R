# Load the libraries we will need for this practical
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(tidyverse)
library(patchwork)
use("SparseArray", c("rowMeans", "rowVars", "colMeans", "colVars"))

# set default ggplot theme
theme_set(theme_classic())

# Read preprocessed and filtered data object
seurat_object <- readRDS("RObjects/Filtered.500.rds")
seurat_object

# Subset seurat object to include only PBMMC-1 sample for demonstration purposes
seurat_pbmmc1 <- subset(seurat_object,
                        subset = SampleName == "PBMMC-1")

# Remove genes that are not expressed in any of the 500 cells in this sample
pbmmc1_expressed_genes <- rowSums(GetAssayData(seurat_pbmmc1, assay = "RNA", layer = "counts")) > 0
seurat_pbmmc1 <- seurat_pbmmc1[pbmmc1_expressed_genes, ]


#### Why normalise? ####

# Reason 1: total UMI differences between cells

# Plot the total UMI counts across cells
# Get the raw counts data
GetAssayData(seurat_pbmmc1, assay = "RNA", layer = "counts") %>%
  # calculate the total UMI counts for each cell (column sums)
  colSums() %>%
  # convert the named vector to a data frame with cell names and UMI counts
  enframe(name = "cell", value = "total_umis") %>%
  # reorder the cells by their total UMI counts for plotting
  mutate(Cell = fct_reorder(cell, total_umis)) %>%
  # make a bar plot of the total UMI counts per cell
  ggplot(aes(x = Cell, y = total_umis)) +
  geom_col() +
  labs(x = "Cell",
       y = "Total cell UMI counts",
       title = "PBMMC-1: Before Normalization") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color = "firebrick")
  )


# Reason 2: mean-variance relationship across genes

# Get the raw counts data from the object
raw_cts <- GetAssayData(seurat_pbmmc1, assay = "RNA", layer = "counts")

# Calculate summary statistics for each gene
gene_raw_stats <- tibble(gene = rownames(raw_cts),
                         mean = rowMeans(raw_cts),
                         variance = rowVars(raw_cts))

# Plot the mean-variance relationship for the raw counts data
ggplot(gene_raw_stats, aes(x = mean, y = variance)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(intercept = 0, slope = 1,
              linewidth = 1, colour = "red") +
  labs(x = "Mean expression (log scale)",
       y = "Variance (log scale)",
       title = "PBMMC-1: Raw Counts - Mean vs Variance")


#### Shifted-log transformation ####

# The default normalisation method (shifted-log transformation)
# 1. counts of each gene are divided by total cell UMIs
# 2. multiplied by a scale factor (default 10,000)
# 3. log-transformed with a pseudocount of 1
seurat_pbmmc1 <- NormalizeData(seurat_pbmmc1)

# A new layer named `data` is added to the RNA assay
seurat_pbmmc1

# Get the log-normalised data from the object
lognorm_cts <- GetAssayData(seurat_pbmmc1, assay = "RNA", layer = "data")
gene_lognorm_stats <- tibble(gene = rownames(lognorm_cts),
                             mean = rowMeans(lognorm_cts),
                             variance = rowVars(lognorm_cts))

# Plot the mean-variance relationship for the log-normalised data
# improved for highly expressed genes
# but still heteroscedastic for lowly expressed genes
ggplot(gene_lognorm_stats, aes(x = mean, y = variance)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1,
              linewidth = 1, colour = "red") +
  labs(x = "Mean expression",
       y = "Variance",
       title = "PBMMC-1: Log Normalized - Mean vs Variance")


#### Selection of highly variable genes ####

# Extract the 1000 most variable genes
# using the variance stabilising transformation method
seurat_pbmmc1 <- FindVariableFeatures(
  seurat_pbmmc1,
  selection.method = "vst",
  nfeatures = 1000
)
VariableFeaturePlot(seurat_pbmmc1)

# Extract the top 10% most variable genes
seurat_pbmmc1 <- FindVariableFeatures(
  seurat_pbmmc1,
  selection.method = "vst",
  nfeatures = round(0.1 * nrow(seurat_pbmmc1))
)
VariableFeaturePlot(seurat_pbmmc1)

# Extract the 3000 most variable genes to match the default in sctransform
seurat_pbmmc1 <- FindVariableFeatures(
  seurat_pbmmc1,
  selection.method = "vst",
  nfeatures = 3000
)
VariableFeaturePlot(seurat_pbmmc1)

# Get the variable genes selected by the vst method
hvgs_vst <- VariableFeatures(seurat_pbmmc1)
head(hvgs_vst)


#### sctransform normalisation ####

# Run sctransform normalisation
# regressing out the percentage of mitochondrial gene expression
seurat_pbmmc1 <- SCTransform(
  seurat_pbmmc1,
  vars.to.regress = "percent.mt",
  verbose = FALSE
)

# A new assay named `SCT` is added to the Seurat object
# 'counts' layer contains the corrected UMI counts
# 'data' layer contains the log-transformed counts data, with a pseudocount of 1 added
# 'scale.data' layer contains the Pearson residuals of the model fit
seurat_pbmmc1

# Get the variable genes selected by sctransform
# which are ranked by residual variance
hvgs_sct <- VariableFeatures(seurat_pbmmc1)
head(hvgs_sct)

# Compare the variable genes selected
# by vst on the raw data and from sctransform residual variance
length(intersect(hvgs_vst, hvgs_sct))


#### Mean-variance relationship after sctransform ####

# Compare the mean-variance relationship for the sctransform-normalised data
# Get summary statistics for the Pearson residuals from the sctransform normalisation
vst_cts <- GetAssayData(seurat_pbmmc1, assay = "SCT", layer = "scale.data")
gene_vst_stats <- tibble(gene = rownames(vst_cts),
                         mean = rowMeans(vst_cts),
                         variance = rowVars(vst_cts))

# Plot the mean-variance relationship for the sctransform-normalised data
ggplot(gene_vst_stats, aes(x = mean, y = variance)) +
  geom_point(alpha = 0.5) +
  labs(x = "Mean expression",
       y = "Variance",
       title = "PBMMC-1: sctransform Normalized - Mean vs Variance")


#### Total UMIs per cell after sctransform ####

# Compare the total UMI counts per cell for the raw counts and sctransform-normalised data

# get the total UMIs for raw UMIs
raw_cts_total <- GetAssayData(seurat_pbmmc1, 
                              assay = "RNA", 
                              layer = "counts") %>%
  colSums() %>%
  enframe(name = "cell", value = "total_counts")

# get the total UMIs for sctransform-normalised counts
sct_cts_total <- GetAssayData(seurat_pbmmc1, 
                              assay = "SCT", 
                              layer = "counts") %>%
  colSums() %>%
  enframe(name = "cell", value = "total_counts")

# Plot the total UMI counts per cell for the PBMMC-1 sample before normalisation
p_raw <- raw_cts_total %>%
  # reorder the cells by their total UMI counts for plotting
  mutate(Cell = fct_reorder(cell, total_counts)) %>%
  # make a bar plot of the total UMI counts per cell
  ggplot(aes(x = Cell, y = total_counts)) +
  geom_col() +
  labs(x = "Cell",
       y = "Total cell UMI counts",
       title = "Before Normalization") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color = "firebrick")
  )

p_sct <- sct_cts_total %>%
  # reorder the cells by their total UMI counts for plotting
  mutate(Cell = fct_reorder(cell, total_counts)) %>%
  # make a bar plot of the total UMI counts per cell
  ggplot(aes(x = Cell, y = total_counts)) +
  geom_col() +
  labs(x = "Cell",
       y = "Total cell UMI counts",
       title = "SCTransform") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color = "firebrick")
  )

# combine the two plots side by side
p_raw + p_sct

#### Exercise: Apply normalisation to the HHD-1 sample ####

# 1. Subset to HHD-1 sample
# YOUR CODE HERE

# 2. Remove genes not expressed in this sample
# YOUR CODE HERE

# 3. Apply log-normalisation
# YOUR CODE HERE

# 4. Select top 3000 variable genes using vst method
# YOUR CODE HERE

# 5. Apply sctransform normalisation
# YOUR CODE HERE

# Compare mean-variance relationships
# YOUR CODE HERE
