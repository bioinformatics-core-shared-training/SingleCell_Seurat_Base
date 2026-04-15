# Load the libraries we will need for this practical
library(Seurat)
library(sctransform)
library(tidyverse)

#### Load data for a single sample ####

# Read 10x data matrices
etv6_runx1_1_data <- Read10X(
  data.dir = "Data/CellRanger_Outputs/SRR9264343/outs/filtered_feature_bc_matrix/"
)

# Create Seurat object from the loaded matrices
etv6_runx1_1 <- CreateSeuratObject(
  counts = etv6_runx1_1_data,
  project = "ETV6_RUNX1_1"
)


# Look at the structure of the Seurat object
etv6_runx1_1

# The metadata can be accessed using [[]]
head(etv6_runx1_1[[]])

# Set the default assay to RNA (raw counts)
DefaultAssay(etv6_runx1_1) <- "RNA"

# Check the active assay was set correctly
DefaultAssay(etv6_runx1_1)

# Pull out the raw counts matrix for the RNA assay for exploratory analysis
raw_counts <- etv6_runx1_1[["RNA"]]$counts

#### Genes detected per cell ####

# Calculate the number of genes detected per cell
genes_per_cell <- colSums(raw_counts > 0)

# Plot the distribution of genes detected per cell
plot(density(genes_per_cell),
     main = "",
     xlab = "Genes per cell")

#### UMIs vs detected cells ####

# Calculate for each average UMI count per cell in which it is expressed
avg_umis_per_gene <- rowSums(raw_counts) / rowSums(raw_counts > 0)

# Calculate the proportion of cells expressing each gene
proportion_cells_expressing <- rowMeans(raw_counts > 0)

# Plot the relationship
plot(
  x = avg_umis_per_gene,
  y = proportion_cells_expressing,
  log = "x",
  xlab = "Expression level per gene",
  ylab = "Proportion of cells expressing the gene"
)

#### Gene count distribution ####

# Calculate the relative expression of each gene in each cell
rel_expression <- t(t(raw_counts) / colSums(raw_counts)) * 100

# Get the 20 most expressed genes
most_expressed <- sort(rowSums(rel_expression), decreasing = TRUE)[20:1]

# Plot the distribution of relative expression for the top 20 genes
plot_data <- as.matrix(t(rel_expression[names(most_expressed), ]))

boxplot(plot_data, cex = 0.1, las = 1,
        xlab = "% total count per cell",
        horizontal = TRUE)

#### Quality control filtering ####

# How many genes are detected in the whole dataset?
table(rowSums(raw_counts) > 0)

# Identify mitochondrial genes
mito_genes <- str_subset(rownames(raw_counts), pattern = "^MT-")
length(mito_genes)

# Calculate the percentage of UMIs mapped to mitochondrial genes 
# and add this as a column to the metadata
etv6_runx1_1[["percent.mt"]] <- PercentageFeatureSet(etv6_runx1_1,
                                                     pattern = "^MT-")

# Check the metadata to see the new column
head(etv6_runx1_1[[]])

# Distribution of number of genes detected per cell
VlnPlot(etv6_runx1_1, features = c("nFeature_RNA"), layer = "counts")

# Distribution of total UMI counts per cell
VlnPlot(etv6_runx1_1, features = c("nCount_RNA"), layer = "counts")

# Distribution of percentage of mitochondrial genes per cell
VlnPlot(etv6_runx1_1, features = c("percent.mt"), layer = "counts")

# Distribution of all three QC metrics together
VlnPlot(etv6_runx1_1,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, layer = "counts")

#### Multiple samples ####

# Read our sample information sheet
sampleinfo <- read_tsv("Data/sample_sheet.tsv")

# Create a named vector of directories for each sample
# We work with a subset of samples for demonstration purposes
samples <- sampleinfo$Sample[c(1, 5, 7, 9)]
list_of_files <- file.path("Data/CellRanger_Outputs/",
                           samples,
                           "/outs/filtered_feature_bc_matrix")
names(list_of_files) <- sampleinfo$SampleName[c(1, 5, 7, 9)]

# Read in the cellranger matrices for all samples
expression_matrix <- Read10X(data.dir = list_of_files)

# Create a new seurat object with the combined expression matrix
multi_seurat_object <- CreateSeuratObject(counts = expression_matrix)
multi_seurat_object

# Pull out the metadata and add sample group and sample name information
temp_metadata <- multi_seurat_object[[]] %>%
  # Add the cell barcodes as a column so we can pull them back in at the end
  rownames_to_column("Cell") %>%
  # Extract the sample group information from the origin identifier
  # The pattern "-.*" matches the first dash and everything after it
  mutate(SampleGroup = str_remove(orig.ident, "-.*")) %>% 
  # Add the sample name as a column
  mutate(SampleName = orig.ident) %>%
  # Put the cell barcodes back as rownames
  column_to_rownames("Cell")

# Add the modified metadata back to the Seurat object
multi_seurat_object[[]] <- temp_metadata

#### QC: remove undetected genes ####

# Get the raw counts to work with
multi_raw_counts <- multi_seurat_object[["RNA"]]$counts

# How many genes are detected in the whole dataset?
table(rowSums(multi_raw_counts) > 0)

# Filter out genes that are not detected in any cell
filtered_multi_seurat_object <- subset(
  multi_seurat_object,
  features = rownames(multi_seurat_object)[rowSums(multi_raw_counts) > 0]
)

# Check how many genes are detected after filtering
filtered_multi_seurat_object

# Confirm all of these now have at least 1 count across all cells
filt_raw_counts <- filtered_multi_seurat_object[["RNA"]]$counts
table(rowSums(filt_raw_counts) > 0)

#### QC: visualising metrics ####

# Add the percentage of mitochondrial genes to the metadata
filtered_multi_seurat_object[["percent.mt"]] <-
  PercentageFeatureSet(filtered_multi_seurat_object, pattern = "^MT-")

# Plot the distribution of QC metrics for the merged samples
VlnPlot(filtered_multi_seurat_object,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  layer = "counts"
)

VlnPlot(
  filtered_multi_seurat_object,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  layer = "counts",
  pt.size = 0
)

#### QC: filtering low-quality droplets ####

# Plot where we would set a threshold
VlnPlot(
  filtered_multi_seurat_object,
  features = "nFeature_RNA",
  group.by = "SampleGroup",
  layer = "counts",
  pt.size = 0
) +
  geom_hline(yintercept = 500, color = "red")

# Pull out the metadata for easier manipulation
metadata <- filtered_multi_seurat_object[[]]

#### Example of how to filter for one sample ####

# Subset the metadata for one of the samples
e_cells <- metadata %>%
  filter(SampleName == "ETV6RUNX1-1")

# calculate the median and MAD for the number of genes detected in this sample
feature_mad <- mad(e_cells$nFeature_RNA)
feature_median <- median(e_cells$nFeature_RNA)

# set a threshold for the number of genes detected as 2 MADs below the median
feature_threshold <- feature_median - 2 * feature_mad

# How many cells would this remove?
length(e_cells$nFeature_RNA) # total cells
sum(e_cells$nFeature_RNA > feature_threshold) # cells above threshold

# The same process can be repeated for the other metrics
umi_threshold <- median(e_cells$nCount_RNA) - 2 * mad(e_cells$nCount_RNA)
sum(e_cells$nCount_RNA > umi_threshold)

# for mitochondrial content we want to set a threshold above the median
mt_threshold <- median(e_cells$percent.mt) + 2 * mad(e_cells$percent.mt)
sum(e_cells$percent.mt < mt_threshold)

# visualise the thresholds on the plot
ggplot(e_cells,
       aes(x = nFeature_RNA,
           y = nCount_RNA,
           color = percent.mt)) +
  geom_point() +
  geom_vline(xintercept = feature_threshold, color = "red") +
  geom_hline(yintercept = umi_threshold, color = "blue") +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(
    title = "QC Metrics for ETV6-RUNX1_1",
    x = "Number of Features (Genes) Detected",
    y = "Number of Counts (UMIs)",
    color = "Percent Mitochondrial"
  )

# We cound apply all 3 filters together to get a vector of cell
# barcodes to keep for this sample
# using standard dplyr filtering and pulling out the cell barcodes
e_keep_cells <- metadata %>%
  rownames_to_column("Cell") %>%
  filter(SampleName == "ETV6RUNX1-1") %>%
  filter(nFeature_RNA > feature_threshold,
         nCount_RNA > umi_threshold,
         percent.mt < mt_threshold) %>%
  pull(Cell)

# check how many cells would be left for this sample
length(e_keep_cells)

#### Apply filtering across all samples ####

all_keep_cells <- metadata %>%
  # Add the cell barcodes as a column so we can pull them back in at the end
  rownames_to_column("cell") %>%
  # Group by sample name to calculate thresholds per sample
  group_by(SampleName) %>%
  # Apply the filtering for each sample
  filter(nFeature_RNA > (median(nFeature_RNA) - 2 * mad(nFeature_RNA)),
         nCount_RNA > (median(nCount_RNA) - 2 * mad(nCount_RNA)),
         percent.mt < (median(percent.mt) + 2 * mad(percent.mt))) %>%
  # Pull out the cell barcodes for the cells that pass the filters
  pull(cell)

# How many cells would be left after filtering?
length(all_keep_cells)

# How many cells per sample
table(metadata[all_keep_cells, "SampleName"])

qc_seurat_object <- subset(
  filtered_multi_seurat_object,
  cells = all_keep_cells
)

VlnPlot(qc_seurat_object,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  layer = "counts",
  pt.size = 0
)


#### Exercises ####

# Apply the same analysis to all 11 samples

# Read in the sample information sheet
sampleinfo <- read_tsv("Data/sample_sheet.tsv")

# Create a named vector of directories for each sample
# YOUR CODE HERE

# Read in the cellranger matrices for all samples
# YOUR CODE HERE

# Create the seurat object
# YOUR CODE HERE



# Sample group and sample name information to the metadata
temp_metadata <- seurat_object[[]] %>%
  rownames_to_column("Cell") %>%
  mutate(SampleGroup = str_remove(orig.ident, "-.*")) %>% 
  mutate(SampleName = orig.ident) %>%
  column_to_rownames("Cell")

seurat_object[[]] <- temp_metadata
head(seurat_object[[]])

# Add percentage of UMIs mapped to mitochondrial genes to the metadata
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,
                                                      pattern = "^MT-")

# Create a vector of colours of each sample group
colours <- c("cyan3", "cyan3", "cyan3", "cyan3",
             "darkgoldenrod1", "darkgoldenrod1",
             "blue", "blue", "blue",
             "lightgreen", "lightgreen")

# Plot the distribution of total UMI counts per cell
VlnPlot(seurat_object, 
        features = c("nCount_RNA"),
        cols = colours,
        layer = "counts",
        pt.size = 0) +
  ggtitle("Total count of UMIs")
  
# Plot the distribution of number of genes detected per cell
# YOUR CODE HERE

# Plot the distribution of percentage of UMIs aligned to mitochondrial genes
# YOUR CODE HERE





# Extract the metadata for filtering
metadata <- seurat_object[[]]

# Apply the filtering at the SampleGroup level
# YOUR CODE HERE

# Number of cells that would be removed
# YOUR CODE HERE
