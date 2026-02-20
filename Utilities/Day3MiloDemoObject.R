# load packages
library(SingleCellExperiment)

set.seed(123)

# Load data
sce <- readRDS("RObjects/Annotated.full.ETV6.PBMMC.sce.rds")

#### T cells ####

# downsample cells in one of the groups
idx_tcells <- which(
  sce$label %in% c("T (c2)", "NK T (c11)") &
  sce$SampleGroup == "ETV6RUNX1"
)
length(idx_tcells)

# keep a few to keep to make it more realistic
idx_tcells <- setdiff(idx_tcells, sample(idx_tcells, 500))
length(idx_tcells)

# swap group labels: downsample from ETV6 and assign them to PBMMC
sce$SampleGroup <- as.character(sce$SampleGroup)
sce$SampleName <- as.character(sce$SampleName)
sce$SampleGroup[idx_tcells] <- "PBMMC"
# balance the samples across the three PBMMC samples
sce$SampleName[sce$SampleGroup == "PBMMC"] <- sample(
  c("PBMMC-1", "PBMMC-2", "PBMMC-3"),
  sum(sce$SampleGroup == "PBMMC"),
  replace = TRUE
)
sce$SampleGroup <- factor(sce$SampleGroup)
sce$SampleName <- factor(sce$SampleName)

#### B cells ####

# balance the samples across the four ETV6 samples for B (c3)
idx_bcells <- which(
  sce$SampleGroup == "ETV6RUNX1" & 
  sce$label == "B (c3)"
)
length(idx_bcells)
sce$SampleGroup <- as.character(sce$SampleGroup)
sce$SampleName <- as.character(sce$SampleName)
sce$SampleName[idx_bcells] <- sample(
  c("ETV6RUNX1-1", "ETV6RUNX1-2", "ETV6RUNX1-3", "ETV6RUNX1-4"),
  length(idx_bcells),
  replace = TRUE
)
sce$SampleGroup <- factor(sce$SampleGroup)
sce$SampleName <- factor(sce$SampleName)

saveRDS(sce, file = "RObjects/MiloDemo.sce.rds")
