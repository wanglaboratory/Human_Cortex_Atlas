# Load necessary libraries
library(SingleCellExperiment)
library(Matrix)
library(SeuratDisk)
library(Seurat)
library(dplyr)
library(tidyverse)
library(pheatmap)

# Load integrated Seurat object
ch.integrated.filter <- readRDS('./data/All_region_cortex_seurat_v1.rds')

# Extract data and metadata
dat <- GetAssayData(object = ch.integrated.filter, slot = "data")
metadata <- ch.integrated.filter@meta.data

# Create SingleCellExperiment object
sce_all <- SingleCellExperiment(assays = list(counts = dat), colData = metadata)

# Check structure of the SingleCellExperiment object
str(sce_all)

# Save SingleCellExperiment object
saveRDS(sce_all, file = "./data/All_region_cortex_v1_sce.rds")

# Load necessary scripts
source("./data/MetaNeighbor/metaneighbor.R")
source("./data/MetaNeighbor/1v1_analysis.R")

# Load gene sets
genesets <- readRDS("./data/MetaNeighbor/hgnc_syngo.rds")

# Prepare for MetaNeighbor analysis
sce_all$study_id <- sce_all$MainLoc
classes <- unique(sce_all$class_label)

# Get variable genes for each class
vgs <- lapply(classes, function(class) {
  get_variable_genes(sce_all[, sce_all$class_label == class])
})
names(vgs) <- classes

# Initialize lists to store results
mn_1v1_clust <- vector("list", length(classes))
mn_1v1_subclass <- vector("list", length(classes))

# Compute best hits for each class
for (i in seq_along(classes)) {
  f <- sce_all$class_label == classes[i]
  mn_1v1_clust[[i]] <- compute_best_hits(sce_all[vgs[[i]], f], sce_all$clusters2_v1[f], one_vs_one = TRUE)
  mn_1v1_subclass[[i]] <- compute_best_hits(sce_all[vgs[[i]], f], sce_all$clusters2_v1[f], one_vs_one = TRUE)
}

# Define base for heatmap indexing
base <- c(1, 23, 12, 45, 34)
base_arr <- c(base, base + 1, base + 2, base + 3, base + 4, base + 5, base + 6, base + 7, base + 8, base + 9, base + 10)

# Generate heatmap
dat <- mn_1v1_clust[[1]][base_arr, base_arr]
hmcols <- colorRampPalette(c("white", "blue"))(100)  # Define heatmap colors
p <- pheatmap(dat, useRaster = TRUE, cluster_cols = FALSE,
               color = hmcols, fontsize_row = 5, show_colnames = FALSE,
               cluster_rows = FALSE)

# Save heatmap as PDF
ggsave('./Figures/Exc_MainLoc_Celltype.pdf', plot = p)
