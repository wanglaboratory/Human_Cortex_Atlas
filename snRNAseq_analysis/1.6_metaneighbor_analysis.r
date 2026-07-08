# 1.6_metaneighbor_analysis.r
#
# Purpose:
#   Convert the final annotated Seurat object into a SingleCellExperiment object
#   and perform MetaNeighbor one-vs-one reproducibility analysis across cortical
#   regions. The final heatmap summarizes cross-region cell-type similarity.
#
# Inputs:
#   ./data/All_region_cortex_seurat_v1.rds
#     Final integrated and annotated Seurat object.
#   ./data/MetaNeighbor/metaneighbor.R
#   ./data/MetaNeighbor/1v1_analysis.R
#     MetaNeighbor helper scripts.
#   ./data/MetaNeighbor/hgnc_syngo.rds
#     Optional gene set resource loaded for MetaNeighbor-related analyses.
#
# Outputs:
#   ./data/All_region_cortex_v1_sce.rds
#     SingleCellExperiment version of the Seurat object.
#   ./Figures/Exc_MainLoc_Celltype.pdf
#     Heatmap of MetaNeighbor one-vs-one similarity scores.

# Load libraries for object conversion, single-cell data representation, data
# manipulation, and heatmap plotting.
library(SingleCellExperiment)
library(Matrix)
library(SeuratDisk)
library(Seurat)
library(dplyr)
library(tidyverse)
library(pheatmap)

# Load the final integrated Seurat object.
ch.integrated.filter <- readRDS('./data/All_region_cortex_seurat_v1.rds')

# Extract normalized expression data and metadata from the Seurat object.
# The `data` slot usually contains log-normalized expression values.
dat <- GetAssayData(object = ch.integrated.filter, slot = "data")
metadata <- ch.integrated.filter@meta.data

# Create a SingleCellExperiment object required by MetaNeighbor functions.
sce_all <- SingleCellExperiment(assays = list(counts = dat), colData = metadata)

# Inspect object structure for troubleshooting.
str(sce_all)

# Save the converted object so the conversion step does not need to be repeated.
saveRDS(sce_all, file = "./data/All_region_cortex_v1_sce.rds")

# Load MetaNeighbor helper scripts. These scripts should define functions such as
# get_variable_genes() and compute_best_hits().
source("./data/MetaNeighbor/metaneighbor.R")
source("./data/MetaNeighbor/1v1_analysis.R")

# Load gene sets. This object is not directly used in the lines below but is kept
# for compatibility with extended MetaNeighbor analyses.
genesets <- readRDS("./data/MetaNeighbor/hgnc_syngo.rds")

# Use anatomical region as the study ID so that MetaNeighbor evaluates whether
# cell-type labels are reproducible across regions.
sce_all$study_id <- sce_all$MainLoc

# Analyze each major cell class separately.
classes <- unique(sce_all$class_label)

# Identify variable genes within each major class. Restricting to class-specific
# variable genes improves the sensitivity of one-vs-one comparisons.
vgs <- lapply(classes, function(class) {
  get_variable_genes(sce_all[, sce_all$class_label == class])
})
names(vgs) <- classes

# Initialize result containers for one-vs-one MetaNeighbor scores.
mn_1v1_clust <- vector("list", length(classes))
mn_1v1_subclass <- vector("list", length(classes))

# Compute best-hit matrices for every major class. Here `clusters2_v1` provides
# the subclass label being compared across regions.
for (i in seq_along(classes)) {
  f <- sce_all$class_label == classes[i]
  mn_1v1_clust[[i]] <- compute_best_hits(sce_all[vgs[[i]], f], sce_all$clusters2_v1[f], one_vs_one = TRUE)
  mn_1v1_subclass[[i]] <- compute_best_hits(sce_all[vgs[[i]], f], sce_all$clusters2_v1[f], one_vs_one = TRUE)
}

# Define row/column indices used for plotting a selected block of the heatmap.
# The values are dataset-specific and assume a particular ordering of labels in
# the MetaNeighbor output matrix.
base <- c(1, 23, 12, 45, 34)
base_arr <- c(base, base + 1, base + 2, base + 3, base + 4, base + 5,
              base + 6, base + 7, base + 8, base + 9, base + 10)

# Select the excitatory-neuron result matrix and subset to the desired rows/cols.
dat <- mn_1v1_clust[[1]][base_arr, base_arr]

# Define the heatmap color scale and plot the similarity matrix.
hmcols <- colorRampPalette(c("white", "blue"))(100)
p <- pheatmap(dat, useRaster = TRUE, cluster_cols = FALSE,
              color = hmcols, fontsize_row = 5, show_colnames = FALSE,
              cluster_rows = FALSE)

# Save the heatmap as a PDF.
ggsave('./Figures/Exc_MainLoc_Celltype.pdf', plot = p)
