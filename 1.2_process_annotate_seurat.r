library(batchelor)
library(Seurat)
library(Matrix)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(SeuratDisk)

# Load data
cls100 <- readRDS('./Data/color.rds')
cortex.data <- Read10X(data.dir = "~/PythonScripts/Adult_cortical/Data/adult_human_cortex/AHCaggr24/filtered_feature_bc_matrix/")
cortex <- CreateSeuratObject(counts = cortex.data, project = "adult_cortex", min.cells = 10, min.features = 200)

# Add metadata
metadata <- read.csv("./Data/HC_HVG2000_metadata.csv", row.names = 1)
cortex <- AddMetaData(cortex, metadata = metadata)

# Calculate percent mitochondria and visualize
cortex[["percent.mt"]] <- PercentageFeatureSet(cortex, pattern = "^MT-")
VlnPlot(cortex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Split, normalize, and find variable features
ch.list <- SplitObject(cortex, split.by = 'location')
ch.list <- lapply(ch.list, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})

# Integration
ch.anchors <- FindIntegrationAnchors(object.list = ch.list, dims = 1:30)
ch.integrated <- IntegrateData(anchorset = ch.anchors, dims = 1:30)

# PCA and UMAP
DefaultAssay(ch.integrated) <- "integrated"
ch.integrated <- ScaleData(ch.integrated, verbose = FALSE) %>%
                 RunPCA(npcs = 30, verbose = FALSE) %>%
                 RunUMAP(dims = 1:30, min.dist = 0.3, reduction = "pca")

# Clustering
ch.integrated <- FindNeighbors(ch.integrated, reduction = "pca", dims = 1:30) %>%
                 FindClusters(resolution = 1)
DimPlot(ch.integrated, reduction = "umap", label = TRUE)

# Subcluster EXC and annotation
ch.Exc <- subset(ch.integrated, idents = c('30','17','20','16','25','31','7','12','14','18','28','34','43','35','9','40','37'))
ch.Exc <- RunUMAP(ch.Exc, dims = 1:30, min.dist = 0.3, reduction = "pca") %>%
          FindNeighbors(reduction = "pca", dims = 1:30) %>%
          FindClusters(resolution = 1)
DimPlot(ch.Exc, reduction = "umap", label = TRUE, group.by = 'seurat_clusters', cols = gastrul_colors_37)

# Cluster re-labeling
df.exc <- ch.Exc@meta.data
cluster_map <- list('4' = 'L2_3_IT', '5' = 'L2_3_IT', '24' = 'L2_3_IT', '2' = 'L2_3_IT',
                    '6' = 'L2_3_IT', '1' = 'L2_3_IT', '3' = 'L2_3_IT', '11' = 'L3_4_IT',
                    '9' = 'L3_4_IT', '13' = 'L3_4_IT', '18' = 'L4_5_IT_1', '0' = 'L4_5_IT_1',
                    '16' = 'L4_5_IT_1', '17' = 'L4_5_IT_2', '7' = 'L4_5_IT_3', 
                    '21' = 'L4_5_IT_3', '14' = 'L4_5_IT_3', '8' = 'L4_5_IT_3',
                    '22' = 'L6_IT_1', '25' = 'L6_IT_1', '10' = 'L6_IT_1',
                    '12' = 'L6_IT_2', '20' = 'L5_6_NP', '15' = 'L6_CT', 
                    '23' = 'L6b', '26' = 'L5_ET')

df.exc$clusters2 <- unlist(lapply(df.exc$seurat_clusters, function(x) cluster_map[[as.character(x)]]))
ch.Exc@meta.data$clusters2 <- df.exc$clusters2
DimPlot(ch.Exc, reduction = "umap", group.by = 'clusters2', label = TRUE, cols = gastrul_colors_37)

# Subcluster finding
DefaultAssay(ch.Exc) <- 'integrated'
subcluster_labels <- c('L2_3_IT', 'L3_4_IT', 'L4_5_IT_1', 'L4_5_IT_2', 'L4_5_IT_3',
                       'L6_IT_1', 'L6_IT_2', 'L5_6_NP', 'L6_CT', 'L6b', 'L5_ET')
for (label in subcluster_labels) {
  ch.Exc <- FindSubCluster(ch.Exc, label, graph.name = 'integrated_snn', resolution = 0.5, subcluster.name = paste0(label, ".sub"))
}

# Visualization and saving
options(repr.plot.width = 10, repr.plot.height = 10)
DimPlot(ch.Exc, group.by = 'clusters2.sub10', pt.size = 0.1, repel = TRUE, cols = cls100)
saveRDS(ch.Exc, './data/All_region_cortex_seurat_Exc.rds')

# Similar processes can be applied to IN and NN subclustering and annotation...

# Combine and save final data
clusters2_meta <- bind_rows(
  ch.Exc@meta.data %>% select(clusters2, subclusters = clusters2.sub10, class_label = 'Glutamatergic'),
  ch.In@meta.data %>% select(clusters2, subclusters = clusters2.sub9, class_label = 'GABAergic'),
  ch.Nn@meta.data %>% select(clusters2, subclusters = clusters2.sub5, class_label = 'Non_neuronal')
)

ch.integrated.filter <- ch.integrated.sub.sub[, rownames(clusters2_meta)]
ch.integrated.filter <- AddMetaData(ch.integrated.filter, metadata = clusters2_meta)

# Final visualization
options(repr.plot.width = 20, repr.plot.height = 6)
DimPlot(ch.integrated.filter, group.by = 'clusters2', cols = gastrul_colors_37)
DimPlot(ch.integrated.filter, group.by = 'subclusters')
DimPlot(ch.integrated.filter, group.by = 'class_label')

saveRDS(ch.integrated.filter, './data/All_region_cortex_seurat_v1.rds')
