# Load necessary libraries
library(Seurat)
library(Matrix)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(SeuratDisk)
library(stringr)
library(grid)
library(pheatmap)

# Load custom color palette and define heat colors
cls100 <- readRDS('./data/color.rds')
heat.colors <- colorRampPalette(c("white", "grey70", "black"))(100)

# Function to reorder matrix for pheatmap
reorder_matrix <- function(matrix1, by.rows = TRUE) {
  if (by.rows) {
    conf.order <- order(apply(matrix1, 1, which.max))
    matrix1.reordered <- matrix1[conf.order, ]
  } else {
    conf.order <- order(apply(matrix1, 2, which.max))
    matrix1.reordered <- matrix1[, conf.order]
  }
  return(matrix1.reordered)
}

# Function to compare cluster correlation across species
compare_cl <- function(cl, ref.cl, plot.title = NA, plot.silent = TRUE,
                       heat.colors = heat.colors, row.cl.num = min(length(unique(cl)), length(unique(ref.cl)))) {
  
  conf1 <- sweep(table(cl, ref.cl), 1, rowSums(table(cl, ref.cl)), "/")
  conf2 <- reorder_matrix(conf1)
  
  # Cluster co-occurrence matrix
  cl.prop.cocl <- apply(conf1, 2, function(x) {
    min.prop <- apply(expand.grid(x, x), 1, min)
  })
  
  cl.prop.cocl.total <- rowSums(cl.prop.cocl)
  cl.prop.cocl.m <- matrix(cl.prop.cocl.total, nrow(conf1), nrow(conf1), 
                           dimnames = list(rownames(conf1), rownames(conf1)))
  
  ph1 <- pheatmap(conf2, cutree_rows = row.cl.num, clustering_method = "ward.D2",
                  color = heat.colors, fontsize = 6, main = plot.title, silent = plot.silent)
  
  return(list(conf = conf2, cocl = cl.prop.cocl.m, ph = ph1))
}

# Function to generate correlation plot between species clusters
CorrelationPlotFun <- function(adata, human_cluster, mouse_cluster, species_cluster, integrated_cluster) {
  
  # Extract human and mouse clusters
  pos.h <- grep('^human$', adata@meta.data[[species_cluster]])
  pos.m <- grep('^mouse$', adata@meta.data[[species_cluster]])
  
  human.cluster <- paste0('human_', adata@meta.data[[human_cluster]][pos.h])
  mouse.cluster <- paste0('mouse_', adata@meta.data[[mouse_cluster]][pos.m])
  
  names(human.cluster) <- colnames(adata@assays$RNA@data)[pos.h]
  names(mouse.cluster) <- colnames(adata@assays$RNA@data)[pos.m]
  
  # Combine clusters into dataframe
  cluster_hm_df <- as.data.frame(c(human.cluster, mouse.cluster))
  adata <- AddMetaData(adata, metadata = cluster_hm_df)
  
  # Compare clustering
  ref.cl <- adata$cluster_hm
  cca.cl <- adata@meta.data[[integrated_cluster]]
  names(cca.cl) <- colnames(adata@assays$RNA@data)
  
  cl.conf <- compare_cl(ref.cl, cca.cl)
  cocl <- cl.conf$cocl
  
  # Extract subset of co-occurrence matrix for human and mouse clusters
  compare.species <- unique(adata@meta.data[[species_cluster]])
  cocl.subset <- cocl[grepl(compare.species[1], row.names(cocl)),
                      grepl(compare.species[2], row.names(cocl))]
  
  row.names(cocl.subset) <- sub(paste0(compare.species[1], "_"), "", row.names(cocl.subset))
  colnames(cocl.subset) <- sub(paste0(compare.species[2], "_"), "", colnames(cocl.subset))
  
  # Reorder matrix and plot heatmap
  cocl.subset2 <- reorder_matrix(cocl.subset, by.rows = FALSE)
  pheatmap(cocl.subset2, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 6)
}

# Load and process mouse brain region data
mouse_neuron <- LoadH5Seurat(file = './data/mouse_sub_neuron.h5Seurat')
mouse_neuron <- mouse_neuron %>%
  mutate(donor_label = external_donor_name_label, 
         cluster_label_v1 = str_replace_all(subclass_label, "[- ]", "_")) %>%
  subset(idents = 'FRP') %>%
  subset(idents = 'Glutamatergic')

# Load and process human brain region data
human_neuron <- readRDS(file = './data/All_region_cortex_seurat_v1.rds')
human_exc <- human_neuron %>%
  subset(idents = 'Frontal lobe') %>%
  subset(idents = 'Glutamatergic') %>%
  mutate(cluster_label = str_replace_all(clusters2_v1, ',', '_'),
         donor_label = paste0('H', gsub('-', '_', location))) %>%
  sample_n(size = ncol(mouse_exc))

# Integrate human and mouse data
human_exc$orig.ident <- "human"
mouse_exc$orig.ident <- "mouse"

all.data <- merge(x = human_exc, y = mouse_exc, add.cell.ids = c("human", "mouse"))
all.data.list <- SplitObject(all.data, split.by = "donor_label")

# Normalize and find variable features
all.data.list <- lapply(all.data.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Integration process
features <- SelectIntegrationFeatures(object.list = all.data.list)
all.data.list <- lapply(all.data.list, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

frl_anchors <- FindIntegrationAnchors(object.list = all.data.list, anchor.features = features, reduction = "rpca")
frl_combined <- IntegrateData(anchorset = frl_anchors, k.weight = 20)

# Perform clustering and visualization
frl_combined <- frl_combined %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.5)

DimPlot(frl_combined, reduction = "umap", group.by = "seurat_clusters", cols = cls100)
DimPlot(frl_combined, reduction = "umap", group.by = "subclass_label", cols = cls100)
DimPlot(frl_combined, reduction = "umap", group.by = "clusters2_v1", cols = cls100)

# Create and save correlation plot
p <- CorrelationPlotFun(frl_combined, 'clusters2_v1', 'subclass_label', 'source', 'seurat_clusters')
pdf("./Figures/Human_mouse_cluster_label_corrplot_FL.pdf", width = 8, height = 10)
print(p)
dev.off()

# Proportion comparison plot
df <- frl_combined@meta.data %>%
  mutate(cluster_comb = ifelse(source == 'human', clusters2_v1, subclass_label),
         cluster_comb_v1 = factor(cluster_comb, levels = cluster_level))

g.exc <- df %>%
  count(source, donor_label, cluster_comb_v1) %>%
  group_by(source, donor_label) %>%
  mutate(freq = n / sum(n)) %>%
  ggplot(aes(x = cluster_comb_v1, y = freq, color = source, fill = source)) +
  stat_summary(fun.data = mean_se, position = position_dodge(0.75), geom = "errorbar", width = 0.5) +
  stat_summary(fun.y = mean, position = position_dodge(0.75), geom = "bar", width = 0.65) +
  xlab("") + ylab("Cell type proportion in brain regions") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

plot(g.exc)
ggsave(g.exc, file = "./Figures/Human_mouse_prop_FL.pdf", width = 7.5, height = 5)
