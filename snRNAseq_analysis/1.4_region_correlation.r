# Load required libraries
library(Seurat)
library(cowplot)
library(patchwork)
library(SeuratDisk)
library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)
library(ggforce)
library(circlize)
library(corrplot)
library(tidyr)
library(viridis)

# Load Seurat object
ch.Exc.sub <- readRDS('./data/All_region_cortex_seurat_Exc.rds')

# Set default assay and cluster identities
DefaultAssay(ch.Exc.sub) <- 'RNA'
Idents(ch.Exc.sub) <- 'clusters2'

# Subsample the data
ch.Exc.sub.subsampled <- ch.Exc.sub[, sample(colnames(ch.Exc.sub), size = 10000, replace = FALSE)]

# Find markers for subclusters
Idents(ch.Exc.sub.subsampled) <- 'clusters2'
exc.markers <- FindAllMarkers(ch.Exc.sub.subsampled)
exc.markers <- exc.markers %>% arrange(desc(avg_log2FC), .by_group = TRUE) # Ensure sorting is correct
top50 <- exc.markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)

# Initialize variables
subclusters <- unique(ch.Exc.sub$clusters2)
region_cor_list <- list()
region_specific_ctp <- c()

# Loop through each subcluster
for (i in seq_along(subclusters)) {
    subcluster <- subclusters[i]
    print(subcluster)
    
    # Subset data and scale
    subcluster_data <- subset(ch.Exc.sub, idents = subcluster)
    Idents(subcluster_data) <- 'MainLoc'
    
    # Select top 50 genes for the subcluster
    subcluster_genes <- top50 %>% filter(cluster == subcluster) %>% pull(gene)
    subcluster_data <- ScaleData(subcluster_data)
    
    # Calculate average expression and correlation
    ave_subcluster_data <- AverageExpression(subcluster_data, features = subcluster_genes, slot = "scale.data")
    ave_subcluster_data.cor <- cor(ave_subcluster_data$RNA)
    ave_subcluster_data.cor[is.na(ave_subcluster_data.cor)] <- 0
    
    # Reorder correlation matrix using hierarchical clustering
    order.hc2 <- corrMatOrder(ave_subcluster_data.cor, order = "hclust", hclust.method = "ward.D")
    ave_subcluster_data.cor <- ave_subcluster_data.cor[order.hc2, order.hc2]
    
    # Subset the correlation matrix for the five regions
    x <- ave_subcluster_data.cor[rownames(ave_subcluster_data.cor) %in% colnames(ave_subcluster_data$RNA),
                                 colnames(ave_subcluster_data.cor) %in% colnames(ave_subcluster_data$RNA)]
    
    # Define color palette for the heatmap
    col2 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
                               "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))
    
    # Plot correlation matrix using corrplot
    corrplot(x, col = col2(10), method = "color", order = "original", 
             hclust.method = "ward.D", col.lim = c(-1, 1), title = subcluster)
    
    # Store correlation values for each region comparison
    regions <- c('Frontal lobe', 'Temporal lobe', 'Parietal lobe', 'Occipital lobe', 'Insula')
    if (identical(sort(regions), sort(rownames(x)))) {
        region_cor_list <- rbind(region_cor_list, c(
            x['Frontal lobe', 'Temporal lobe'],
            x['Frontal lobe', 'Parietal lobe'],
            x['Frontal lobe', 'Occipital lobe'],
            x['Frontal lobe', 'Insula'],
            x['Temporal lobe', 'Parietal lobe'],
            x['Temporal lobe', 'Occipital lobe'],
            x['Temporal lobe', 'Insula'],
            x['Parietal lobe', 'Occipital lobe'],
            x['Parietal lobe', 'Insula'],
            x['Occipital lobe', 'Insula']
        ))
    } else {
        region_specific_ctp <- append(region_specific_ctp, subcluster)
    }
}

# Convert region correlation list to dataframe
colnames(region_cor_list) <- c('Frontal vs Temporal', 'Frontal vs Parietal', 'Frontal vs Occipital', 
                               'Frontal vs Insula', 'Temporal vs Parietal', 'Temporal vs Occipital', 
                               'Temporal vs Insula', 'Parietal vs Occipital', 'Parietal vs Insula', 
                               'Occipital vs Insula')
rownames(region_cor_list) <- setdiff(subclusters, region_specific_ctp)
region_cor_list_df <- as.data.frame(region_cor_list)

# Reshape the dataframe for plotting
region_cor_list_df$celltype <- rownames(region_cor_list_df)
long_data <- pivot_longer(region_cor_list_df, -celltype, names_to = "cluster", values_to = "Correlation")

# Set factors for plotting
long_data$celltype <- as.factor(long_data$celltype)
long_data$cluster <- as.factor(long_data$cluster)

# Generate heatmap with ggplot
col_use <- viridis(10)[c(1, 2, 5, 6, 7, 8, 10)]
ggplot(long_data, aes(x = celltype, y = cluster, fill = Correlation)) +
    geom_tile(width = 1, height = 1, size = 0.1, color = "black") +
    scale_fill_gradientn(colors = col_use) +
    theme_classic() +
    coord_fixed() +
    RotatedAxis() + 
    theme(legend.position = "right",
          axis.text.x = element_text(size = rel(0.6), angle = 45, hjust = 1),
          axis.text.y = element_text(size = rel(0.6)),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank())
