library(tidyverse)
library(forcats)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(dendsort)

# Define file path and load results
files <- list.files(path = "./data/CELLECT-MAGMA/out/prioritization", 
                    pattern = "^MainLoc_clusters2_v1.*.txt$", 
                    full.names = TRUE)

# Read and combine data from files into a single data frame
results <- bind_rows(lapply(files, function(file) {
    df <- read.csv(file, sep = '\t')
    df$Trait <- gsub('.cell_type_results.txt', '', strsplit(basename(file), '__')[[1]][2])
    return(df)
}))

# Clean up the Cluster names
results$Cluster <- gsub('MainLoc_clusters2_v1_human_brain_cells__', '', results$Name)

# Define traits of interest
traits <- c('Sleep', 'Major_depressive_disorder', 'unipolar_depression', 'Schizophrenia', 
            'Education', 'Neuroticism', 'Tourette_syndrome', 'Tiredness', 'ADHD', 
            'Autism_disease', 'aggressive_behavior', 'Intelligence', 'Allergic', 
            'Alzheimer', 'anxiety_disorder', 'Insomnia', 'drug_dependence', 'panic_disorder')

# Get unique clusters associated with the traits
cls <- results %>%
    filter(Trait %in% traits) %>%
    distinct(Cluster) %>%
    pull()

# Prepare enrichment data
Cluster_traits <- results %>%
    filter(Trait %in% traits & Cluster %in% cls) %>%
    mutate(FDR = p.adjust(Coefficient_P_value, method = "fdr")) %>%
    arrange(FDR)

# Create enrichment matrix
mat <- results %>%
    filter(Cluster %in% cls & Trait %in% traits) %>%
    mutate(enrich = Coefficient / Coefficient_std_error) %>%
    select(Cluster, Trait, enrich) %>%
    filter(!grepl("adj_BMI", Trait)) %>%
    pivot_wider(names_from = Trait, values_from = enrich)

rownames(mat) <- mat$Cluster
mat <- mat %>% select(-Cluster)

# Create p-value matrix
mat_p <- Cluster_traits %>%
    filter(Cluster %in% cls & Trait %in% traits) %>%
    mutate(enrich = FDR) %>%
    select(Cluster, Trait, enrich) %>%
    filter(!grepl("adj_BMI", Trait)) %>%
    pivot_wider(names_from = Trait, values_from = enrich)

rownames(mat_p) <- mat_p$Cluster
mat_p <- mat_p %>% select(-Cluster)

# Define clusters and traits for subsetting
clusters <- c('Frontal lobe', 'Temporal lobe', 'Parietal lobe', 'Occipital lobe', 'Insula')
traits_detail <- c('L2_3_IT', 'L3_4_IT', 'L4_5_IT_1', 'L4_5_IT_2', 'L4_5_IT_3', 
                   'L6_IT_1', 'L6_IT_2', 'L5_6_NP', 'L5_ET', 'L6_CT', 'L6b',
                   'PVALB_Chc', 'PVALB', 'SST', 'LAMP5_LHX6', 'LAMP5', 
                   'LAMP5_RELN', 'ADARB2_KCNG1', 'VIP', 'ODC', 'OPC', 'Ast', 
                   'MG', 'Vascu')

# Filter enrichment matrices
mat <- mat[as.vector(outer(clusters, traits_detail, paste, sep = "_")), ]
mat_p <- mat_p[as.vector(outer(clusters, traits_detail, paste, sep = "_")), ]

# Create annotation matrix for p-values
anno <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
colnames(anno) <- colnames(mat)
rownames(anno) <- rownames(mat)
anno[mat_p < 0.1] <- "-"
anno[mat_p < 0.05] <- "*"
anno[mat_p < 0.01] <- "**"
anno[mat_p < 0.001] <- "***"

# Clustering
hclust_rows <- sort_hclust(hclust(dist(mat), method = "ward.D2"))
hclust_cols <- hclust(dist(t(mat)), method = "ward.D2")

# Generate heatmap
options(repr.plot.width = 10, repr.plot.height = 14)
p <- pheatmap(mat,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(80),
               breaks = seq(-2, 6, 0.1),
               fontsize_row = 8,
               fontsize_col = 8,
               display_numbers = anno,
               fontsize_number = 8,
               labels_col = gsub("_", " ", colnames(mat)),
               labels_row = rownames(mat),
               cluster_rows = FALSE,
               cluster_cols = hclust_cols,
               angle_col = 45,
               main = "GWAS on human cortex region Clusters"
)

# Save heatmap
ggsave("./Figures/GWAS_heatmap.pdf", plot = p, width = 8, height = 12)
