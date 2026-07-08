# 1.8_GWAS_2.r
#
# Purpose:
#   Summarize CELLECT-MAGMA prioritization outputs for multiple neuropsychiatric,
#   behavioral, immune, and neurological GWAS traits. The script calculates
#   enrichment Z-like statistics, adjusts P values, annotates significance, and
#   visualizes region-by-cell-subclass trait enrichments as a heatmap.
#
# Input:
#   ./data/CELLECT-MAGMA/out/prioritization/MainLoc_clusters2_v1*.txt
#     CELLECT-MAGMA cell-type prioritization output files.
#
# Output:
#   ./Figures/GWAS_heatmap.pdf
#     Heatmap of GWAS enrichment signals across cortical regions and cell types.

# Load data manipulation, factor handling, color, heatmap, and clustering tools.
library(tidyverse)
library(forcats)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(dendsort)

# -----------------------------------------------------------------------------
# Step 1. Load CELLECT-MAGMA prioritization results
# -----------------------------------------------------------------------------

# Find all CELLECT-MAGMA prioritization files matching the region-by-cell-type
# prefix generated from CELLEX specificity scores.
files <- list.files(
  path = "./data/CELLECT-MAGMA/out/prioritization",
  pattern = "^MainLoc_clusters2_v1.*.txt$",
  full.names = TRUE
)

# Read all trait-specific result files and combine them into one data frame. The
# trait name is parsed from each file name and stored in the `Trait` column.
results <- bind_rows(lapply(files, function(file) {
  df <- read.csv(file, sep = '\t')
  df$Trait <- gsub('.cell_type_results.txt', '', strsplit(basename(file), '__')[[1]][2])
  return(df)
}))

# Clean cell-type/cluster names by removing the CELLEX file prefix.
results$Cluster <- gsub('MainLoc_clusters2_v1_human_brain_cells__', '', results$Name)

# -----------------------------------------------------------------------------
# Step 2. Select traits and prepare enrichment statistics
# -----------------------------------------------------------------------------

# Define GWAS traits included in the manuscript heatmap.
traits <- c(
  'Sleep', 'Major_depressive_disorder', 'unipolar_depression', 'Schizophrenia',
  'Education', 'Neuroticism', 'Tourette_syndrome', 'Tiredness', 'ADHD',
  'Autism_disease', 'aggressive_behavior', 'Intelligence', 'Allergic',
  'Alzheimer', 'anxiety_disorder', 'Insomnia', 'drug_dependence', 'panic_disorder'
)

# Collect all region-by-cell-type labels represented in the selected traits.
cls <- results %>%
  filter(Trait %in% traits) %>%
  distinct(Cluster) %>%
  pull()

# Calculate FDR-adjusted P values across selected tests and sort by significance.
Cluster_traits <- results %>%
  filter(Trait %in% traits & Cluster %in% cls) %>%
  mutate(FDR = p.adjust(Coefficient_P_value, method = "fdr")) %>%
  arrange(FDR)

# Create an enrichment matrix. Here `Coefficient / Coefficient_std_error` is used
# as a Z-like statistic to represent the direction and magnitude of enrichment.
mat <- results %>%
  filter(Cluster %in% cls & Trait %in% traits) %>%
  mutate(enrich = Coefficient / Coefficient_std_error) %>%
  select(Cluster, Trait, enrich) %>%
  filter(!grepl("adj_BMI", Trait)) %>%
  pivot_wider(names_from = Trait, values_from = enrich)

rownames(mat) <- mat$Cluster
mat <- mat %>% select(-Cluster)

# Create a matching matrix of FDR values for significance annotation.
mat_p <- Cluster_traits %>%
  filter(Cluster %in% cls & Trait %in% traits) %>%
  mutate(enrich = FDR) %>%
  select(Cluster, Trait, enrich) %>%
  filter(!grepl("adj_BMI", Trait)) %>%
  pivot_wider(names_from = Trait, values_from = enrich)

rownames(mat_p) <- mat_p$Cluster
mat_p <- mat_p %>% select(-Cluster)

# -----------------------------------------------------------------------------
# Step 3. Order rows by cortical region and cell subclass
# -----------------------------------------------------------------------------

# Define anatomical region order.
clusters <- c('Frontal lobe', 'Temporal lobe', 'Parietal lobe', 'Occipital lobe', 'Insula')

# Define cell subclass order used for the manuscript heatmap.
traits_detail <- c(
  'L2_3_IT', 'L3_4_IT', 'L4_5_IT_1', 'L4_5_IT_2', 'L4_5_IT_3',
  'L6_IT_1', 'L6_IT_2', 'L5_6_NP', 'L5_ET', 'L6_CT', 'L6b',
  'PVALB_Chc', 'PVALB', 'SST', 'LAMP5_LHX6', 'LAMP5',
  'LAMP5_RELN', 'ADARB2_KCNG1', 'VIP', 'ODC', 'OPC', 'Ast',
  'MG', 'Vascu'
)

# Subset and reorder rows to show every region-by-subclass combination in a
# consistent order.
mat <- mat[as.vector(outer(clusters, traits_detail, paste, sep = "_")), ]
mat_p <- mat_p[as.vector(outer(clusters, traits_detail, paste, sep = "_")), ]

# -----------------------------------------------------------------------------
# Step 4. Create significance labels and plot heatmap
# -----------------------------------------------------------------------------

# Initialize an empty annotation matrix. Symbols indicate FDR thresholds:
# -     FDR < 0.1
# *     FDR < 0.05
# **    FDR < 0.01
# ***   FDR < 0.001
anno <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
colnames(anno) <- colnames(mat)
rownames(anno) <- rownames(mat)
anno[mat_p < 0.1] <- "-"
anno[mat_p < 0.05] <- "*"
anno[mat_p < 0.01] <- "**"
anno[mat_p < 0.001] <- "***"

# Cluster traits to group GWAS phenotypes with similar enrichment patterns.
hclust_rows <- sort_hclust(hclust(dist(mat), method = "ward.D2"))
hclust_cols <- hclust(dist(t(mat)), method = "ward.D2")

# Generate heatmap. Rows are not clustered because the region/cell-type order is
# manually defined above; columns are clustered to group related traits.
options(repr.plot.width = 10, repr.plot.height = 14)
p <- pheatmap(
  mat,
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

# Save the heatmap as a PDF.
ggsave("./Figures/GWAS_heatmap.pdf", plot = p, width = 8, height = 12)
