# 1.3_cluster_fractions_across_regions.r
#
# Purpose:
#   Calculate broad cell-class proportions across adult human cortex regions and
#   visualize regional differences as a grouped bar plot.
#
# Input:
#   ./data/All_region_cortex_metadata.csv
#     Cell-level metadata exported after Seurat integration and annotation.
#     Required columns include MainLoc, location, and class_label.
#
# Output:
#   ./Fig1/Class_prop.pdf
#     Bar plot showing mean class-level proportions across anatomical regions.

# Load data manipulation and plotting libraries.
library(dplyr)
library(ggplot2)

# Read cell-level metadata. Row names are expected to be cell barcodes.
df <- read.csv('./data/All_region_cortex_metadata.csv', row.names = 1)

# Set plot size for notebook or interactive R sessions.
options(repr.plot.width = 8, repr.plot.height = 6)

# Set the anatomical region order used in the final figure.
# This ensures that ggplot displays regions in a biologically interpretable order
# rather than alphabetical order.
df$MainLoc <- factor(
  df$MainLoc,
  levels = c('Frontal lobe', 'Temporal lobe', 'Parietal lobe', 'Occipital lobe', 'Insula lobe')
)

# Define colors for the five cortical regions. The order must match the factor
# levels defined above.
MainLoc_color <- c('#EB7369', '#2AA4DE', '#27B076', '#9C9E23', '#B273AE')

# Define the order of broad cell classes.
cluster_level <- c("Glutamatergic", 'GABAergic', 'Non_neuronal')

# Count cells per region, sample/location, and broad class. Then convert counts
# into within-sample proportions so that samples with different cell numbers are
# comparable.
g <- df %>%
  mutate(class_label = factor(class_label, levels = cluster_level)) %>%
  count(MainLoc, location, class_label) %>%
  group_by(MainLoc, location) %>%
  mutate(freq = n / sum(n)) %>%
  ggplot(aes(x = class_label, y = freq, fill = MainLoc)) +
  # Show mean +/- standard error across samples/locations within each region.
  stat_summary(fun.data = mean_se, geom = "errorbar",
               position = position_dodge(0.75), width = 0.4, color = "black") +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(0.75),
               width = 0.65) +
  scale_fill_manual(values = MainLoc_color) +
  labs(x = "", y = "Cell type proportion in brain regions") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())

# Display the plot in the current graphics device.
plot(g)

# Save the figure as a PDF for manuscript use.
ggsave(filename = "./Fig1/Class_prop.pdf", plot = g, width = 4.5, height = 3)
