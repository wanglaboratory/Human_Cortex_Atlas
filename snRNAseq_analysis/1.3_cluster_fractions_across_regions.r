# Load required libraries
library(dplyr)
library(ggplot2)

# Read data
df <- read.csv('./data/All_region_cortex_metadata.csv', row.names = 1)

# Set plot size
options(repr.plot.width = 8, repr.plot.height = 6)

# Set factor levels for MainLoc
df$MainLoc <- factor(df$MainLoc, levels = c('Frontal lobe', 'Temporal lobe', 'Parietal lobe', 'Occipital lobe', 'Insula lobe'))

# Set color palette for MainLoc
MainLoc_color <- c('#EB7369', '#2AA4DE', '#27B076', '#9C9E23', '#B273AE')

# Define levels for cluster class
cluster_level <- c("Glutamatergic", 'GABAergic', 'Non_neuronal')

# Data transformation and plot creation
g <- df %>% 
  mutate(class_label = factor(class_label, levels = cluster_level)) %>% 
  count(MainLoc, location, class_label) %>% 
  group_by(MainLoc, location) %>% 
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = class_label, y = freq, fill = MainLoc)) +
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

# Plot the figure
plot(g)

# Save the figure to a PDF file
ggsave(filename = "./Fig1/Class_prop.pdf", plot = g, width = 4.5, height = 3)
