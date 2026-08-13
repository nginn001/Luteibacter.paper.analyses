# Luteibacter genome structure analysis
# Nichole Ginnan (nginn001@ucr.edu)
# April 2, 2026
# Figure 3b in the manuscript

##### Load libs ####
library(readr)
library(tidyverse)
library(vegan)
library(pheatmap)
library(RColorBrewer)
library(writexl)

#### Load Data ####
setwd("/Users/nicholeginnan/Documents/KU/Luteibacter/april.26.genome.analyses/")
d<- read_delim("structure/output_files/mash_dist_matrix.tsv")
meta <- read_delim("genome_metadata_with_files.tsv")

#### PCA to look at Core genome structure from MASH ####
# Convert to matrix
d_mat <- d %>% column_to_rownames(var = "...1") %>%   # use strain names as rownames
  as.matrix()
# Ensure it's numeric
d_mat <- apply(d_mat, 2, as.numeric)
rownames(d_mat) <- colnames(d_mat)   # enforce symmetry naming
d_dist <- as.dist(d_mat) # Convert to dist object
isSymmetric(d_mat)
range(d_mat)
#### Unconstrained ordination ####
pcoa <- cmdscale(d_dist, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa$points)
colnames(pcoa_df) <- c("PCo1", "PCo2") # update column names
pcoa_df$strain <- rownames(pcoa_df)
# Merge metadata
plot_df <- left_join(pcoa_df, meta, by = "strain")
# fix lat, long column names (mistake)
plot_df <- plot_df %>% rename(lat = longitude)
plot_df <- plot_df %>% rename(long = latitude)

#### plot ####
ggplot(plot_df, aes(x = PCo1, y = PCo2, color = region)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(title = "Mash distance PCoA",
       color = "Region")

ggplot(plot_df, aes(x = PCo1, y = PCo2, color = precip_in)) +
  geom_point(size = 5, alpha=0.3, position=position_jitter(width = 0.002, height = 0.002)) +
  scale_color_gradientn(colors = c("#aa4465","#3d5a80"))+
  theme_bw() +
  xlab("Axis 1") + ylab("Axis 2") +
  theme(panel.grid = element_blank())
ggsave(filename= "all.mash.dist.PCoA.pdf", device="pdf", units="mm", dpi=300, width=100, height=70, path="structure/plots")