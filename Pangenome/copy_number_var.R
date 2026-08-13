# Nichole Ginnan (nginn001@ucr.ed)
# drafted April 2026

# load libs
library(ggVennDiagram)
#######################################################
#######################################################
#### Single Copy Gene Venn Diagram  ###################
#######################################################
#######################################################
gene_df <- read_delim("structure/output_files/gene_presence_absence_binary.tsv")
meta <- "genome_metadata_with_files.tsv"
out_pdf  <- "structure/output_files/gene_presence_absence_heatmap.pdf"
to_drop <- c("KNZ2-5D","KNZ2-6D","KNZ12-10F","KNZ12-1A","KNZ12-7H","KNZ1-1H",
             "TLI6-11F","TLI6-IE","TLI6-1G","TLI6-5H","TLI6-8H","TLI6-4A","TLI7-6A","TLI8-4H",
             "TLI8-9D","TLI8-10B","TLI9-1C","TLI9-11B","TLI9-11C","SVR3-5F")

# Remove clone strains from gene presence/absence matrix
gene_df <- gene_df[, !colnames(gene_df) %in% to_drop]
# Remove clone strains from metadata
meta2 <- meta %>% filter(strain %in% colnames(gene_df),!strain %in% to_drop)

# Helper: get genes present in any strain from a site
genes_by_site <- split(meta2$strain, meta2$site) %>%
  lapply(function(strains) {
    rownames(gene_df)[rowSums(gene_df[, strains, drop = FALSE]) > 0]
  })
# Check sizes
sapply(genes_by_site, length)
ggVennDiagram(genes_by_site[c("SVR", "HAY", "TLI", "KNZ")], label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "grey50") +
  theme_void()
ggsave(filename= "gene.venn.diagram.four.sites.pdf", device="pdf", units="mm", dpi=300, width=100, height=70, path="structure/plots")

# region
genes_by_region <- split(meta2$strain, meta2$region) %>%
  lapply(function(strains) {
    rownames(gene_df)[rowSums(gene_df[, strains, drop = FALSE]) > 0]
  })
# Check sizes
sapply(genes_by_region, length)
ggVennDiagram(genes_by_region[c("WK", "TLI", "KNZ")], label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "grey50") +
  theme_void()
ggsave(filename= "gene.venn.diagram.three.regions.pdf", device="pdf", units="mm", dpi=300, width=100, height=70, path="structure/plots")

# species_site
genes_by_species_site <- split(meta2$strain, meta2$species_site) %>%
  lapply(function(strains) {
    rownames(gene_df)[rowSums(gene_df[, strains, drop = FALSE]) > 0]
  })
# Check sizes
sapply(genes_by_species_site, length)
ggVennDiagram(genes_by_species_site[c("III_SVR", "II_SVR", "IV_HAY", "I_TLI", "I_KNZ")], label_alpha = 0, label = "count",label_size =4) +
  scale_fill_gradient(low = "white", high = "grey50") +
  theme_void()
ggsave(filename= "gene.venn.diagram.three.species_sites.pdf", device="pdf", units="mm", dpi=300, width=100, height=80, path="structure/plots")

##############################################################################
##############################################################################
# Gene copy number analysis ##################################################
##############################################################################
##############################################################################
cnv<-read_delim("structure/output_files/roary_nosplit_CNV/gene_copy_number_matrix.tsv")
to_drop <- c("KNZ2-5D","KNZ2-6D","KNZ12-10F","KNZ12-1A","KNZ12-7H","KNZ1-1H",
             "TLI6-11F","TLI6-IE","TLI6-1G","TLI6-5H","TLI6-8H","TLI6-4A","TLI7-6A","TLI8-4H",
             "TLI8-9D","TLI8-10B","TLI9-1C","TLI9-11B","TLI9-11C","SVR3-5F")
# Remove clone strains from gene presence/absence matrix
cnv <- cnv[, !colnames(cnv) %in% to_drop]
# keep only genes with copy number > 1 in at least one strain
cnv_multi <- cnv %>%filter(if_any(-Gene, ~ .x > 1))
dim(cnv_multi)
write_tsv(cnv_multi, "structure/output_files/gene_copy_number_matrix_multicopy_only.tsv")
# convert to long format
cnv_long <- cnv_multi %>% pivot_longer(cols = -Gene, names_to = "strain", values_to = "copy_number") %>%
  left_join(meta2, by = "strain")
glimpse(cnv_long) # look at summary

# summarize average copy number by site
cnv_site_summary <- cnv_long %>% group_by(Gene, site) %>%
  summarise(
    mean_copy = mean(copy_number, na.rm = TRUE),
    median_copy = median(copy_number, na.rm = TRUE),
    prop_multicopy = mean(copy_number > 1, na.rm = TRUE), # fraction of strains that have a copy number greater than 1
    n = n(), .groups = "drop")
# summarize average copy number by species
cnv_species_summary <- cnv_long %>% group_by(Gene, putative_species) %>%
  summarise(
    mean_copy = mean(copy_number, na.rm = TRUE),
    median_copy = median(copy_number, na.rm = TRUE),
    prop_multicopy = mean(copy_number > 1, na.rm = TRUE),
    n = n(), .groups = "drop")
# summarize average copy number by species_site
cnv_species_site_summary <- cnv_long %>% group_by(Gene, species_site) %>%
  summarise(
    mean_copy = mean(copy_number, na.rm = TRUE),
    median_copy = median(copy_number, na.rm = TRUE),
    prop_multicopy = mean(copy_number > 1, na.rm = TRUE),
    n = n(), .groups = "drop")

# Stats - Kruskal-Wallis per gene - site
site_tests <- cnv_long %>%
  group_by(Gene) %>%
  summarise(
    p_site = tryCatch(
      kruskal.test(copy_number ~ site)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(padj_site = p.adjust(p_site, method = "BH")) %>%
  arrange(padj_site)

# Stats - Kruskal-Wallis per gene - species
species_tests <- cnv_long %>%
  group_by(Gene) %>%
  summarise(
    p_species = tryCatch(
      kruskal.test(copy_number ~ putative_species)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(padj_species = p.adjust(p_species, method = "BH")) %>%
  arrange(padj_species)


# Binary multicopy model: does multicopy frequency differ by site after accounting for species?
cnv_long_filt <- cnv_long %>%
  mutate(multicopy = as.integer(copy_number > 1))

glm_results <- cnv_long_filt %>% group_by(Gene) %>% group_modify(~{
    fit <- tryCatch(glm(multicopy ~ site + putative_species, data = .x, family = binomial),
      error = function(e) NULL)
    if (is.null(fit)) {
      tibble(term = NA_character_, p = NA_real_)
    } else {broom::tidy(fit) %>%
        dplyr::select(term, p.value) %>%
        rename(p = p.value)}}) %>%ungroup()

####### heat maps of mean gene copy number #######
site_mean <- cnv_long %>%
  group_by(Gene, site) %>%
  summarise(mean_copy = mean(copy_number, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = site, values_from = mean_copy) %>%
  as.data.frame()
rownames(site_mean) <- site_mean$Gene
site_mean <- site_mean[, -1]
site_mean <- as.matrix(site_mean)

# species
species_mean <- cnv_long %>%
  group_by(Gene, putative_species) %>%
  summarise(mean_copy = mean(copy_number, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = putative_species, values_from = mean_copy) %>%
  as.data.frame()
rownames(species_mean) <- species_mean$Gene
species_mean <- species_mean[, -1]
species_mean <- as.matrix(species_mean)

# species_site
#only keep genes where at least 50% of at least one group has multple copies
genes_keep <- cnv_long %>%
  group_by(Gene, species_site) %>%
  summarise(
    prop_multicopy = mean(copy_number > 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Gene) %>%
  summarise(
    max_prop_multicopy = max(prop_multicopy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(max_prop_multicopy >= 0.5) %>%
  pull(Gene)
# only keep threshold genes
cnv_long_filt <- cnv_long %>%
  filter(Gene %in% genes_keep)

## all multiple copy genes
species_site_mean_filt <- cnv_long_filt %>%
  group_by(Gene, species_site) %>%
  summarise(mean_copy = mean(copy_number, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = species_site, values_from = mean_copy) %>%
  as.data.frame()
rownames(species_site_mean_filt) <- species_site_mean_filt$Gene
species_site_mean_filt <- species_site_mean_filt[, -1]
species_site_mean_filt <- as.matrix(species_site_mean_filt)

## all multiple copy genes
species_site_mean <- cnv_long %>%
  group_by(Gene, species_site) %>%
  summarise(mean_copy = mean(copy_number, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = species_site, values_from = mean_copy) %>%
  as.data.frame()
rownames(species_site_mean) <- species_site_mean$Gene
species_site_mean <- species_site_mean[, -1]
species_site_mean <- as.matrix(species_site_mean)

# order cols
site_order <- c("SVR", "HAY", "TLI", "KNZ")
site_mean <- site_mean[, site_order]
site_mean <- site_mean[, site_order[site_order %in% colnames(site_mean)]]

species_order <- c("III", "II", "IV", "I")
species_mean <- species_mean[, species_order]
species_mean <- species_mean[, species_order[species_order %in% colnames(species_mean)]]

species_site_order <- c("III_SVR", "II_SVR", "IV_HAY", "I_TLI", "I_KNZ")
species_site_mean <- species_site_mean[, species_site_order]
species_site_mean <- species_site_mean[, species_site_order[species_site_order %in% colnames(species_site_mean)]]

species_site_order <- c("III_SVR", "II_SVR", "IV_HAY", "I_TLI", "I_KNZ")
species_site_mean_filt <- species_site_mean_filt[, species_site_order]
species_site_mean_filt <- species_site_mean_filt[, species_site_order[species_site_order %in% colnames(species_site_mean_filt)]]

##### Site heatmap ####
pheatmap(
  site_mean,
  #scale = "row", 
  cutree_cols =2, 
  cutree_rows = 8,
  clustering_distance_rows = "euclidean", 
  cluster_cols = FALSE, 
  clustering_method = "complete", 
  show_rownames = TRUE, 
  show_colnames = TRUE,
  color= brewer.pal(7,"Greys"),
  border_color = "black",
  cellheight=8, cellwidth = 15,
  #treeheight_col= 10,
  #treeheight_row =20,
  fontsize_row = 8,
  fontsize_col = 11,
  angle_col = 45,
  filename="structure/plots/gene_copy_num_site.heatmap.pdf")

##### species_site_filt heatmap ####
pheatmap(
  species_site_mean_filt,
  scale = "row", 
  cutree_cols =1, 
  cutree_rows = 6,
  clustering_distance_rows = "euclidean", 
  cluster_cols = FALSE, 
  clustering_method = "complete", 
  show_rownames = TRUE, 
  show_colnames = TRUE,
  #color= brewer.pal(7,"Greys"),
  color= brewer.pal(7,"RdGy"), # use with z-scores
  border_color = "black",
  cellheight=10, cellwidth = 18,
  #treeheight_col= 10,
  #treeheight_row =20,
  fontsize_row = 8,
  display_numbers = round(species_site_mean_filt, 1),
  number_color = "white",
  fontsize_col = 11,
  angle_col = 45,
  filename="structure/plots/gene_copy_num_species_site.heatmap_filt_z.scores.pdf"
  #filename="structure/plots/gene_copy_num_species_site.heatmap_filt.pdf"
  )

##### species_site heatmap ####
pheatmap(
  species_mean,
  #scale = "row", 
  cutree_cols =2, 
  cutree_rows = 8,
  clustering_distance_rows = "euclidean", 
  cluster_cols = FALSE, 
  clustering_method = "complete", 
  show_rownames = TRUE, 
  show_colnames = TRUE,
  #color= brewer.pal(7,"Greys"),
  color= brewer.pal(7,"RdGy"), # use with z-scores
  border_color = "black",
  cellheight=10, cellwidth = 18,
  #treeheight_col= 10,
  #treeheight_row =20,
  fontsize_row = 8,
  display_numbers = round(species_mean, 1),
  number_color = "white",
  fontsize_col = 11,
  angle_col = 45,
  filename="structure/plots/gene_copy_num_species.heatmap.z.score.rows.pdf")

# Stats - Kruskal-Wallis per gene - species_site
species_site_tests <- cnv_long %>%
  group_by(Gene) %>%
  summarise(
    p_species_site = tryCatch(
      kruskal.test(copy_number ~ species_site)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(padj_species_site = p.adjust(p_species_site, method = "BH")) %>%
  arrange(padj_species_site)

species_site_tests <- cnv_long_filt %>%
  group_by(Gene) %>%
  summarise(
    p = tryCatch(
      kruskal.test(copy_number ~ species_site)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(padj = p.adjust(p, method = "BH")) %>%
  arrange(padj)
write_tsv(species_site_tests, "structure/R_stats/KW.species_sites_filt_gene_CNV.tsv")

################################################
##### Genome size #############################
###############################################
meta
ggplot(aes(x=species_site,y=(Genome_size_bp/1000), color=Dry.Wet, fill=Dry.Wet), data=meta)+
  geom_boxplot(alpha=0.5, outliers = FALSE)+
  geom_point(position=position_jitter(0.3))+
  xlab("Lineage and Site")+
  ylab("Assembled Genome Size (Kb)")+
  theme_classic() +
  scale_color_manual(values=c("#aa4465","#3d5a80"),limits=c("dry", "wet")) +
  scale_fill_manual(values=c("#aa4465","#3d5a80"),limits=c("dry", "wet")) 
ggsave(filename= "gene.venn.diagram.three.species_sites.pdf", device="pdf", units="mm", dpi=300, width=100, height=80, path="structure/plots")

kw.gen<-kruskal.test(Genome_size_bp~species_site,data=meta)
dunn.test::dunn.test(meta$Genome_size_bp, g=meta$species_site, method="bonferroni")

meta_avg <- meta %>%group_by(species_site) %>%summarise(
    mean_genome_size_kb = mean(Genome_size_bp, na.rm = TRUE)/1000,n = n())

ggplot(aes(x=log10_relative_growth_nacl,y=Genome_size_bp), data=meta)+
  geom_point()+
  geom_smooth(method=lm)+
  theme_classic()
