# Nichole Ginnan, University of California, Riverside (nginn001@ucr.edu)
# Created May 5, 2026
# Last updated May 5, 2026
# Analysis and visualization of predicted Kansas Luteibacter phage genes

# Load Libraries ####
library(ggplot2)
library(tidyverse)
library(ggnewscale)
library(stringr)
library(ape)
library(phytools)
library(ggrepel)


# Load data (phage gene predictions and phage completeness were done using Phaster) ####
setwd("/Users/nicholeginnan/Documents/KU/Luteibacter/phage/clones_removed_5_2026")
phage.df<- read_excel("phaster_phage_prediction_results.xlsx")
md<- read.delim("/Users/nicholeginnan/Documents/KU/Luteibacter/april.26.genome.analyses/genome_metadata_with_files.tsv")
# combine phage results with metadata
phage.df <- left_join(phage.df, md, by = "strain")
# Remove obvious clones
to_drop <- c("KNZ2-5D","KNZ2-6D","KNZ12-10F","KNZ12-1A","KNZ12-7H","KNZ1-1H",
             "TLI6-11F","TLI6-IE","TLI6-1G","TLI6-5H","TLI6-8H","TLI6-4A","TLI7-6A","TLI8-4H",
             "TLI8-9D","TLI8-10B","TLI9-1C","TLI9-11B","TLI9-11C","SVR3-5F")
phage.df.nc <- phage.df %>% filter(!strain %in% to_drop) # 76 strains


# Stacked barplot of phage completeness species_site ####
phage_long <- phage.df.nc %>% pivot_longer(cols = c(phage_1, phage_2), names_to="Phage", values_to="Status") %>% filter(Status != "")
unique(phage_long$species_site)

ggplot(phage_long, aes(x = species_site, fill = Status)) +
  geom_bar(position="fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(limits=c("III_SVR", "II_SVR","IV_HAY", "I_TLI","I_KNZ")) +
  labs(y="Proportion of predicted phage", x="Lineage and Site", fill="Status") +
  scale_fill_manual(values=c("#4E9F3D", "#8E7CC3", "skyblue2","grey20"),limits=c("intact","questionable", "incomplete", "none")) +
  theme_bw()+
  theme(legend.text=element_text(size=10,family="sans"),
        legend.box.spacing = unit(0.0, "pt"),
        axis.title.x= element_text(size=10, family="sans", vjust = 0.4),
        strip.text.x= element_text(size=10),  
        axis.title.y=element_text(size=10, family="sans",vjust = 1), 
        strip.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing=unit(0.5, "lines"), 
        axis.text.x= element_text(colour="black", size=10, family="sans", angle=30, hjust=0.8, vjust=1), 
        axis.text.y= element_text(colour="black", size=10, family="sans"))
ggsave(filename= "phage.status.stacked.bar.plot.species.site.pdf", device="pdf", units="mm", dpi=300, width=110, height=90, path="plots/")

##########################################################
# gene count boxplot ####
ggplot(phage.df.nc, aes(x = Dry.Wet, y = phage_genes, fill = Dry.Wet, color=Dry.Wet)) +
  geom_boxplot(outlier.shape = NA, alpha=0.1) +
  geom_point(position=position_jitter(0.3),alpha=0.7, size=2.4, aes(shape=species_site)) +
  scale_shape_manual(values=c("III_SVR"=16, "II_SVR"=15,"IV_HAY"=2,"I_KNZ"=17, "I_TLI"=1)) +
  scale_x_discrete(limits=c("dry","wet")) +
  scale_fill_manual(values=c("#aa4465","#3d5a80"),limits=c("dry", "wet")) +
  scale_color_manual(values=c("#aa4465","#3d5a80"),limits=c("dry", "wet")) +
  labs(y="Number of predicted phage genes", x="Region") +
  theme_bw() +
  theme(legend.text=element_text(size=10,family="sans"),
        legend.box.spacing = unit(0.0, "pt"),
        axis.title.x= element_text(size=10, family="sans", vjust = 0.4),
        strip.text.x= element_text(size=10),  
        axis.title.y=element_text(size=10, family="sans",vjust = 1), 
        strip.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing=unit(0.5, "lines"), 
        axis.text.x= element_text(colour="black", size=10, family="sans"), 
        axis.text.y= element_text(colour="black", size=10, family="sans"))
ggsave(filename= "phage_gene_count_species_site.pdf", device="pdf", units="mm", dpi=300, width=90, height=80, path="plots/")

hist(phage.df.nc$phage_genes)
kruskal.test(phage.df.nc$phage_genes~phage.df.nc$Dry.Wet)
# Kruskal-Wallis chi-squared = 20, df = 1, p-value = 7.745e-06
##########################################################
# Phage genes and tolerance correlation ####
phage.df.nc.v2<-phage.df.nc %>% filter(!strain == "KNZ2-4C") # PEG outlier

ggplot(phage.df.nc.v2, aes(x = phage_genes, y = log10_relative_growth_peg #,label=strain
                        )) +
  geom_point(position=position_jitter(0.3),alpha=0.7, size=2.4, aes(shape=species_site)) +
  scale_shape_manual(values=c("III_SVR"=16, "II_SVR"=15,"IV_HAY"=2,"I_KNZ"=17, "I_TLI"=1)) +
  geom_smooth(method="lm") +theme_bw() +
  #geom_text_repel() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# no correlation
ggsave(filename= "phage_gene_vs.log10.relative.peg.growth.no.labels.pdf", device="pdf", units="mm", dpi=300, width=125, height=100, path="plots/")

anova(lm(log10_relative_growth_peg~phage_genes, data = phage.df.nc.v2))
#             Df Sum Sq Mean Sq F value Pr(>F)
#phage_genes  1  0.268 0.26788  0.4326 0.5133
#Residuals   58 35.918 0.61927 
anova(lm(log10_relative_growth_nacl~phage_genes, data = phage.df.nc))
#             Df Sum Sq Mean Sq F value Pr(>F)
#phage_genes  1  0.023 0.02313  0.0358 0.8505
#Residuals   74 47.860 0.64676               

# NaCl and PEG tolerance ####
phage.df.nc.v2<-phage.df.nc %>% filter(!strain == "KNZ2-4C")
ggplot(phage.df.nc.v2, aes(x = log10_relative_growth_nacl, y = log10_relative_growth_peg , label=strain
                        )) +
  geom_point(position=position_jitter(0.3),alpha=0.7, size=2.4, aes(shape=species_site)) +
  scale_shape_manual(values=c("III_SVR"=16, "II_SVR"=15,"IV_HAY"=2,"I_KNZ"=17, "I_TLI"=1)) +
  geom_smooth(method="lm") +theme_linedraw() + coord_cartesian(xlim=c(-2,1)) +
  geom_text_repel() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename= "log10_relative_nacl_vs.relative.peg.growth.pdf", device="pdf", units="mm", dpi=300, width=125, height=100, path="plots/")

anova(lm(log10_relative_growth_peg~log10_relative_growth_nacl, data = phage.df.nc.v2))
#                             Df  Sum Sq  Mean Sq  F value    Pr(>F)    
# log10_relative_growth_nacl  1  22.759  22.7588   98.31   4.252e-14 ***
# Residuals                  58  13.427  0.2315                      
############################################################################
########## Phage region and Luteibacter whole genome co-phylogeny ##########
############################################################################
# load tree files
host_tree_file  <- "luteibacter_core.treefile" # built with single-copy core genes
phage_tree_file <- "full.phage.region.treefile"
# output tree
out_pdf <- "plots/Luteibacter_core_v_phage_tanglegram.nc.pdf"
# --- obvious clones to remove ---
to_drop <- c("KNZ2-5D","KNZ2-6D","KNZ12-10F","KNZ12-1A","KNZ12-7H","KNZ1-1H",
             "TLI6-11F","TLI6-IE","TLI6-1G","TLI6-5H","TLI6-8H","TLI6-4A","TLI7-6A","TLI8-4H",
             "TLI8-9D","TLI8-10B","TLI9-1C","TLI9-11B","TLI9-11C","SVR3-5F")
# --- read trees ---
host  <- read.tree(host_tree_file)
phage <- read.tree(phage_tree_file)
# --- clean host names ---
clean_host_id <- function(x) {x %>% str_replace("\\.uniq$", "") %>% str_replace("__[0-9a-f]{12}$", "")}
# --- extract Luteibacter strain ID from phage tip name ---
# Handles:
# KNZ2-5D|1
# RKNZ2-5D|1
# KNZ2-5D_R_|1
# KNZ2-5D.uniq|1
# KNZ2-5D__abcdef123456|1
clean_phage_to_host <- function(x) {x %>%
    str_replace("^R", "") %>%                  # remove leading R for reverse phage
    str_replace("_R_", "") %>%                 # remove internal _R_ if present
    str_replace("\\|.*$", "") %>%              # remove |1 or |2 and anything after
    str_replace("\\.uniq$", "") %>%            # remove .uniq
    str_replace("__[0-9a-f]{12}$", "")         # remove hash suffix
}
host$tip.label <- clean_host_id(host$tip.label)
# map every phage tip to host strain
phage_map <- tibble(phage_tip = phage$tip.label, host_id   = clean_phage_to_host(phage$tip.label))
# join metadata
phage_map <- phage_map %>% left_join(phage.df, by = c("host_id" = "strain"))
# --- remove clone hosts from host tree ---
host <- drop.tip(host, intersect(host$tip.label, to_drop))
# --- remove phage tips belonging to clone hosts ---
phage_drop_clones <- phage_map %>%filter(host_id %in% to_drop) %>%pull(phage_tip)
phage <- drop.tip(phage, intersect(phage$tip.label, phage_drop_clones))
# update phage map after dropping clone-associated phages
phage_map <- phage_map %>% filter(phage_tip %in% phage$tip.label)
# --- keep only phages whose host exists in host tree ---
phage_keep <- phage_map %>% filter(host_id %in% host$tip.label) %>% pull(phage_tip)
phage <- drop.tip(phage, setdiff(phage$tip.label, phage_keep))
phage_map <- phage_map %>% filter(phage_tip %in% phage$tip.label)
# --- association table for cophylo ---
assoc_pairs <- phage_map %>%transmute(host = host_id, parasite = phage_tip) %>% as.matrix()
# sanity checks
cat("Host tips kept:", Ntip(host), "\n") #76
cat("Phage tips kept:", Ntip(phage), "\n") #99
cat("Hosts with phage:", length(unique(assoc_pairs[, "host"])), "\n") #68
cat("Hosts with 0 phage:", length(setdiff(host$tip.label, assoc_pairs[, "host"])), "\n") #8
# --- ladderize trees ---
host  <- ladderize(host)
phage <- ladderize(phage)
# --- get phage tree tip y positions ---
par(xpd = NA)
plot(phage, show.tip.label = FALSE, no.margin = TRUE)
lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)

yy <- lp$yy[seq_len(Ntip(phage))]
names(yy) <- phage$tip.label
# --- order host tree by mean phage position ---
ph_by_host <- split(assoc_pairs[, "parasite"], assoc_pairs[, "host"])
host_ymean <- sapply(ph_by_host, function(tips) {mean(yy[tips], na.rm = TRUE)})
host_ymean <- host_ymean[names(host_ymean) %in% host$tip.label]
host_order <- names(sort(host_ymean))
# rotate host tree toward phage order
#host <- rotateConstr(host, host_order)
# --- order phage tree grouped by host ---
phage_order <- unlist(lapply(host_order, function(h) {
  tips_h <- ph_by_host[[h]]
  if (!is.null(tips_h)) names(sort(yy[tips_h])) else character(0)
}), use.names = FALSE)

phage_order <- intersect(phage_order, phage$tip.label)
#phage <- rotateConstr(phage, phage_order)
# --- link colors by host prefix ---
host_first <- substr(assoc_pairs[, "host"], 1, 1)
link_cols <- setNames(
  c("#3d5a80", "#3d5a80", "#aa4465", "#aa4465"),
  c("K", "T", "S", "H"))[host_first]
link_cols[is.na(link_cols)] <- "gray50"
# Assign phage status to each tip:
phage_map <- phage_map %>%mutate(
    phage_number = case_when(
      str_detect(phage_tip, "\\|1$") ~ "1",
      str_detect(phage_tip, "\\|2$") ~ "2",
      TRUE ~ NA_character_),
    Status = case_when(
      phage_number == "1" ~ phage_1,
      phage_number == "2" ~ phage_2,
      TRUE ~ NA_character_),
    Status = factor(Status,
      levels = c("intact", "questionable", "incomplete", "none")))
# assign status colors
#status_cols <- c(
#  "intact"       = "#4E9F3D",
#  "questionable" = "#8E7CC3",
#  "incomplete"   = "skyblue2",
#  "none"         = "grey20")
status_vec <- setNames(phage_map$Status, phage_map$phage_tip)
# --- plot cophylogeny ---
set.seed(1)
pdf(out_pdf, width = 5, height = 11)
cx <- cophylo(host, phage, assoc = assoc_pairs, rotate = TRUE)
plot(cx, assoc = assoc_pairs, fsize = 0.6, link.lwd = 1.3, link.col = link_cols,
  link.type = "curved", link.lty = "solid")

pp <- get("last_plot.cophylo", envir = .PlotPhyloEnv)
names(pp)
right_pp <- pp$right
phy <- cx$trees[[2]]
xx <- right_pp$xx[seq_len(Ntip(phy))]
yy <- right_pp$yy[seq_len(Ntip(phy))]
names(xx) <- phy$tip.label
names(yy) <- phy$tip.label

tip_status <- status_vec[phy$tip.label]
#tip_cols <- status_cols[as.character(tip_status)]
points(xx, yy, pch = 21, #bg = tip_cols,
       col = "black",
       cex = 0.9)

dev.off()
