# Luteibacter functional annotations enrichment analyses
# Nichole Ginnan (nginn001@ucr.edu)
# April 20, 2026
##### Load libs ####
library(readr)
library(tidyverse)
library(RColorBrewer)
library(writexl)
library(emmeans)
library(lme4)
library(janitor)
library(multcomp)
library(stringr)
library(ggbreak)

#### Load Data ####
setwd("/Users/nicholeginnan/Documents/KU/Luteibacter/april.26.genome.analyses/functional_annotations/")
# 
md<- read_delim("~/Documents/KU/Luteibacter/april.26.genome.analyses/genome_metadata_with_files.tsv")
all<- read_delim("output_files/roary_clusters_with_all_annotations.tsv")
# Remove obvious clones
to_drop <- c("KNZ2-5D","KNZ2-6D","KNZ12-10F","KNZ12-1A","KNZ12-7H","KNZ1-1H",
             "TLI6-11F","TLI6-IE","TLI6-1G","TLI6-5H","TLI6-8H","TLI6-4A","TLI7-6A","TLI8-4H",
             "TLI8-9D","TLI8-10B","TLI9-1C","TLI9-11B","TLI9-11C","SVR3-5F")
md.nc <- md %>% filter(!strain %in% to_drop)
all.nc <- all %>% dplyr::select(-all_of(to_drop))
# How many gene clusters were assigned a KO?
all.nc.ko <- all.nc %>%filter(!is.na(ko), ko != "")
nrow(all.nc.ko) #10852
# Add more annotation confident metrics
all.nc.ko <- all.nc.ko %>%mutate(
  score = as.numeric(score),
  threshold = as.numeric(threshold),
  evalue = as.numeric(evalue),
  score_minus_threshold = score - threshold,
  score_threshold_ratio = score / threshold)
# Implement moderate KO confidence threshold
all.nc.ko_conf.mod <- all.nc.ko %>%filter(trusted_hit == "yes"|(evalue <= 1e-50 & score_threshold_ratio >= 0.75))
nrow(all.nc.ko_conf.mod) #6383
# convert to long format
genome_cols <- colnames(all.nc)[29:104] # Genomes cols
all_long <- all.nc.ko_conf.mod %>% pivot_longer(cols = all_of(genome_cols), names_to = "sample", values_to = "present") %>%
  filter(present == 1)

meta_cols <- c(
  "sample","species_site","site","putative_species","region","Dry.Wet",
  "precip_in","longitude","latitude","fasta","gff",
  "relative_growth_nacl","relative_growth_peg",
  "log10_relative_growth_nacl","log10_relative_growth_peg",
  "has_fasta","has_gff","has_genome_files","total_clusters")

############################################################
##### BRITE B enrichment analysis on mod conf filtered KOs #
############################################################
#Separate BRITE A, B, and C into individual IDs to names mapping files
brite_parsed <- read_tsv("mapping_files/brite_ko00001.parsed.tsv")
brite_B_lookup <- brite_parsed %>%dplyr::select(brite_B_id, brite_B_name) %>%
  filter(!is.na(brite_B_id), brite_B_id != "") %>%distinct() %>%arrange(brite_B_id)
brite_B_lookup <- brite_B_lookup %>%mutate(brite_B_id = sub("^0+", "", brite_B_id))
# remove BRITE B categories not related to bacteria/phage
remove_brite_B <- c(
  "Cancer: overview","Cardiovascular disease","Cellular community - eukaryotes",
  "Development and regeneration","Digestive system","Drug resistance: antineoplastic",
  "Endocrine and metabolic disease","Endocrine system","Excretory system",
  "Immune system","Infectious disease: parasitic","Nervous system",
  "Neurodegenerative disease","Sensory system","Signaling molecules and interaction")
# filter, long-format, and count BRITE B IDs per strain
all_brite_long <- all.nc.ko_conf.mod %>% filter(!is.na(brite_B_name),brite_B_name != "",
    !brite_B_name %in% remove_brite_B,!is.na(brite_B_id),brite_B_id != "") %>%
  pivot_longer(cols = all_of(intersect(md.nc$strain, colnames(.))),names_to = "sample",values_to = "present") %>%
  filter(present == 1)
counts_wide.brite_B_id <- all_brite_long %>%count(sample, brite_B_id, name = "n_clusters") %>%
  pivot_wider(names_from = brite_B_id,values_from = n_clusters,values_fill = 0) %>%
  left_join(md.nc, by = c("sample" = "strain"))
counts_wide.brite_B_id %>% count(Dry.Wet)
# feature columns
feature_cols.brite_B_id <- setdiff(colnames(counts_wide.brite_B_id), meta_cols)
# normalize by total BRITE-B annotated cluster counts
counts_norm.brite_B_id <- counts_wide.brite_B_id %>%mutate(
    total_clusters = rowSums(across(all_of(feature_cols.brite_B_id), as.numeric), na.rm = TRUE)) %>%
  mutate(across(all_of(feature_cols.brite_B_id),~ as.numeric(.x) / total_clusters,.names = "{.col}_norm"))
# KW test loop
results_clean.brite_B_id <- map_dfr(feature_cols.brite_B_id, function(i) {
  dat <- counts_norm.brite_B_id %>%
    dplyr::select(Dry.Wet, raw_count = all_of(i), feature = all_of(paste0(i, "_norm"))) %>%
    filter(!is.na(Dry.Wet), !is.na(feature)) %>% mutate(Dry.Wet = factor(Dry.Wet))
  if (nrow(dat) == 0 || n_distinct(dat$feature) < 2) {
    return(NULL)}
  kt <- kruskal.test(feature ~ Dry.Wet, data = dat)
  summary_wide <- dat %>%
    group_by(Dry.Wet) %>% summarise(
      mean_raw_count = mean(raw_count, na.rm = TRUE),
      median_raw_count = median(raw_count, na.rm = TRUE),
      mean_norm = mean(feature, na.rm = TRUE),
      median_norm = median(feature, na.rm = TRUE),
      percent_present = mean(raw_count > 0, na.rm = TRUE) * 100,
      n_strains = n(),.groups = "drop") %>%
    mutate(Dry.Wet = as.character(Dry.Wet)) %>%pivot_wider(
      names_from = Dry.Wet,
      values_from = c(mean_raw_count, median_raw_count, mean_norm,
                      median_norm, percent_present, n_strains),
      names_glue = "{.value}_{Dry.Wet}")
  tibble(test_annotation = i,chi.squared = as.numeric(kt$statistic),df = as.numeric(kt$parameter),p_value = kt$p.value) %>%
    bind_cols(summary_wide)}) %>% mutate(
    p_adj_fdr = p.adjust(p_value, method = "fdr"),
    mean_norm_diff_dry_minus_wet = mean_norm_dry - mean_norm_wet,
    mean_raw_diff_dry_minus_wet = mean_raw_count_dry - mean_raw_count_wet,
    prevalence_diff_dry_minus_wet = percent_present_dry - percent_present_wet) %>%
  left_join(brite_B_lookup, by = c("test_annotation" = "brite_B_id")) %>%arrange(p_adj_fdr)
# final filter
results_filtered.brite_B_id <- results_clean.brite_B_id %>%filter(p_adj_fdr <= 0.05,percent_present_dry >= 75 | percent_present_wet >= 75)
# write table
write_xlsx(list(all_brite_B_id.results = results_clean.brite_B_id,filtered_brite_B_id.results = results_filtered.brite_B_id),
  path = "results/kw.results.brite_B_id.bact.phage.only.mod.conf.filt.v2.xlsx")

############################################################
##### Module enrichment analysis on mod conf filtered KOs #
############################################################
# module ID to module name lookup
module_lookup <- read_tsv("mapping_files/list_module.2026_04_15.tsv",col_names = c("module_id", "module_name"))
# split module IDs, convert to long format, and count module IDs per strain
all_module_long <- all.nc.ko_conf.mod %>%filter(!is.na(module_ids), module_ids != "") %>%
  separate_rows(module_ids, sep = ";\\s*") %>%filter(!is.na(module_ids), module_ids != "") %>%
  pivot_longer(cols = all_of(intersect(md.nc$strain, colnames(.))),names_to = "sample",values_to = "present") %>%
  filter(present == 1)
counts_wide.module_ids <- all_module_long %>%count(sample, module_ids, name = "n_clusters") %>%
  pivot_wider(names_from = module_ids,values_from = n_clusters,values_fill = 0) %>%
  left_join(md.nc, by = c("sample" = "strain"))
counts_wide.module_ids %>% count(Dry.Wet)
# feature columns
feature_cols.module_ids <- setdiff(colnames(counts_wide.module_ids), meta_cols)
# normalize by total module-annotated cluster counts
counts_norm.module_ids <- counts_wide.module_ids %>%
  mutate(total_clusters = rowSums(across(all_of(feature_cols.module_ids), as.numeric),na.rm = TRUE)) %>%
  mutate(across(all_of(feature_cols.module_ids),~ as.numeric(.x) / total_clusters,.names = "{.col}_norm"))
# KW test loop
results_clean.module_ids <- map_dfr(feature_cols.module_ids, function(i) {
  dat <- counts_norm.module_ids %>%dplyr::select(Dry.Wet,raw_count = all_of(i),feature = all_of(paste0(i, "_norm"))) %>%
    filter(!is.na(Dry.Wet), !is.na(feature)) %>% mutate(Dry.Wet = factor(Dry.Wet))
  if (nrow(dat) == 0 || n_distinct(dat$feature) < 2) {
    return(NULL)}
  kt <- kruskal.test(feature ~ Dry.Wet, data = dat)
  summary_wide <- dat %>%group_by(Dry.Wet) %>%summarise(
      mean_raw_count = mean(raw_count, na.rm = TRUE),
      median_raw_count = median(raw_count, na.rm = TRUE),
      mean_norm = mean(feature, na.rm = TRUE),
      median_norm = median(feature, na.rm = TRUE),
      percent_present = mean(raw_count > 0, na.rm = TRUE) * 100,
      n_strains = n(),.groups = "drop") %>% mutate(Dry.Wet = as.character(Dry.Wet)) %>%
    pivot_wider(names_from = Dry.Wet,values_from = c(mean_raw_count, median_raw_count,mean_norm, median_norm,
                                                     percent_present, n_strains),names_glue = "{.value}_{Dry.Wet}")
  tibble(test_annotation = i,chi.squared = as.numeric(kt$statistic),df = as.numeric(kt$parameter),p_value = kt$p.value) %>%
    bind_cols(summary_wide)}) %>% mutate(p_adj_fdr = p.adjust(p_value, method = "fdr"),
    mean_norm_diff_dry_minus_wet = mean_norm_dry - mean_norm_wet,
    mean_raw_diff_dry_minus_wet = mean_raw_count_dry - mean_raw_count_wet,
    prevalence_diff_dry_minus_wet = percent_present_dry - percent_present_wet) %>%
  left_join(module_lookup, by = c("test_annotation" = "module_id")) %>% arrange(p_adj_fdr)
# final filter
results_filtered.module_ids <- results_clean.module_ids %>% filter(p_adj_fdr <= 0.05,percent_present_dry >= 75 | percent_present_wet >= 75)
# write table
write_xlsx(list(all_module_id.results = results_clean.module_ids,filtered_module_id.results = results_filtered.module_ids),
  path = "results/kw.results.module_ids.mod.conf.filt.v2.xlsx")

############################################################
##### KO enrichment analysis on mod conf filtered KOs #
############################################################
# KO ID to KO name lookup
ko_lookup <- read_tsv("mapping_files/list_ko.2026_04_15.tsv",col_names = c("ko_id", "ko_name"))
# convert KO-annotated clusters to long format and count KO IDs per strain
all_ko_long <- all.nc.ko_conf.mod %>%filter(!is.na(ko), ko != "") %>%pivot_longer(
    cols = all_of(intersect(md.nc$strain, colnames(.))),names_to = "sample",values_to = "present") %>%
  filter(present == 1)
counts_wide.ko <- all_ko_long %>%
  count(sample, ko, name = "n_clusters") %>%pivot_wider(
    names_from = ko,
    values_from = n_clusters,
    values_fill = 0) %>%
  left_join(md.nc, by = c("sample" = "strain"))
counts_wide.ko %>% count(Dry.Wet)
# feature columns
feature_cols.ko <- setdiff(colnames(counts_wide.ko), meta_cols)
# normalize by total KO-annotated cluster counts
counts_norm.ko <- counts_wide.ko %>%mutate(total_clusters = rowSums(across(all_of(feature_cols.ko), as.numeric),na.rm = TRUE)) %>%
  mutate(across(all_of(feature_cols.ko),~ as.numeric(.x) / total_clusters,.names = "{.col}_norm"))
# KW test loop
results_clean.ko <- map_dfr(feature_cols.ko, function(i) {
  dat <- counts_norm.ko %>%dplyr::select(Dry.Wet,raw_count = all_of(i),feature = all_of(paste0(i, "_norm"))) %>%
    filter(!is.na(Dry.Wet), !is.na(feature)) %>%mutate(Dry.Wet = factor(Dry.Wet))
  if (nrow(dat) == 0 || n_distinct(dat$feature) < 2) {
    return(NULL)}
  kt <- kruskal.test(feature ~ Dry.Wet, data = dat)
  summary_wide <- dat %>%group_by(Dry.Wet) %>%summarise(
      mean_raw_count = mean(raw_count, na.rm = TRUE),
      median_raw_count = median(raw_count, na.rm = TRUE),
      mean_norm = mean(feature, na.rm = TRUE),
      median_norm = median(feature, na.rm = TRUE),
      percent_present = mean(raw_count > 0, na.rm = TRUE) * 100,
      n_strains = n(),.groups = "drop") %>%mutate(Dry.Wet = as.character(Dry.Wet)) %>%pivot_wider(
      names_from = Dry.Wet,values_from = c(mean_raw_count, median_raw_count,mean_norm, median_norm,percent_present, n_strains),
      names_glue = "{.value}_{Dry.Wet}")
  tibble(test_annotation = i,chi.squared = as.numeric(kt$statistic),df = as.numeric(kt$parameter),p_value = kt$p.value) %>%
    bind_cols(summary_wide)}) %>%mutate(
    p_adj_fdr = p.adjust(p_value, method = "fdr"),
    mean_norm_diff_dry_minus_wet = mean_norm_dry - mean_norm_wet,
    mean_raw_diff_dry_minus_wet = mean_raw_count_dry - mean_raw_count_wet,
    prevalence_diff_dry_minus_wet = percent_present_dry - percent_present_wet) %>%
  left_join(ko_lookup, by = c("test_annotation" = "ko_id")) %>%arrange(p_adj_fdr)
# final filter
results_filtered.ko <- results_clean.ko %>%filter(p_adj_fdr <= 0.05,percent_present_dry >= 75 | percent_present_wet >= 75)
# write table
write_xlsx(list(all_ko.results = results_clean.ko,filtered_ko.results = results_filtered.ko),path = "results/kw.results.ko.mod.conf.filt.v2.xlsx")

############################################################
############# Sorting and Plotting Results #################
############################################################
############## Figure 4a ###################################
# BRITE B bubble significantly different number of gene clusters in a category enriched ####
df <- read_excel("results/kw.results.brite_B_id.bact.phage.only.mod.conf.filt.v2.xlsx",sheet = "filtered_brite_B_id.results")
# choose a single prevalence metric
df <- df %>% mutate(prevalence = pmax(percent_present_dry, percent_present_wet),
  direction = ifelse(mean_norm_diff_dry_minus_wet > 0, "Arid-enriched", "Humid-enriched"),
  neg_log10_p = -log10(p_adj_fdr),log2FC_dry_vs_wet = log2((mean_norm_dry) / (mean_norm_wet)))

#df <- df %>% arrange(log2FC_dry_vs_wet) %>% mutate(brite_B_name = factor(brite_B_name, levels = brite_B_name))
df <- df %>% arrange(mean_norm_diff_dry_minus_wet) %>% mutate(brite_B_name = factor(brite_B_name, levels = brite_B_name))
# plot #
ggplot(df, aes(x = mean_norm_diff_dry_minus_wet, y = brite_B_name, color=direction)) +
  geom_point(size=4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values=c("#aa4465", "#3d5a80"),limits=c("Arid-enriched", "Humid-enriched")) +
  labs(
    #x = "log2 fold-change (Dry / Wet)",
    x= "Effect size: Normalized abundance difference",
    y = "BRITE B Functional Category",
    color = "Enriched") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),axis.text.x = element_text(size = 8, angle = 20, hjust = 0.8),
        panel.grid.minor = element_blank())
ggsave(filename= "brite.b.bact.phage.only_enrichment_bubbleplot.effect.size.mod.conf.filt.pdf", device="pdf", units="mm", dpi=300, width=180, height=125, path="plots")

####################################################################
## What KOs are significantly enriched in those BRITE B? ###########
####################################################################
# 3. BRITE B -> KO map, keeping both module-linked and no-module KOs
brite_ko_map_all <- all.nc.ko_conf.mod %>%filter(brite_B_name %in% sig_brite_B,!is.na(ko), ko != "") %>%
  mutate(module_ids = ifelse(is.na(module_ids) | module_ids == "", "no_module", module_ids)) %>%
  separate_rows(module_ids, sep = ";\\s*") %>%rename(module_id = module_ids) %>%
  distinct(source_brite_subset = brite_B_name,module_id,ko)
# module IDs that passed your significant module filter
sig_module_ids <- sig_modules_from_sig_brite %>% distinct(module_id = test_annotation)
# keep KOs that are either in significant modules OR have no module
brite_ko_map_keep <- brite_ko_map_all %>% filter(module_id == "no_module" | module_id %in% sig_module_ids$module_id)
# 4. Pull KO stats for significant KOs in significant BRITE B categories
sig_kos_from_sig_modules <- results_filtered.ko %>%left_join(brite_ko_map_keep,by = c("test_annotation" = "ko"),
    relationship = "many-to-many") %>%filter(!is.na(source_brite_subset)) %>%
  mutate(module_name = ifelse(module_id == "no_module","No significant KEGG module",NA_character_)) %>%
  left_join(results_filtered.module_ids %>%dplyr::select(module_id = test_annotation, module_name_lookup = module_name),by = "module_id") %>%
  mutate(module_name = ifelse(module_id == "no_module",module_name,module_name_lookup)) %>%
  dplyr::select(-module_name_lookup) %>%relocate(source_brite_subset, module_id, module_name, .before = everything()) %>%
  arrange(source_brite_subset, module_id, p_adj_fdr)
# 5. Collapse module table as before
sig_modules_collapsed <- sig_modules_from_sig_brite %>%
  group_by(test_annotation) %>% summarise(source_brite_subset = paste(sort(unique(source_brite_subset)), collapse = "; "),
    across(everything(), first),.groups = "drop") %>%relocate(source_brite_subset, .before = everything())
# high-level categories for modules
modules_collapsed <- sig_modules_collapsed %>%mutate(high_level_category = case_when(
    str_detect(module_name, regex(
      "betaine|proline|ectoine|trehalose|glutathione|polyamine|sulfate|sulfo|sulfur|DMSP|formaldehyde|glycine cleavage|osm|stress",
      ignore_case = TRUE)) ~ "Stress tolerance & osmoprotection",
    str_detect(module_name, regex(
      "glycan|lipid A|KDO2|UDP|dTDP|CDP-|CMP-KDO|cell wall|membrane|peptidoglycan|phosphatidyl|fatty acid|resistance|nodulation|pathogenicity|CAMP|porin|LPS|EPS",
      ignore_case = TRUE)) ~ "Cell envelope & interactions",
    str_detect(module_name, regex(
      "oxidase|dehydrogenase|electron|respiration|quinone|methane|sulfide|cytochrome|heme|ATPase|succinate dehydrogenase|formate dehydrogenase",
      ignore_case = TRUE)) ~ "Respiration & redox metabolism",
    str_detect(module_name, regex(
      "glycolysis|TCA|citrate cycle|pentose|entner|glyoxylate|gluconeogenesis|pyruvate|acetyl-CoA|beta-oxidation|Wood-Ljungdahl|methylcitrate|methylaspartate",
      ignore_case = TRUE)) ~ "Core carbon & energy metabolism",
    str_detect(module_name, regex(
      "biosynthesis|nucleotide|amino acid|coenzyme|vitamin|purine|pyrimidine|lysine|tryptophan|thiamine|NAD|cobalamin|biotin|riboflavin|folate|pantothenate",
      ignore_case = TRUE)) ~ "Biosynthesis & growth",
    TRUE ~ "Specialized / secondary metabolism"))
# 5b. Collapse KO table, now including no_module KOs
sig_kos_collapsed <- sig_kos_from_sig_modules %>%group_by(test_annotation) %>%summarise(
    source_brite_subset = paste(sort(unique(source_brite_subset)), collapse = "; "),
    module_id = paste(sort(unique(module_id)), collapse = "; "),
    module_name = paste(sort(unique(module_name)), collapse = "; "),
    across(everything(), first),
    .groups = "drop") %>%
  mutate(high_level_category = ifelse(str_detect(module_id, "no_module"),"No significant KEGG module",NA_character_)) %>%
  relocate(source_brite_subset, module_id, module_name, high_level_category, .before = everything())
# 6. Direct BRITE B -> significant KO table, regardless of module assignment
sig_kos_direct_from_sig_brite <- results_filtered.ko %>%left_join(brite_ko_map_all %>%distinct(source_brite_subset, ko),
    by = c("test_annotation" = "ko"),relationship = "many-to-many") %>%
  filter(!is.na(source_brite_subset)) %>%
  relocate(source_brite_subset, .before = everything()) %>%
  arrange(source_brite_subset, p_adj_fdr)
# 7. Collapsed version: one row per KO, with all BRITE B categories combined
sig_kos_direct_collapsed <- sig_kos_direct_from_sig_brite %>%group_by(test_annotation) %>%
  summarise(source_brite_subset = paste(sort(unique(source_brite_subset)), collapse = "; "),across(everything(), first),.groups = "drop") %>%
  relocate(source_brite_subset, .before = everything())

write_xlsx(list(
    modules_long = sig_modules_from_sig_brite,
    modules_collapsed = modules_collapsed,
    kos_long_module_filtered = sig_kos_from_sig_modules,
    kos_collapsed_module_filtered = sig_kos_collapsed,
    kos_long_direct_brite = sig_kos_direct_from_sig_brite,
    kos_collapsed_direct_brite = sig_kos_direct_collapsed),
  "results/significant_modules_and_kos_with_brite_sources_no.mod.too.5_1_26.xlsx")

#####################################################################################
### Identify sig KOs within sig BRITEs that are related to stress and competition ##
#####################################################################################
# Sort KOs into groups, order matters for this sorting proccess
sig_kos_direct_collapsed <- read_excel("results/significant_modules_and_kos_with_brite_sources_no.mod.too.5_1_26.xlsx",sheet = "kos_collapsed_direct_brite")

kos_collapsed_cat <- sig_kos_direct_collapsed %>%mutate(ko_high_level_category = case_when(
      # 1. Phage related
      str_detect(ko_name, regex(
        "phage|prophage|capsid|tail|tape measure|integrase|recombinase|portal|terminase|Straboviridae",
        ignore_case = TRUE)) ~"Phage related",
      # 2. Secretion systems
      str_detect(ko_name, regex(
        "secretion system|type II secretion|type III secretion|type IV secretion|type VI secretion|T2SS|T3SS|T4SS|T6SS|general secretion pathway|\\bgsp[A-Z]?\\b|\\bhcp\\b|\\bvgrG\\b|\\bimpA\\b|\\bimpB\\b|\\bimpC\\b|\\bimpF\\b|\\bimpG\\b|\\bimpH\\b|\\bimpI\\b|\\bimpJ\\b|\\bimpK\\b|\\bimpL\\b|\\bimpM\\b|\\bimpN\\b|\\bvasA\\b|\\bvasB\\b|\\bvasC\\b|\\bvasE\\b|\\bvasF\\b|\\bvasK\\b|\\bicmF\\b|\\bdotU\\b|\\bsecA\\b|\\bsecB\\b|\\bsecD\\b|\\bsecE\\b|\\bsecF\\b|\\bsecG\\b|\\bsecY\\b|\\btatA\\b|\\btatB\\b|\\btatC\\b|\\btolC\\b|\\bhasF\\b|\\bprtF\\b|hemolysin activation/secretion",
        ignore_case = TRUE)) ~"Secretion systems",
      # 3. Osmotic homeostasis / water availability
      str_detect(ko_name, regex(
        "betaine|trehalose|proline|ectoine|osm|aquaporin|mechanosensitive|mscL|mscS|mscK|K\\+|potassium|kdp|kef|antiporter|Na\\+|H\\+ antiporter|chloride channel|glycerol uptake|\\bglp\\b|\\btrk\\b|\\bKdp\\b|\\bKup\\b|\\bopu\\b",
        ignore_case = TRUE)) ~"Osmotic homeostasis",
      # 4. Oxidative stress
      str_detect(ko_name, regex(
        "catalase|peroxidase|peroxiredoxin|thioredoxin|glutaredoxin|glutathione|superoxide|methionine sulfoxide|msrA|msrB|msrP|Dps|DNA protection|oxidative|redox",
        ignore_case = TRUE)) ~"Oxidative stress",
      # 5. Cytochrome / electron transport
      str_detect(ko_name, regex(
        "cytochrome|cytochrome c|cytochrome bd|cytochrome o|cytochrome c oxidase|ubiquinol oxidase|quinol oxidase|coxA|coxB|coxC|coxD|ctaA|ctaB|ctaC|ctaD|ctaE|ctaF|cyoA|cyoB|cyoC|cyoD|cydA|cydB|cydX|petB|CYTB|ccmA|ccmB|ccmC|ccmD|ccmF|ccmG",
        ignore_case = TRUE)) ~"Cytochrome / electron transport",
      # 6. Nitrogen metabolism
      str_detect(ko_name, regex(
        "ammonium|ammonia|nitrogen|nitrate|nitrite|glnA|glnG|ntrC|glnL|ntrB|glnK|glnD|ntrY|nitric oxide|norR|hmp|urease|fixJ|fixL|nitrogen fixation",
         ignore_case = TRUE)) ~"Nitrogen metabolism",
      # 7. Metal / toxic ion stress
      str_detect(ko_name, regex(
        "arsen|chromate|cobalt|zinc|cadmium|mercur|copper|heavy metal|czc|cus|cnr|chrA|arsC|arsR|merR|nickel",
        ignore_case = TRUE)) ~"Metal / toxic ion stress",
      # 8. Resource competition, especially iron / phosphorus / Fe-S
      str_detect(ko_name, regex(
        "iron|siderophore|enterobactin|heme transport|Fe\\(3\\+\\)|Fe3|ferrous|ferric|fec|feo|fiu|afu|fbp|hmu|phu|bhu|TonB|ExbB|ExbD|phosphate transport|phosphonate transport|pst|phnD|molybdate|modA|modB|modC|modE|sufA|sufB|sufC|sufD|sufE|sufS|Fe-S|iron-sulfur|ferredoxin",
        ignore_case = TRUE)) ~"Resource competition",
      # 9. Antimicrobial resistance / defense
      str_detect(ko_name, regex(
        "beta-lactamase|bla|penicillin|D-Ala-D-Ala|vancomycin|imipenem|porin OprD|multidrug|drug resistance|efflux|MFS transporter.*resistance|RND|mdt|emr|bcr|marC|fosmidomycin|antibiotic|CAMP resistance|cationic antimicrobial",
        ignore_case = TRUE)) ~
        "Antimicrobial resistance / defense",
      # 10a. Motility
      str_detect(ko_name, regex(
        "flagellar|flagellin|\\bflg[A-Z]?\\b|\\bfli[A-Z]?\\b|\\bflh[A-Z]?\\b|chemotaxis|\\bche[A-Z]?\\b|\\bmotA\\b|\\bmotB\\b|twitching motility|\\bpilT\\b|\\bpilU\\b|\\biolB\\b",
        ignore_case = TRUE)) ~
        "Motility",
      # 10b. Biofilm / surface attachment
      str_detect(ko_name, regex(
        "biofilm|filamentous hemagglutinin|c-di-GMP|diguanylate|\\bdgc\\b|phosphodiesterase|exopolysaccharide|EPS|polysaccharide export|fimbr|adhesin|autotransporter|\\bpilA\\b|\\bpilB\\b|\\bpilC\\b|\\bpilQ\\b",
        ignore_case = TRUE)) ~
        "Biofilm / surface attachment",
      # 11a. Cell sensory/regulatory
      str_detect(ko_name, regex(
        "two-component system.*(OmpR|PhoB|PhoP|PhoQ|PhoR|ChvG|ChvI|RpfC|TctD|TctE)|sensor histidine kinase.*(PhoQ|PhoR|ChvG|TctE|RpfC)|response regulator.*(PhoB|PhoP|ChvI|TctD)",
        ignore_case = TRUE)) ~"Cell environmental sensing",
      # 11b. Cell envelope structure
      str_detect(ko_name, regex(
        "lipopolysaccharide|\\bLPS\\b|lipid A|peptidoglycan|murein|\\bmur[A-Z]?\\b|\\bmraY\\b|\\bmurJ\\b|\\blpt[A-Z]?\\b|\\blpx[A-Z]?\\b|outer membrane protein assembly|\\bbam[A-Z]?\\b|\\bSAM50\\b|\\bmla[A-Z]?\\b|cardiolipin|cyclopropane-fatty-acyl|phosphatidyl|lipoprotein|UDP-N-acetylglucosamine|N-acetylglucosamine|\\bglmS\\b|\\bglmM\\b|\\bglmU\\b",
        ignore_case = TRUE)) ~
        "Cell envelope structure",
      # 12. Microbe-microbe and plant-microbe interactions
      str_detect(ko_name, regex(
        "toxin|antitoxin|quorum|autoinducer|AHL|homoserine lactone|bacteriocin|colicin",
        ignore_case = TRUE)) ~"Microbe-microbe interactions",
      # 14. General stress / repair
      str_detect(ko_name, regex(
        "heat shock|cold shock|chaperone|Clp|Lon|Hsl|Hsp|DnaK|DnaJ|GroEL|GroES|starvation|carbon starvation|stringent|SpoT|LexA|SOS|DNA repair|mismatch repair|photolyase|restriction enzyme",
        ignore_case = TRUE)) ~"General stress / repair",
      TRUE ~ "Other / core metabolism"))

kos_collapsed_cat %>%count(ko_high_level_category, sort = TRUE)
# filter biologically irrelevant significant KOs, that ended up being significant because of normalization

kos_collapsed_cat.filt <- kos_collapsed_cat %>%filter(abs(mean_raw_diff_dry_minus_wet) > 0.5)
kos_collapsed_cat.filt.counts<- kos_collapsed_cat.filt %>%count(ko_high_level_category, sort = TRUE)

write_xlsx(kos_collapsed_cat, "results/sig.ko.sorted.in.cats.5_1_26.xlsx")
write_xlsx(kos_collapsed_cat.filt, "results/sig.ko.sorted.in.cats.filt.5_1_26.xlsx")

############## Figure 4b ###########################################
# Plot the counts of KO in each cats
# plot
ggplot(kos_collapsed_cat.filt.counts,aes(y=reorder(ko_high_level_category, n),x=n)) +
  geom_col()+
  theme_bw()+
  geom_text(aes(label = n), hjust = 1.5, size = 3, color="white") +
  scale_x_break(c(35, 215), scales = 0.3, expand =expansion(mult = c(0.02, 0.2))) +
  labs(x = "Number of significantly enriched KOs",y = "KO Groups")+
  theme(axis.text.y = element_text(size = 10),axis.text.x = element_text(size = 8, angle = 40, hjust = 0.8))
ggsave(filename= "bar.count.ko.in.groups.mod.con.filt.v2.pdf", device="pdf", units="mm", dpi=300, width=150, height=100, path="plots")

############################################
# Plot KOs grouped by catagories ###########
############################################
df<-read_excel("results/sig.ko.sorted.in.cats.filt.5_1_26.unit.names.xlsx")
unique(df$ko_high_level_category)
df.top <- df %>%filter(ko_high_level_category %in% c(
    "Secretion systems", #"Cell envelope structure",
    #"Resource competition",# "Antimicrobial resistance / defense",
    "Phage related",
    #"Metal / toxic ion stress",
    "Oxidative stress",# "Cell environmental sensing",
    "Osmotic homeostasis",
    #"Nitrogen metabolism",
    "Motility",
    "Microbe-microbe interactions",
    "General stress / repair",
    "Cytochrome / electron transport",
    "Biofilm / surface attachment"
   # "Antimicrobial resistance / defense"
   ))

# summary bubbleplot ####
############## Figure 4c ###########################################
#df.top <- df.top %>% arrange(mean_norm_diff_dry_minus_wet) %>% mutate(ko_name = factor(ko_name, levels = ko_name))
df.top.avg <- df.top %>%group_by(ko_collapsed_units) %>%summarise(
    mean_norm_diff_dry_minus_wet = mean(mean_norm_diff_dry_minus_wet, na.rm = TRUE),
    across(where(is.numeric) & !matches("mean_norm_diff_dry_minus_wet"),~ mean(.x, na.rm = TRUE)),
    across(where(is.character),~ paste(unique(.x), collapse = "; ")),.groups = "drop")
df.top.avg <- df.top.avg %>%arrange(mean_norm_diff_dry_minus_wet) %>%mutate(
    ko_collapsed_units = factor(ko_collapsed_units, levels = ko_collapsed_units))
df.top.avg <- df.top.avg %>% mutate(prevalence = pmax(percent_present_dry, percent_present_wet),
  direction = ifelse(mean_norm_diff_dry_minus_wet > 0, "Dry-enriched", "Wet-enriched"),
  neg_log10_p = -log10(p_adj_fdr), log2FC_dry_vs_wet = log2((mean_norm_dry) / (mean_norm_wet)))

ggplot(df.top.avg , aes(x = mean_norm_diff_dry_minus_wet, y = ko_collapsed_units, color=direction)) +
  geom_point(size=3) + geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(ko_high_level_category~., scales = "free_y",space = "free_y",
             labeller = as_labeller(c("Secretion systems"="Secretion systems","Phage related"="Phage",
                                      "Oxidative stress"="Oxidative stress","Osmotic homeostasis"="Osmotic homeostasis",
                                      "General stress / repair"="General stress","Motility"="Motility","Microbe-microbe interactions"="M-M",
                                      "Cytochrome / electron transport"="Cytochr.","Biofilm / surface attachment"="Biofilm/attachment")))+
  scale_color_manual(values=c("#aa4465", "#3d5a80"),limits=c("Dry-enriched", "Wet-enriched")) +
  labs(x = "Effect size: Normalized abundance difference (Dry-Wet)",y = "KEGG Ortholog",color = "Enriched") +
  theme_bw() +theme(axis.text.y = element_text(size = 8),axis.text.x = element_text(size = 8, angle = 20, hjust = 0.8),panel.grid.minor = element_blank(),
                    strip.background = element_rect(fill = "white"),strip.text = element_text(size = 8),panel.spacing = unit(0.2, "lines"))
ggsave(filename= "top.cat.sig.ko.in.groups.bubbleplot.mod.conf.filt.pdf", device="pdf", units="mm", dpi=300, width=200, height=250, path="plots/ko_groups")

####### extended figure ####
############## Supp Figure S5 ###########################################
df.ext <- df %>%filter(ko_high_level_category %in% c(#"Secretion systems", 
  "Cell envelope structure","Resource competition","Antimicrobial resistance / defense",
  #"Phage related", "Oxidative stress", "Osmotic homeostasis",
  "Metal / toxic ion stress","Cell environmental sensing","Nitrogen metabolism"
  #"Motility",#"Microbe-microbe interactions","General stress / repair","Cytochrome / electron transport","Biofilm / surface attachment"
))
df.ext.avg <- df.ext %>%group_by(ko_collapsed_units) %>%summarise(
  mean_norm_diff_dry_minus_wet = mean(mean_norm_diff_dry_minus_wet, na.rm = TRUE),
  across(where(is.numeric) & !matches("mean_norm_diff_dry_minus_wet"),~ mean(.x, na.rm = TRUE)),
  across(where(is.character),~ paste(unique(.x), collapse = "; ")),.groups = "drop")
df.ext.avg <- df.ext.avg %>%arrange(mean_norm_diff_dry_minus_wet) %>%mutate(
  ko_collapsed_units = factor(ko_collapsed_units, levels = ko_collapsed_units))
df.ext.avg <- df.ext.avg %>% mutate(prevalence = pmax(percent_present_dry, percent_present_wet),
                                    direction = ifelse(mean_norm_diff_dry_minus_wet > 0, "Dry-enriched", "Wet-enriched"),
                                    neg_log10_p = -log10(p_adj_fdr), log2FC_dry_vs_wet = log2((mean_norm_dry) / (mean_norm_wet)))
ggplot(df.ext.avg, aes(x = mean_norm_diff_dry_minus_wet, y = ko_collapsed_units, color=direction)) +
  geom_point(size=3) + geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(ko_high_level_category~., scales = "free_y",space = "free_y")+
  scale_color_manual(values=c("#aa4465", "#3d5a80"),limits=c("Dry-enriched", "Wet-enriched")) +
  labs(x = "Effect size: Normalized abundance difference (Dry-Wet)",y = "KEGG Ortholog",color = "Enriched") +
  theme_bw() +theme(axis.text.y = element_text(size = 8),axis.text.x = element_text(size = 8, angle = 20, hjust = 0.8),panel.grid.minor = element_blank(),
                    strip.background = element_rect(fill = "white"),strip.text = element_text(size = 8),panel.spacing = unit(0.2, "lines"))
ggsave(filename= "extended.cat.sig.ko.in.groups.bubbleplot.mod.conf.filt.pdf", device="pdf", units="mm", dpi=300, width=200, height=175, path="plots/ko_groups")
