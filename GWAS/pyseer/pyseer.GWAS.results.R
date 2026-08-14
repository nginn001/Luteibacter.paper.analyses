# Luteibacter pyseer GWAS results sorting and visualization
# Nichole Ginnan (nginn001@ucr.edu)
# April 20, 2026
##### Load libs ####
library(readr)
library(tidyverse)
library(RColorBrewer)
library(writexl)

#### Load Data ####
setwd("/Users/nicholeginnan/Documents/KU/Luteibacter/april.26.genome.analyses/gwas_pyseer/")
# SVR no clones
snp.peg.svr.nc<- read_delim("output_files/GWAS_results/06_SVR_noclone_peg.snps_pyseer_annotated.tsv")
snp.nacl.svr.nc<- read_delim("output_files/GWAS_results/06_SVR_noclone_nacl.snps_pyseer_annotated.tsv")
gene.peg.svr.nc<- read_delim("output_files/GWAS_results/06_SVR_noclone_peg.genes_pyseer.tsv")
gene.nacl.svr.nc<- read_delim("output_files/GWAS_results/06_SVR_noclone_nacl.genes_pyseer.tsv")
# SVR
snp.peg.svr<- read_delim("output_files/GWAS_results/05_SVR_peg.snps_pyseer_annotated.tsv")
gene.peg.svr<- read_delim("output_files/GWAS_results/05_SVR_peg.genes_pyseer.tsv")
gene.nacl.svr<- read_delim("output_files/GWAS_results/05_SVR_nacl.genes_pyseer.tsv")
snp.nacl.svr<- read_delim("output_files/GWAS_results/05_SVR_nacl.snps_pyseer_annotated.tsv")
# Sp1 no clones
gene.peg.sp1.nc<- read_delim("output_files/GWAS_results/04_sp1_noclone_peg.genes_pyseer.tsv")
snp.peg.sp1.nc<- read_delim("output_files/GWAS_results/04_sp1_noclone_peg.snps_pyseer_annotated.tsv")
gene.nacl.sp1.nc<- read_delim("output_files/GWAS_results/04_sp1_noclone_nacl.genes_pyseer.tsv")
snp.nacl.sp1.nc<- read_delim("output_files/GWAS_results/04_sp1_noclone_nacl.snps_pyseer_annotated.tsv")
# Sp1
snp.peg.sp1<- read_delim("output_files/GWAS_results/03_sp1_peg.snps_pyseer_annotated.tsv")
gene.peg.sp1<- read_delim("output_files/GWAS_results/03_sp1_peg.genes_pyseer.tsv")
snp.nacl.sp1<- read_delim("output_files/GWAS_results/03_sp1_nacl.snps_pyseer_annotated.tsv")
gene.nacl.sp1<- read_delim("output_files/GWAS_results/03_sp1_nacl.genes_pyseer.tsv")
# All no clones
snp.peg.all.nc<- read_delim("output_files/GWAS_results/02_all_noclone_peg.snps_pyseer_annotated.tsv")
gene.peg.all.nc<- read_delim("output_files/GWAS_results/02_all_noclone_peg.genes_pyseer.tsv")
snp.nacl.all.nc<- read_delim("output_files/GWAS_results/02_all_noclone_nacl.snps_pyseer_annotated.tsv")
gene.nacl.all.nc<- read_delim("output_files/GWAS_results/02_all_noclone_nacl.genes_pyseer.tsv")
# All
snp.peg.all<- read_delim("output_files/GWAS_results/01_all_peg.snps_pyseer_annotated.tsv")
gene.peg.all<- read_delim("output_files/GWAS_results/01_all_peg.genes_pyseer.tsv")
snp.nacl.all<- read_delim("output_files/GWAS_results/01_all_nacl.snps_pyseer_annotated.tsv")
gene.nacl.all<- read_delim("output_files/GWAS_results/01_all_nacl.genes_pyseer.tsv")
##################################################################################
##### FDR adjustment of the `lrt-pvalue` for each results file ###################
##################################################################################
# create function #
add_gwas_stats <- function(df) {
  df %>%
    mutate(p = as.numeric(`lrt-pvalue`),
      fdr = p.adjust(p, method = "BH"),
      neglog10p = -log10(p),
      neglog10fdr = -log10(fdr),
      beta = as.numeric(beta),
      abs_beta = abs(beta),
      maf = pmin(af, 1 - af),
      is_rare = maf < 0.05,
      unstable_beta = abs(beta) > 100 | is.na(beta),
      high_se = `beta-std-err` > 10) %>% arrange(fdr, p)}
# apply the function to the list of files #
gwas_list <- list(snp.peg.svr.nc, snp.nacl.svr.nc, gene.peg.svr.nc, gene.nacl.svr.nc,
  snp.peg.svr, gene.peg.svr, gene.nacl.svr, snp.nacl.svr,gene.peg.sp1.nc,
  snp.peg.sp1.nc, gene.nacl.sp1.nc, snp.nacl.sp1.nc, snp.peg.sp1, gene.peg.sp1,
  snp.nacl.sp1, gene.nacl.sp1, snp.peg.all.nc, gene.peg.all.nc, snp.nacl.all.nc,
  gene.nacl.all.nc, snp.peg.all, gene.peg.all, snp.nacl.all, gene.nacl.all)

gwas_list <- lapply(gwas_list, add_gwas_stats)

# Reassign names #
list2env(setNames(gwas_list, c("snp.peg.svr.nc","snp.nacl.svr.nc","gene.peg.svr.nc","gene.nacl.svr.nc",
  "snp.peg.svr","gene.peg.svr","gene.nacl.svr","snp.nacl.svr","gene.peg.sp1.nc","snp.peg.sp1.nc",
  "gene.nacl.sp1.nc","snp.nacl.sp1.nc","snp.peg.sp1","gene.peg.sp1","snp.nacl.sp1","gene.nacl.sp1",
  "snp.peg.all.nc","gene.peg.all.nc","snp.nacl.all.nc","gene.nacl.all.nc",
  "snp.peg.all","gene.peg.all","snp.nacl.all","gene.nacl.all")), envir = .GlobalEnv)
##################################################################################
##### Add metadata, and combine results #########################################
#################################################################################
# create function
add_metadata <- function(df, run, phenotype, variant_type) {
  df %>% mutate(run = run, phenotype = phenotype, variant_type = variant_type)}

snp.peg.svr.nc <- add_metadata(snp.peg.svr.nc, "svr.nc", "peg", "snp")
snp.nacl.svr.nc <- add_metadata(snp.nacl.svr.nc, "svr.nc", "nacl", "snp")
gene.peg.svr.nc <- add_metadata(gene.peg.svr.nc, "svr.nc", "peg", "gene")
gene.nacl.svr.nc <- add_metadata(gene.nacl.svr.nc, "svr.nc", "nacl", "gene")
write_xlsx(snp.peg.svr.nc,"output_files/GWAS_results/FDR.snp.peg.svr.nc.xlsx")
write_xlsx(snp.nacl.svr.nc,"output_files/GWAS_results/FDR.snp.nacl.svr.nc.xlsx")
write_xlsx(gene.peg.svr.nc,"output_files/GWAS_results/FDR.gene.peg.svr.nc.xlsx")
write_xlsx(gene.nacl.svr.nc,"output_files/GWAS_results/FDR.gene.nacl.svr.nc.xlsx")

snp.peg.svr <- add_metadata(snp.peg.svr, "svr", "peg", "snp")
gene.peg.svr <- add_metadata(gene.peg.svr, "svr", "peg", "gene")
gene.nacl.svr <- add_metadata(gene.nacl.svr, "svr", "nacl", "gene")
snp.nacl.svr <- add_metadata(snp.nacl.svr, "svr", "nacl", "snp")
write_xlsx(snp.peg.svr,"output_files/GWAS_results/FDR.snp.peg.svr.xlsx")
write_xlsx(snp.nacl.svr,"output_files/GWAS_results/FDR.snp.nacl.svr.xlsx")
write_xlsx(gene.peg.svr,"output_files/GWAS_results/FDR.gene.peg.svr.xlsx")
write_xlsx(gene.nacl.svr,"output_files/GWAS_results/FDR.gene.nacl.svr.xlsx")

gene.peg.sp1.nc <- add_metadata(gene.peg.sp1.nc, "sp1.nc", "peg", "gene")
snp.peg.sp1.nc <- add_metadata(snp.peg.sp1.nc, "sp1.nc", "peg", "snp")
gene.nacl.sp1.nc <- add_metadata(gene.nacl.sp1.nc, "sp1.nc", "nacl", "gene")
snp.nacl.sp1.nc <- add_metadata(snp.nacl.sp1.nc, "sp1.nc", "nacl", "snp")
write_xlsx(snp.peg.sp1.nc,"output_files/GWAS_results/FDR.snp.peg.sp1.nc.xlsx")
write_xlsx(snp.nacl.sp1.nc,"output_files/GWAS_results/FDR.snp.nacl.sp1.nc.xlsx")
write_xlsx(gene.peg.sp1.nc,"output_files/GWAS_results/FDR.gene.peg.sp1.nc.xlsx")
write_xlsx(gene.nacl.sp1.nc,"output_files/GWAS_results/FDR.gene.nacl.sp1.nc.xlsx")

snp.peg.sp1 <- add_metadata(snp.peg.sp1, "sp1", "peg", "snp")
gene.peg.sp1 <- add_metadata(gene.peg.sp1, "sp1", "peg", "gene")
snp.nacl.sp1 <- add_metadata(snp.nacl.sp1, "sp1", "nacl", "snp")
gene.nacl.sp1 <- add_metadata(gene.nacl.sp1, "sp1", "nacl", "gene")
write_xlsx(snp.peg.sp1,"output_files/GWAS_results/FDR.snp.peg.sp1.xlsx")
write_xlsx(snp.nacl.sp1,"output_files/GWAS_results/FDR.snp.nacl.sp1.xlsx")
write_xlsx(gene.peg.sp1,"output_files/GWAS_results/FDR.gene.peg.sp1.xlsx")
write_xlsx(gene.nacl.sp1,"output_files/GWAS_results/FDR.gene.nacl.sp1.xlsx")

snp.peg.all.nc <- add_metadata(snp.peg.all.nc, "all.nc", "peg", "snp")
gene.peg.all.nc <- add_metadata(gene.peg.all.nc, "all.nc", "peg", "gene")
snp.nacl.all.nc <- add_metadata(snp.nacl.all.nc, "all.nc", "nacl", "snp")
gene.nacl.all.nc <- add_metadata(gene.nacl.all.nc, "all.nc", "nacl", "gene")
write_xlsx(snp.peg.all.nc,"output_files/GWAS_results/FDR.snp.peg.all.nc.xlsx")
write_xlsx(snp.nacl.all.nc,"output_files/GWAS_results/FDR.snp.nacl.all.nc.xlsx")
write_xlsx(gene.peg.all.nc,"output_files/GWAS_results/FDR.gene.peg.all.nc.xlsx")
write_xlsx(gene.nacl.all.nc,"output_files/GWAS_results/FDR.gene.nacl.all.nc.xlsx")

snp.peg.all <- add_metadata(snp.peg.all, "all", "peg", "snp")
gene.peg.all <- add_metadata(gene.peg.all, "all", "peg", "gene")
snp.nacl.all <- add_metadata(snp.nacl.all, "all", "nacl", "snp")
gene.nacl.all <- add_metadata(gene.nacl.all, "all", "nacl", "gene")
write_xlsx(snp.peg.all,"output_files/GWAS_results/FDR.snp.peg.all.xlsx")
write_xlsx(snp.nacl.all,"output_files/GWAS_results/FDR.snp.nacl.all.xlsx")
write_xlsx(gene.peg.all,"output_files/GWAS_results/FDR.gene.peg.all.xlsx")
write_xlsx(gene.nacl.all,"output_files/GWAS_results/FDR.gene.nacl.all.xlsx")

##### look at top-hits and sort in excel
# plot candidates ####
cand<-read_excel("GWAS.results.table.xlsx", sheet = "candidate_metadata")
gene_order <- c(
  "epsF","fadN","uvrC","ybbN",
  "atpI","DLAT, aceF, pdhC","xylB",
  "espI","pilQ","pilR, pehR",
  "gp04_gc_4853","gp05_gc_vcpB","nicK_gc_pri",
  "gp04_gc_4854","gp05_gc_4852","nicK_gc_4859",
  "group_12331","group_12335","group_247","group_81","group_9106",
  "group_1374","group_7399","group_9821","group_7655")

cand <- cand %>% mutate(Kofam_gene_name = factor(Kofam_gene_name, levels = gene_order))

ggplot(cand, aes(x = Kofam_gene_name, y = beta, color = GWAS_strains, shape = stressor)) +
  geom_pointrange(aes(ymin = beta - `beta-std-err`, ymax = beta + `beta-std-err`),
    position = position_dodge(width = 0.5), fatten = 8, linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Humid" = "#3d5a80","Semi-arid" = "#aa4465","All" = "#f7b557")) +
  facet_grid(. ~ variant_type, scales = "free_x", space = "free_x",
             labeller = labeller(variant_type = c("snp" = "SNPs","gene" = "Gene Presence/Absence"))) +
  ylab("Effect size: Beta coefficient") + xlab("Gene name") +
  theme_bw() + theme(axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 40, hjust = 0.9),
    panel.grid.minor = element_blank(),strip.background = element_rect(fill = "white"))
ggsave(filename= "GWAS.candidates.point.range.beta.pdf", device="pdf", units="mm", dpi=300, width=275, height=100, path="plots/")



