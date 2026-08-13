ani_file <- "/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/ani/ani_within.tsv"
meta_file <- "/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/metadata/genome_metadata_with_files.tsv"
out_file <- "/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/results/ani_pair_summary.tsv"

ani <- read.delim(ani_file, header = FALSE,
                  col.names = c("query", "reference", "ANI", "fragments_mapped", "total_fragments"))
meta <- read.delim(meta_file)

clean_name <- function(x) {
  x <- basename(x)
  x <- sub("\\.fasta$", "", x)
  x <- sub("\\.fna$", "", x)
  x
}

ani$query <- clean_name(ani$query)
ani$reference <- clean_name(ani$reference)

meta_q <- meta[, c("strain", "site", "region", "precip_in")]
colnames(meta_q) <- c("query", "query_site", "query_region", "query_precip")

meta_r <- meta[, c("strain", "site", "region", "precip_in")]
colnames(meta_r) <- c("reference", "ref_site", "ref_region", "ref_precip")

ani <- merge(ani, meta_q, by = "query")
ani <- merge(ani, meta_r, by = "reference")

ani$pair_type <- ifelse(
  ani$query_site == ani$ref_site,
  paste0(ani$query_site, "_within"),
  apply(cbind(ani$query_site, ani$ref_site), 1, function(x) paste(sort(x), collapse = "_vs_"))
)

write.table(ani, out_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nANI by pair type:\n")
print(aggregate(ANI ~ pair_type, data = ani,
                FUN = function(x) c(n = length(x), min = min(x), median = median(x), max = max(x))))