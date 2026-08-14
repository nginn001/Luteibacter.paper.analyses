# === CONFIG ===
DATA_DIR="/rhome/nginn001/bigdata/KU_Luteibacter/phage/all_phage_region_seq"   # where your *.fna live
OUT_DIR="${DATA_DIR}/../cluster_work"
ID=0.90        # identity threshold (try 0.95 too)
THREADS=2      # keep small on head node

mkdir -p "$OUT_DIR"
cd "$OUT_DIR" || exit 1

# 0) Load module(s)
module purge
module load cd-hit   # or cd-hit/4.8.1 (whatever your HPCC provides)

# 1) Combine your files and normalize headers to Strain|idx (unique)
cat /dev/null > all_regions.fna
idx=0
for f in "${DATA_DIR}"/*.fna; do
  strain=$(basename "$f"); strain=${strain%%.*}; strain=${strain%_phage_regions}
  awk -v S="$strain" '
    BEGIN{OFS=""}
    /^>/ {h++; print ">" S "|" h; next}
    {print}
  ' "$f" >> all_regions.fna
done

# 2) Cluster with cd-hit-est
#    For ID >=0.90, use -n 8; for 0.95 use -n 10.
NWORD=8; awk -v id="$ID" 'BEGIN{if (id>=0.95) print 10; else print 8}' >/dev/null && NWORD=$(awk -v id="$ID" 'BEGIN{print (id>=0.95?10:8)}')
cd-hit-est -i all_regions.fna -o clusters_${ID}.fna \
  -c "$ID" -n "$NWORD" -T "$THREADS" -M 0 -d 0

# Outputs:
#   clusters_${ID}.fna        (cluster centroids)
#   clusters_${ID}.fna.clstr  (cluster membership)

# 3) Split into per-cluster ID lists
awk '
  /^>Cluster/ {cl++}
  /, >/ {
    if (match($0,/>[^ ]+/)) {
      id=substr($0,RSTART+1,RLENGTH-1);
      printf "%s\n", id > sprintf("cluster_%05d.ids", cl)
    }
  }
' clusters_${ID}.fna.clstr

# 4) Make per-cluster FASTAs
# If you have seqkit, this is easiest:
if module avail 2>/dev/null | grep -qi seqkit; then
  module load seqkit
  for ids in cluster_*.ids; do
    c=${ids%.ids}
    seqkit grep -f "$ids" all_regions.fna > "${c}.fna"
  done
else
  # Pure-awk fallback to extract sequences by ID list
  # (reads all ids into a map, then streams the big FASTA once per cluster)
  for ids in cluster_*.ids; do
    c=${ids%.ids}
    awk -v IDS="$ids" '
      BEGIN{
        while ((getline < IDS)>0) keep[$0]=1;
        close(IDS)
      }
      /^>/{
        h=substr($0,2);
        printseq=(h in keep);
      }
      { if (printseq) print }
    ' all_regions.fna > "${c}.fna"
  done
fi

echo "Done. Per-cluster FASTAs: $(ls cluster_*.fna | wc -l)"

