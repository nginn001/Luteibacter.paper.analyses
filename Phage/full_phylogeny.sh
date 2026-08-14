#!/usr/bin/env bash
#SBATCH --job-name=phage_full_phylo
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --mail-user=nginn001@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -p epyc
set -euo pipefail

########################
#        CONFIG        #
########################
# Folder containing ALL your per-strain FASTAs (e.g., TLI6-IE_phage_regions.fna, etc.)
DATA_DIR="/rhome/nginn001/bigdata/KU_Luteibacter/phage/all_phage_region_seq"   # <-- CHANGE ME
# Output folder
OUT_DIR="/rhome/nginn001/bigdata/KU_Luteibacter/phage/phylo_full"                          # <-- CHANGE ME (or leave; it will be created)

# Modules (adjust names to your HPCC)
MAFFT_MOD="mafft/7.505"
TRIMAL_MOD="trimal"       # e.g. trimal/1.4 on some clusters
IQTREE_MOD="iqtree/2"     # e.g. iqtree/2.2.x

# Alignment preset: SAFE | BALANCED | THOROUGH
PRESET="${PRESET:-SAFE}"
# Strand flip: on|off  (turn off if memory is still tight)
ADJ="${ADJ:-on}"
################################

THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_NUM_THREADS="${THREADS}"

echo "[INFO] Job $SLURM_JOB_ID  threads=$THREADS  mem=${SLURM_MEM_PER_NODE:-$SLURM_MEM_PER_CPU}  preset=$PRESET  adj=$ADJ"
echo "[INFO] DATA_DIR=$DATA_DIR"
echo "[INFO] OUT_DIR=$OUT_DIR"
mkdir -p "$OUT_DIR"

module purge || true

have_module(){ module avail 2>&1 | grep -qiE "$1"; }
load_or_die(){ have_module "$1" && module load "$1" || { echo "[ERR] missing module: $1"; exit 1; }; }

load_or_die "$MAFFT_MOD"
have_module "$TRIMAL_MOD" && module load "$TRIMAL_MOD" && trimal_ok=1 || trimal_ok=0
load_or_die "$IQTREE_MOD"

# 0) Combine sequences with unique headers
shopt -s nullglob
mapfile -t FASTAS < <(ls -1 "$DATA_DIR"/*.fna 2>/dev/null || true)
(( ${#FASTAS[@]} > 0 )) || { echo "[ERR] no .fna in $DATA_DIR"; exit 1; }

COMBINED="$OUT_DIR/all_regions.fna"
awk '
  BEGIN{OFS=""}
  FILENAME!=prev {
    prev=FILENAME
    split(FILENAME,a,"/"); fname=a[length(a)]
    strain=fname; sub(/\..*$/,"",strain); sub(/_phage_regions$/,"",strain)
    idx=0
  }
  /^>/ { idx++; print ">" strain "|" idx; next }
  { print }
' "${FASTAS[@]}" > "$COMBINED"
echo "[INFO] sequences: $(grep -c '^>' "$COMBINED")"

# 1) Choose MAFFT flags by preset
# SAFE:   very low memory (FFT-NS-2 + parttree, no iter. refinement)
# BALANCED: a bit slower, still reasonable RAM (some refinement)
# THOROUGH: DO NOT use unless you have >64–128G; closer to your previous G-INS-i style
case "$PRESET" in
  SAFE)
    MAFFT_FLAGS=(--parttree --retree 2 --maxiterate 0) ;;    # FFT-NS-2
  BALANCED)
    MAFFT_FLAGS=(--parttree --retree 4 --maxiterate 2) ;;    # a touch more refinement
  THOROUGH)
    MAFFT_FLAGS=(--retree 4 --maxiterate 100) ;;             # can be expensive; avoid on 16G
  *)
    echo "[WARN] unknown PRESET=$PRESET; defaulting SAFE"
    MAFFT_FLAGS=(--parttree --retree 2 --maxiterate 0) ;;
esac

# Strand direction check
if [[ "$ADJ" == "on" ]]; then
  MAFFT_FLAGS=(--adjustdirectionaccurately "${MAFFT_FLAGS[@]}")
fi

ALN="$OUT_DIR/all.aln.fasta"
echo "[INFO] MAFFT flags: ${MAFFT_FLAGS[*]}"
echo "[INFO] Aligning -> $ALN"
mafft --thread "$THREADS" "${MAFFT_FLAGS[@]}" "$COMBINED" > "$ALN"

# 2) Trim
ALN_FOR_TREE="$ALN"
if (( trimal_ok == 1 )); then
  TRIM="$OUT_DIR/all.trim.fasta"
  echo "[INFO] trimAl -> $TRIM"
  trimal -in "$ALN" -out "$TRIM" -automated1
  ALN_FOR_TREE="$TRIM"
else
  echo "[INFO] trimAl not loaded; using untrimmed alignment."
fi

# 3) IQ-TREE
echo "[INFO] IQ-TREE2 on $ALN_FOR_TREE"
(
  cd "$OUT_DIR"
  iqtree2 -s "$(basename "$ALN_FOR_TREE")" \
          -m MFP \
          -B 1000 -alrt 1000 \
          -nt "$THREADS" \
          --prefix full
)

echo "[INFO] Done. Tree: $OUT_DIR/full.treefile"
