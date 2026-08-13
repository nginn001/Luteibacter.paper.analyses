#!/bin/bash
#SBATCH --job-name=mk_roary_gffs
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/logs/%x_%j.out
#SBATCH --error=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/logs/%x_%j.err

set -eo pipefail
set -x

PROJECT="/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS"
cd "$PROJECT"

OUTDIR="structure/phylogeny/roary_input_gff"
mkdir -p "$OUTDIR"

count=0

for gff in annotations/*.gff; do
    [ -e "$gff" ] || continue

    strain=$(basename "$gff" .gff)
    fasta="assemblies/${strain}.fasta"
    out="${OUTDIR}/${strain}.gff"

    if [ ! -e "$fasta" ]; then
        echo "Skipping ${strain}: missing ${fasta}"
        continue
    fi

    # Remove any existing ##FASTA section from the annotation file if present
    awk '
        BEGIN{fasta=0}
        /^##FASTA/ {fasta=1}
        fasta==0 {print}
    ' "$gff" > "$out"

    echo "##FASTA" >> "$out"
    cat "$fasta" >> "$out"

    count=$((count + 1))
done

echo "Created ${count} Roary-compatible GFF files in ${OUTDIR}"