#!/bin/bash
#SBATCH --job-name=pfam_roary_rep
#SBATCH --output=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/pfam/pfam_roary_rep.%j.out
#SBATCH --error=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/pfam/pfam_roary_rep.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

set -euo pipefail

INPUT=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/kofamscan/roary_cluster_representatives.faa
OUTDIR=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/pfam
PFAM_DB=/bigdata/roperlab/nginn001/databases/pfam/Pfam-A.hmm

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "INPUT=${INPUT}"
echo "PFAM_DB=${PFAM_DB}"

test -f "${INPUT}"
test -f "${PFAM_DB}"

# If not already pressed:
# hmmpress "${PFAM_DB}"

hmmscan \
  --cpu "${SLURM_CPUS_PER_TASK}" \
  --domtblout "${OUTDIR}/pfam.domtblout" \
  "${PFAM_DB}" \
  "${INPUT}" \
  > "${OUTDIR}/pfam.hmmscan.stdout.txt"

echo "Pfam hmmscan finished."