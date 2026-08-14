#!/bin/bash
#SBATCH --job-name=snippy_all
#SBATCH --partition=epyc
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/logs/%x_%j.out
#SBATCH --error=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/logs/%x_%j.err

set -eo pipefail
set -x

ROOT="/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS"
REF="${ROOT}/assemblies/KNZ12-1B.fasta"
SAMPLES="${ROOT}/gwas/inputs/snippy_all_samples.tsv"
OUTBASE="${ROOT}/gwas/snippy/all"

CONDA_BASE="/opt/linux/rocky/8.x/x86_64/pkgs/miniconda3/py39_4.12.0"
SNIPPY_ENV="/opt/linux/rocky/8.x/x86_64/pkgs/snippy/4.6.0/env"

mkdir -p "${OUTBASE}"

# Clean inherited conda state
unset -f conda || true
unset -f __conda_activate || true
unset -f __conda_exe || true
unset -f __conda_hashr || true
unset CONDA_EXE CONDA_PREFIX CONDA_PROMPT_MODIFIER CONDA_DEFAULT_ENV CONDA_SHLVL _CE_M _CE_CONDA || true

export PATH="${CONDA_BASE}/bin:${PATH}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${SNIPPY_ENV}"

which snippy
snippy --version || true

[ -f "${REF}" ] || { echo "Missing reference: ${REF}"; exit 1; }
[ -f "${SAMPLES}" ] || { echo "Missing samples file: ${SAMPLES}"; exit 1; }

tail -n +2 "${SAMPLES}" | while IFS=$'\t' read -r STRAIN ASSEMBLY; do
    if [ "${STRAIN}" = "KNZ12-1B" ]; then
        echo "Skipping reference strain ${STRAIN}"
        continue
    fi

    [ -f "${ASSEMBLY}" ] || { echo "Missing assembly for ${STRAIN}: ${ASSEMBLY}"; exit 1; }

    OUTDIR="${OUTBASE}/${STRAIN}"
    rm -rf "${OUTDIR}"

    snippy \
      --outdir "${OUTDIR}" \
      --ref "${REF}" \
      --ctgs "${ASSEMBLY}" \
      --cpus "${SLURM_CPUS_PER_TASK}"
done