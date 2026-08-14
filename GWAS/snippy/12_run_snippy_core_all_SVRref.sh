#!/bin/bash
#SBATCH --job-name=SVR_snippycore_all
#SBATCH --partition=epyc
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/logs/%x_%j.out
#SBATCH --error=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/logs/%x_%j.err

set -eo pipefail
set -x

ROOT="/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS"
INDIR="${ROOT}/gwas/snippy/all_SVRref"
OUTDIR="${ROOT}/gwas/vcf/SVR_genomes_svr3_8d"

CONDA_BASE="/opt/linux/rocky/8.x/x86_64/pkgs/miniconda3/py39_4.12.0"
SNIPPY_ENV="/opt/linux/rocky/8.x/x86_64/pkgs/snippy/4.6.0/env"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

unset -f conda || true
unset -f __conda_activate || true
unset -f __conda_exe || true
unset -f __conda_hashr || true
unset CONDA_EXE CONDA_PREFIX CONDA_PROMPT_MODIFIER CONDA_DEFAULT_ENV CONDA_SHLVL _CE_M _CE_CONDA || true

export PATH="${CONDA_BASE}/bin:${PATH}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${SNIPPY_ENV}"

which snippy-core
which bgzip
which tabix

snippy-core --ref "${ROOT}/assemblies/SVR3-8D.fasta" "${INDIR}"/*

if [ -f core.vcf ]; then
    bgzip -f core.vcf
    tabix -f -p vcf core.vcf.gz
fi