#!/bin/bash
#SBATCH --job-name=annot_all_snps
#SBATCH --partition=epyc
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/logs/%x_%j.out
#SBATCH --error=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/logs/%x_%j.err

set -eo pipefail
set -x

CONDA_BASE="/opt/linux/rocky/8.x/x86_64/pkgs/miniconda3/py39_4.12.0"
ENV_NAME="/rhome/nginn001/.conda/envs/gwas_py"

unset -f conda || true
unset -f __conda_activate || true
unset -f __conda_exe || true
unset -f __conda_hashr || true
unset CONDA_EXE CONDA_PREFIX CONDA_PROMPT_MODIFIER CONDA_DEFAULT_ENV CONDA_SHLVL _CE_M _CE_CONDA || true

export PATH="${CONDA_BASE}/bin:${PATH}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

python /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/scripts/04_annotate_pyseer_snps.py