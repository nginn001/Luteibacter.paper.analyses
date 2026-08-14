#!/bin/bash
#SBATCH --job-name=pyseer_subset
#SBATCH --partition=epyc
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/logs/%x_%j.out
#SBATCH --error=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/logs/%x_%j.err

set -eo pipefail
set -x

RUNDIR="$1"

CONDA_BASE="/opt/linux/rocky/8.x/x86_64/pkgs/miniconda3/py39_4.12.0"
PYSEER_ENV="/rhome/nginn001/.conda/envs/pyseer_env"

unset -f conda || true
unset -f __conda_activate || true
unset -f __conda_exe || true
unset -f __conda_hashr || true
unset CONDA_EXE CONDA_PREFIX CONDA_PROMPT_MODIFIER CONDA_DEFAULT_ENV CONDA_SHLVL _CE_M _CE_CONDA || true

export PATH="${CONDA_BASE}/bin:${PATH}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${PYSEER_ENV}"

which pyseer
pyseer --version

cd "${RUNDIR}"

# Gene presence/absence GWAS
pyseer \
  --phenotypes phenotypes.tsv \
  --pres genes.Rtab \
  --distances mash.tsv \
  --continuous \
  --min-af 0.01 \
  --max-af 0.99 \
  --max-missing 0.10 \
  --cpu "${SLURM_CPUS_PER_TASK}" \
  > genes_pyseer.tsv

# SNP GWAS
pyseer \
  --phenotypes phenotypes.tsv \
  --vcf variants.vcf.gz \
  --distances mash.tsv \
  --continuous \
  --min-af 0.01 \
  --max-af 0.99 \
  --max-missing 0.10 \
  --cpu "${SLURM_CPUS_PER_TASK}" \
  > snps_pyseer.tsv