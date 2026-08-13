#!/bin/bash
#SBATCH --job-name=lutei_roary
#SBATCH --partition=epyc
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/logs/%x_%j.out
#SBATCH --error=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/logs/%x_%j.err

set -eo pipefail
set -x

PROJECT="/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS"
CONDA_BASE="/opt/linux/rocky/8.x/x86_64/pkgs/miniconda3/py39_4.12.0"
ROARY_ENV="/opt/linux/rocky/8.x/x86_64/pkgs/roary/3.13.0"

cd "$PROJECT"
mkdir -p structure/phylogeny/roary

unset -f conda || true
unset -f __conda_activate || true
unset -f __conda_exe || true
unset -f __conda_hashr || true
unset CONDA_EXE CONDA_PREFIX CONDA_PROMPT_MODIFIER CONDA_DEFAULT_ENV CONDA_SHLVL _CE_M _CE_CONDA || true

export PATH="${CONDA_BASE}/bin:${PATH}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ROARY_ENV}"

which roary

mapfile -t gff_files < <(find structure/phylogeny/roary_input_gff -maxdepth 1 -name "*.gff" | sort)

echo "Found ${#gff_files[@]} Roary-input GFF files"

if [ "${#gff_files[@]}" -lt 2 ]; then
    echo "ERROR: Need at least 2 GFF files for Roary"
    exit 1
fi

roary \
  -e \
  -n \
  -z \
  -p "${SLURM_CPUS_PER_TASK}" \
  -f structure/phylogeny/roary \
  "${gff_files[@]}"