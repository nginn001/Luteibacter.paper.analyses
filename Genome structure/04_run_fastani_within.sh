#!/bin/bash
#SBATCH --job-name=lute_fastani
#SBATCH --partition=epyc
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/logs/%x_%j.out
#SBATCH --error=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/logs/%x_%j.err

set -euo pipefail

module purge
module load fastani

PROJECT="/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS"
cd "$PROJECT"

mkdir -p structure/ani

fastANI \
  --ql structure/mash/assemblies_list.txt \
  --rl structure/mash/assemblies_list.txt \
  --matrix \
  -o structure/ani/ani_within.tsv \
  -t 16