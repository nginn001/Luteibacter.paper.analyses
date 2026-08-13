#!/bin/bash
#SBATCH --job-name=lutei_mash
#SBATCH --partition=short
#SBATCH --time=01:58:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/logs/%x_%j.out
#SBATCH --error=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/logs/%x_%j.err

set -euo pipefail

module purge
module load mash/2.3

PROJECT="/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS"
cd "$PROJECT"

mkdir -p structure/mash

mash sketch \
  -l structure/mash/assemblies_list.txt \
  -o structure/mash/luteibacter_all \
  -p 8 \
  -s 10000

mash dist \
  structure/mash/luteibacter_all.msh \
  structure/mash/luteibacter_all.msh \
  > structure/mash/mash_dist.tsv