#!/bin/bash
#SBATCH --job-name=iqtree_luteibacter
#SBATCH --output=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/logs/iqtree_%j.out
#SBATCH --error=/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/logs/iqtree_%j.err
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00

set -eo pipefail
set -x

module purge
module load IQ-TREE/2.2.2.6

cd /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/phylogeny/roary

which iqtree2 || which iqtree

if command -v iqtree2 >/dev/null 2>&1; then
    IQTREE_BIN="iqtree2"
elif command -v iqtree >/dev/null 2>&1; then
    IQTREE_BIN="iqtree"
else
    echo "ERROR: neither iqtree2 nor iqtree found after loading module"
    exit 1
fi

"$IQTREE_BIN" \
  -s core_gene_alignment.aln \
  -m MFP \
  -B 1000 \
  -T AUTO \
  --prefix luteibacter_core