#!/bin/bash
set -eo pipefail

ROOT="/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas"
SCRIPT="${ROOT}/scripts/02_run_pyseer_subset.sh"

for run in \
  01_all_nacl \
  01_all_peg \
  02_all_noclone_nacl \
  02_all_noclone_peg \
  03_sp1_nacl \
  03_sp1_peg \
  04_sp1_noclone_nacl \
  04_sp1_noclone_peg \
  05_SVR_nacl \
  05_SVR_peg \
  06_SVR_noclone_nacl \
  06_SVR_noclone_peg
do
  sbatch "${SCRIPT}" "${ROOT}/runs/${run}"
done