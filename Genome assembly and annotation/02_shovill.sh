#!/bin/bash
#SBATCH --job-name=shovill
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 1
#SBATCH --time=2:00:00
#SBATCH --mem=60G
#SBATCH --output=%x_%A_%a.out
#SBATCH --array=1-96

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p sample_list.csv)
PREFIX=$(echo "$LINE" | cut -f 1 -d ',')
FASTQ1=$(echo "$LINE" | cut -f 2 -d ',')
FASTQ2=$(echo "$LINE" | cut -f 3 -d ',' | tr -d '\r\n')

CONTAINERDIR=/kuhpc/work/sjmac/observer/containers/
SCRATCHTMP=/kuhpc/scratch/sjmac/observer/tmp

apptainer exec ${CONTAINERDIR}/staphb-shovill-1.1.0-2022Dec.img shovill \
    --cpus 1 \
    --ram 60 \
    --opts "--cov-cutoff auto" \
    --outdir assembly/${PREFIX}/ \
    --tmpdir ${SCRATCHTMP} \
    --R1 reads/${PREFIX}.trimmed.R1.fastq.gz \
    --R2 reads/${PREFIX}.trimmed.R2.fastq.gz
