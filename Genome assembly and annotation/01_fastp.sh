#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 6
#SBATCH --time=1:00:00
#SBATCH --mem=64G
#SBATCH --output=%x_%A_%a.out
#SBATCH --array=1-96

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p sample_list.csv)
PREFIX=$(echo "$LINE" | cut -f 1 -d ',')
FASTQ1=$(echo "$LINE" | cut -f 2 -d ',')
FASTQ2=$(echo "$LINE" | cut -f 3 -d ',' | tr -d '\r\n')

CONTAINERDIR=/kuhpc/work/sjmac/observer/containers/

apptainer exec ${CONTAINERDIR}/quay.io-biocontainers-fastp-0.23.4--hadf994f_2.img fastp \
    -i ${FASTQ1} \
    -I ${FASTQ2} \
    -o reads/${PREFIX}.trimmed.R1.fastq.gz \
    -O reads/${PREFIX}.trimmed.R2.fastq.gz \
    -q 15 \
    -u 40 \
    --trim_poly_g \
    -w 6 \
    --detect_adapter_for_pe \
    -j logs/fastp/${PREFIX}_fastp.json \
    -h logs/fastp/${PREFIX}_fastp.html \
    -R "${PREFIX} fastp report" \
    -V
