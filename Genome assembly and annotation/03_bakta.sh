#!/bin/bash
#SBATCH --job-name=bakta
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 16
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --output=%x_%A_%a.out
#SBATCH --array=1-96

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p sample_list.csv)
PREFIX=$(echo "$LINE" | cut -f 1 -d ',')
FASTQ1=$(echo "$LINE" | cut -f 2 -d ',')
FASTQ2=$(echo "$LINE" | cut -f 3 -d ',' | tr -d '\r\n')

CONTAINERDIR=/kuhpc/work/sjmac/observer/containers/
SCRATCHTMP=/kuhpc/scratch/sjmac/observer/tmp

apptainer exec ${CONTAINERDIR}/quay.io-biocontainers-bakta-1.9.3--pyhdfd78af_0.img bakta \
    --db bakta_db/db \
    --prefix ${PREFIX} \
    --output annotation/${PREFIX} \
    --genus Luteibacter \
    --species spp \
    --strain ${PREFIX} \
    --locus-tag $(echo ${PREFIX} | sed 's/-//g') \
    --threads 16 \
    --tmp-dir ${SCRATCHTMP} \
    --compliant \
    assembly/${PREFIX}/contigs.fa
