#!/bin/bash
#SBATCH --job-name=sim-fixed
#SBATCH --array=1-50
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=slurm/fixed_sim_%A_%a.out
#SBATCH --error=slurm/fixed_sim_%A_%a.err

set -euo pipefail

SAMPLE_SIZE="${SAMPLE_SIZE:-600}"
OUTPUT_DIR="${OUTPUT_DIR:-results/fixed-parameter-simulation-n${SAMPLE_SIZE}}"

mkdir -p slurm "${OUTPUT_DIR}"

module load r

Rscript scripts/run-fixed-parameter-simulation-replicate.R \
  --replicate-index "${SLURM_ARRAY_TASK_ID}" \
  --data-n "${SAMPLE_SIZE}" \
  --n-data-sets 10 \
  --output-dir "${OUTPUT_DIR}"
