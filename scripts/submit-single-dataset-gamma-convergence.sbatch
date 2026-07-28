#!/bin/bash
#SBATCH --job-name=gamma-conv
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=slurm/gamma_convergence_%j.out
#SBATCH --error=slurm/gamma_convergence_%j.err

set -euo pipefail

SAMPLE_SIZE="${SAMPLE_SIZE:-600}"
SEED="${SEED:-20260728}"
SAMPLES="${SAMPLES:-20000}"
BURN_IN="${BURN_IN:-10000}"
THIN="${THIN:-10}"
OUTPUT_DIR="${OUTPUT_DIR:-results/single-dataset-gamma-convergence-n${SAMPLE_SIZE}}"

mkdir -p slurm "${OUTPUT_DIR}"

module load r

Rscript scripts/run-single-dataset-gamma-convergence.R \
  --data-n "${SAMPLE_SIZE}" \
  --seed "${SEED}" \
  --samples "${SAMPLES}" \
  --burn-in "${BURN_IN}" \
  --thin "${THIN}" \
  --output-dir "${OUTPUT_DIR}"
