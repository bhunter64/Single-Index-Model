#!/bin/bash
#SBATCH --job-name=sim-fixed-100
#SBATCH --array=1-10
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=slurm/fixed_sim_%A_%a.out
#SBATCH --error=slurm/fixed_sim_%A_%a.err

set -euo pipefail

mkdir -p slurm results/fixed-parameter-simulation

module load r

Rscript scripts/run-fixed-parameter-simulation-replicate.R \
  --replicate-index "${SLURM_ARRAY_TASK_ID}" \
  --n-data-sets 10 \
  --output-dir results/fixed-parameter-simulation
