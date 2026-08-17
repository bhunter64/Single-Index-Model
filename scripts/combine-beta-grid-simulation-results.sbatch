#!/bin/bash
#SBATCH --job-name=sim-grid-combine
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=slurm/beta_grid_combine_%j.out
#SBATCH --error=slurm/beta_grid_combine_%j.err

set -euo pipefail

BASE_OUTPUT_DIR="${BASE_OUTPUT_DIR:-results/beta-grid-simulation}"

mkdir -p slurm "${BASE_OUTPUT_DIR}"

module load r

Rscript scripts/combine-beta-grid-simulation-results.R \
  "${BASE_OUTPUT_DIR}" \
  "${BASE_OUTPUT_DIR}/beta-grid-parameter-summary.csv" \
  "${BASE_OUTPUT_DIR}/beta-grid-convergence-summary.csv"
