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
BETA_1_HR="${BETA_1_HR:-1.5}"
BETA_2_HR="${BETA_2_HR:-1.2}"
BETA_3_HR="${BETA_3_HR:-2}"
SCENARIO_LABEL="${SCENARIO_LABEL:-n${SAMPLE_SIZE}}"
OUTPUT_DIR="${OUTPUT_DIR:-results/fixed-parameter-simulation-${SCENARIO_LABEL}}"

mkdir -p slurm "${OUTPUT_DIR}"

module load r

Rscript scripts/run-fixed-parameter-simulation-replicate.R \
  --replicate-index "${SLURM_ARRAY_TASK_ID}" \
  --data-n "${SAMPLE_SIZE}" \
  --n-data-sets 10 \
  --true-beta-1 "$(Rscript -e "cat(log(${BETA_1_HR}))")" \
  --true-beta-2 "$(Rscript -e "cat(log(${BETA_2_HR}))")" \
  --true-beta-3 "$(Rscript -e "cat(log(${BETA_3_HR}))")" \
  --scenario-label "${SCENARIO_LABEL}" \
  --output-dir "${OUTPUT_DIR}"
