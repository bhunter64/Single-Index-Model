#!/bin/bash
#SBATCH --job-name=sim-uniform
#SBATCH --array=1-150%50
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=slurm/uniform_sim_%A_%a.out
#SBATCH --error=slurm/uniform_sim_%A_%a.err

set -euo pipefail

# There are 50 array tasks per sample size and 10 successful datasets per
# task, giving 500 successful fitted datasets at each n when all targets finish.
SAMPLE_SIZES=(600 1000 1400)
REPLICATES_PER_SAMPLE_SIZE=50
DATA_SETS_PER_REPLICATE=10

TASK_OFFSET=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE_SIZE_INDEX=$((TASK_OFFSET / REPLICATES_PER_SAMPLE_SIZE))
REPLICATE_INDEX=$((TASK_OFFSET % REPLICATES_PER_SAMPLE_SIZE + 1))
SAMPLE_SIZE="${SAMPLE_SIZES[${SAMPLE_SIZE_INDEX}]}"

# Match the mean and variance of the normal N(0.5, 1) prior used in the
# original study. A uniform distribution with mean 0.5 and SD 1 has bounds
# 0.5 +/- sqrt(3).
GAMMA_MIN="${GAMMA_MIN:--1.2320508075688772}"
GAMMA_MAX="${GAMMA_MAX:-2.2320508075688772}"
SCENARIO_LABEL="uniform-prior-n${SAMPLE_SIZE}"
OUTPUT_DIR="results/uniform-prior-simulation-n${SAMPLE_SIZE}"

mkdir -p slurm "${OUTPUT_DIR}"

module load r

Rscript scripts/run-fixed-parameter-simulation-replicate.R \
  --replicate-index "${REPLICATE_INDEX}" \
  --data-n "${SAMPLE_SIZE}" \
  --n-data-sets "${DATA_SETS_PER_REPLICATE}" \
  --max-data-set-attempts 100 \
  --samples 20000 \
  --burn-in 10000 \
  --thin 10 \
  --true-beta-1 "$(Rscript -e 'cat(log(1.5))')" \
  --true-beta-2 "$(Rscript -e 'cat(log(1.2))')" \
  --true-beta-3 "$(Rscript -e 'cat(log(2))')" \
  --gamma-prior uniform \
  --gamma-min "${GAMMA_MIN}" \
  --gamma-max "${GAMMA_MAX}" \
  --scenario-label "${SCENARIO_LABEL}" \
  --output-dir "${OUTPUT_DIR}"
