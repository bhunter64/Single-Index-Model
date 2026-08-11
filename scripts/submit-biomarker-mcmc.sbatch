#!/bin/bash
#SBATCH --job-name=biomarker-mcmc
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=slurm/biomarker_mcmc_%j.out
#SBATCH --error=slurm/biomarker_mcmc_%j.err

set -euo pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
DEFAULT_DATA_FILE="${PROJECT_DIR}/data/cleaned_data.rda"
if [[ -f /data/cleaned_data.rda ]]; then
  DEFAULT_DATA_FILE="/data/cleaned_data.rda"
fi
DATA_FILE="${DATA_FILE:-${DEFAULT_DATA_FILE}}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DIR}/results/biomarker-mcmc-${SLURM_JOB_ID:-local}}"
CHAINS="${CHAINS:-4}"
SAMPLES="${SAMPLES:-5000}"
BURN_IN="${BURN_IN:-5000}"
THIN="${THIN:-5}"
SEED="${SEED:-20260811}"

REFERENCE_ARGS=()
if [[ -n "${REFERENCE_GAMMA:-}" ]]; then
  REFERENCE_ARGS+=(--reference-gamma "${REFERENCE_GAMMA}")
fi
if [[ -n "${REFERENCE_BETA:-}" ]]; then
  REFERENCE_ARGS+=(--reference-beta "${REFERENCE_BETA}")
fi

mkdir -p "${PROJECT_DIR}/slurm" "${OUTPUT_DIR}"
cd "${PROJECT_DIR}"

module load r

Rscript scripts/run-biomarker-mcmc.R \
  --data "${DATA_FILE}" \
  --output-dir "${OUTPUT_DIR}" \
  --biomarker-cols "bio30,bio3,bio37" \
  --biomarker-labels "HER-2,CA19-9,Axl" \
  --gamma-start "-0.6,0.5,0.5" \
  --beta-start "${BETA_START:-auto}" \
  --chains "${CHAINS}" \
  --workers "${SLURM_CPUS_PER_TASK:-1}" \
  --samples "${SAMPLES}" \
  --burn-in "${BURN_IN}" \
  --thin "${THIN}" \
  --seed "${SEED}" \
  "${REFERENCE_ARGS[@]}"
