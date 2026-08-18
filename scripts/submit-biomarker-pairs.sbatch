#!/bin/bash
#SBATCH --job-name=biomarker-pairs
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=slurm/biomarker_pairs_%j.out
#SBATCH --error=slurm/biomarker_pairs_%j.err

set -euo pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
DEFAULT_DATA_FILE="${PROJECT_DIR}/data/cleaned_data.rda"
if [[ -f /data/cleaned_data.rda ]]; then
  DEFAULT_DATA_FILE="/data/cleaned_data.rda"
fi

DATA_FILE="${DATA_FILE:-${DEFAULT_DATA_FILE}}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DIR}/results/biomarker-pairs-${SLURM_JOB_ID:-local}}"
CHAINS="${CHAINS:-4}"
SAMPLES="${SAMPLES:-5000}"
BURN_IN="${BURN_IN:-5000}"
THIN="${THIN:-5}"
SEED="${SEED:-20260818}"
WORKERS="${SLURM_CPUS_PER_TASK:-1}"

PAIR_NAMES=(
  "HER2_CA19-9"
  "HER2_Axl"
  "CA19-9_Axl"
)
PAIR_COLUMNS=(
  "bio30,bio3"
  "bio30,bio37"
  "bio3,bio37"
)
PAIR_LABELS=(
  "HER-2,CA19-9"
  "HER-2,Axl"
  "CA19-9,Axl"
)
PAIR_GAMMA_STARTS=(
  "-0.6,0.5"
  "-0.6,0.5"
  "0.5,0.5"
)

mkdir -p "${PROJECT_DIR}/slurm" "${OUTPUT_DIR}"
cd "${PROJECT_DIR}"

if command -v module >/dev/null 2>&1; then
  module load r
fi

for index in "${!PAIR_NAMES[@]}"; do
  pair_name="${PAIR_NAMES[$index]}"
  pair_output="${OUTPUT_DIR}/${pair_name}"
  pair_seed="$((SEED + index * 1000))"

  echo "Starting ${pair_name}"
  Rscript scripts/run-biomarker-mcmc.R \
    --data "${DATA_FILE}" \
    --output-dir "${pair_output}" \
    --biomarker-cols "${PAIR_COLUMNS[$index]}" \
    --biomarker-labels "${PAIR_LABELS[$index]}" \
    --gamma-start "${PAIR_GAMMA_STARTS[$index]}" \
    --gamma-prior-mean "0,0" \
    --gamma-prior-sd "1,1" \
    --gamma-proposal-sd "0.05,0.05" \
    --beta-start "auto" \
    --chains "${CHAINS}" \
    --workers "${WORKERS}" \
    --samples "${SAMPLES}" \
    --burn-in "${BURN_IN}" \
    --thin "${THIN}" \
    --seed "${pair_seed}"
  echo "Finished ${pair_name}"
done

combined_summary="${OUTPUT_DIR}/pairwise-parameter-summary.csv"
first_summary="${OUTPUT_DIR}/${PAIR_NAMES[0]}/parameter-summary.csv"
{
  IFS= read -r header
  printf '"biomarker_pair",%s\n' "${header}"
} < "${first_summary}" > "${combined_summary}"

for pair_name in "${PAIR_NAMES[@]}"; do
  awk -v pair="${pair_name}" 'NR > 1 { print "\"" pair "\"," $0 }' \
    "${OUTPUT_DIR}/${pair_name}/parameter-summary.csv" >> "${combined_summary}"
done

cat > "${OUTPUT_DIR}/README-pairwise-results.txt" <<EOF
This job fit all three two-biomarker combinations.

Pair directories:
  HER2_CA19-9: bio30 (HER-2), bio3 (CA19-9); gamma start -0.6, 0.5
  HER2_Axl: bio30 (HER-2), bio37 (Axl); gamma start -0.6, 0.5
  CA19-9_Axl: bio3 (CA19-9), bio37 (Axl); gamma start 0.5, 0.5

Each directory contains parameter-summary.csv, posterior samples, and PNG/PDF
trace and autocorrelation plots. pairwise-parameter-summary.csv combines all
three parameter summaries. Empirical bias and reference coverage remain blank
for observed data because true parameter values are unknown.
EOF

echo "All pairwise analyses completed"
echo "Results: ${OUTPUT_DIR}"
