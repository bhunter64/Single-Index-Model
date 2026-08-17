#!/bin/bash

set -euo pipefail

sample_sizes=(600 1000 1400)
beta_1_hrs=(1 2)
beta_2_hrs=(2.5 2.8)
beta_3_hrs=(2.5 2.8)
base_output_dir="${BASE_OUTPUT_DIR:-results/beta-grid-simulation}"
submission_log="${SUBMISSION_LOG:-${base_output_dir}/submitted-jobs.csv}"
combine_after="${COMBINE_AFTER:-1}"
job_ids=()

format_hr() {
  printf "%s" "$1" | tr "." "p"
}

mkdir -p "${base_output_dir}"
printf "job_id,data_n,beta_1_hr,beta_2_hr,beta_3_hr,output_dir\n"
printf "job_id,data_n,beta_1_hr,beta_2_hr,beta_3_hr,output_dir\n" > "${submission_log}"

for sample_size in "${sample_sizes[@]}"; do
  for beta_1_hr in "${beta_1_hrs[@]}"; do
    for beta_2_hr in "${beta_2_hrs[@]}"; do
      for beta_3_hr in "${beta_3_hrs[@]}"; do
        scenario_label="n${sample_size}-b1_$(format_hr "${beta_1_hr}")-b2_$(format_hr "${beta_2_hr}")-b3_$(format_hr "${beta_3_hr}")"
        output_dir="${base_output_dir}/${scenario_label}"

        job_id=$(SAMPLE_SIZE="${sample_size}" \
          BETA_1_HR="${beta_1_hr}" \
          BETA_2_HR="${beta_2_hr}" \
          BETA_3_HR="${beta_3_hr}" \
          SCENARIO_LABEL="${scenario_label}" \
          OUTPUT_DIR="${output_dir}" \
          sbatch --parsable scripts/submit-fixed-parameter-simulation-replicates.sbatch)
        job_ids+=("${job_id}")

        row=$(printf "%s,%s,%s,%s,%s,%s" \
          "${job_id}" \
          "${sample_size}" \
          "${beta_1_hr}" \
          "${beta_2_hr}" \
          "${beta_3_hr}" \
          "${output_dir}")
        printf "%s\n" "${row}"
        printf "%s\n" "${row}" >> "${submission_log}"
      done
    done
  done
done

printf "wrote %s\n" "${submission_log}" >&2

if [[ "${combine_after}" == "1" ]]; then
  dependency=$(IFS=:; printf "%s" "${job_ids[*]}")
  combine_job_id=$(BASE_OUTPUT_DIR="${base_output_dir}" \
    sbatch --parsable \
    --dependency="afterok:${dependency}" \
    scripts/combine-beta-grid-simulation-results.sbatch)
  printf "combine_job_id,%s\n" "${combine_job_id}" | tee -a "${submission_log}"
fi
