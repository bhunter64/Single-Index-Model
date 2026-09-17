#!/usr/bin/env Rscript

sample_sizes <- c(600L, 1000L, 1400L)
project_root <- normalizePath(".", mustWork = TRUE)
combine_script <- file.path(project_root, "scripts", "combine-fixed-parameter-simulation-results.R")
results_root <- file.path(project_root, "results")

parameter_summaries <- vector("list", length(sample_sizes))
convergence_summaries <- vector("list", length(sample_sizes))

for (index in seq_along(sample_sizes)) {
  data_n <- sample_sizes[index]
  input_dir <- file.path(results_root, sprintf("uniform-prior-simulation-n%d", data_n))
  combined_file <- file.path(input_dir, "uniform-prior-simulation-combined.csv")
  parameter_file <- file.path(input_dir, "uniform-prior-simulation-parameter-summary.csv")
  convergence_file <- file.path(input_dir, "uniform-prior-simulation-convergence-summary.csv")

  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(combine_script, input_dir, combined_file, parameter_file, convergence_file)
  )
  if (!identical(status, 0L)) {
    stop("failed to combine uniform-prior results for n = ", data_n, call. = FALSE)
  }

  parameter_summaries[[index]] <- utils::read.csv(parameter_file, check.names = FALSE)
  convergence_summaries[[index]] <- utils::read.csv(convergence_file, check.names = FALSE)
}

parameter_summary <- do.call(rbind, parameter_summaries)
convergence_summary <- do.call(rbind, convergence_summaries)

parameter_output <- file.path(results_root, "uniform-prior-simulation-parameter-summary.csv")
convergence_output <- file.path(results_root, "uniform-prior-simulation-convergence-summary.csv")
utils::write.csv(parameter_summary, parameter_output, row.names = FALSE)
utils::write.csv(convergence_summary, convergence_output, row.names = FALSE)

if (any(convergence_summary$n_successful != 500L)) {
  warning(
    "at least one sample size does not yet have exactly 500 successful fits; ",
    "inspect the convergence summary before reporting results",
    call. = FALSE
  )
}

print(parameter_summary)
print(convergence_summary)
message("wrote ", parameter_output)
message("wrote ", convergence_output)
