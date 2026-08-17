#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input_root <- if (length(args) >= 1L) args[[1L]] else "results/beta-grid-simulation"
parameter_summary_file <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  file.path(input_root, "beta-grid-parameter-summary.csv")
}
convergence_summary_file <- if (length(args) >= 3L) {
  args[[3L]]
} else {
  file.path(input_root, "beta-grid-convergence-summary.csv")
}

scenario_dirs <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
replicate_files_for <- function(input_dir) {
  list.files(
    input_dir,
    pattern = "^fixed-parameter-simulation-replicate-[0-9]+[.]csv$",
    full.names = TRUE
  )
}
scenario_dirs <- scenario_dirs[vapply(scenario_dirs, function(input_dir) {
  length(replicate_files_for(input_dir)) > 0L
}, logical(1))]

if (length(scenario_dirs) == 0L) {
  stop("no beta-grid simulation replicate CSV files found in ", input_root, call. = FALSE)
}

coverage_mean <- function(fits, coverage_column) {
  if (!coverage_column %in% names(fits)) {
    return(NA_real_)
  }

  mean(fits[[coverage_column]], na.rm = TRUE)
}

censoring_mean <- function(fits) {
  if (nrow(fits) == 0L) {
    return(NA_real_)
  }
  if (!"censoring_proportion" %in% names(fits)) {
    return(NA_real_)
  }

  mean(fits$censoring_proportion, na.rm = TRUE)
}

scenario_value <- function(fits, column) {
  if (nrow(fits) == 0L) {
    return(NA)
  }
  if (!column %in% names(fits)) {
    return(NA)
  }

  unique_values <- unique(fits[[column]])
  if (length(unique_values) == 0L) {
    return(NA)
  }
  unique_values[1]
}

summary_for <- function(fits, parameters, bias_columns, coverage_columns) {
  data.frame(
    parameter = parameters,
    bias = vapply(fits[, bias_columns, drop = FALSE], mean, numeric(1)),
    absolute_bias = vapply(fits[, bias_columns, drop = FALSE], function(bias) {
      mean(abs(bias))
    }, numeric(1)),
    rmse = vapply(fits[, bias_columns, drop = FALSE], function(bias) {
      sqrt(mean(bias^2))
    }, numeric(1)),
    coverage_probability = vapply(coverage_columns, function(coverage_column) {
      coverage_mean(fits, coverage_column)
    }, numeric(1)),
    bias_mc_se = vapply(fits[, bias_columns, drop = FALSE], function(bias) {
      if (length(bias) < 2L) {
        return(NA_real_)
      }
      stats::sd(bias) / sqrt(length(bias))
    }, numeric(1))
  )
}

summarize_scenario <- function(input_dir) {
  files <- replicate_files_for(input_dir)
  results <- do.call(rbind, lapply(files, utils::read.csv))

  combined_file <- file.path(input_dir, "fixed-parameter-simulation-combined.csv")
  utils::write.csv(results, combined_file, row.names = FALSE)

  successful_fits <- results[results$converged, , drop = FALSE]
  fitted_rows <- results[!results$skipped, , drop = FALSE]

  convergence_summary <- data.frame(
    scenario_label = scenario_value(results, "scenario_label"),
    data_n = scenario_value(results, "data_n"),
    true_beta_1 = scenario_value(results, "true_beta_1"),
    true_beta_2 = scenario_value(results, "true_beta_2"),
    true_beta_3 = scenario_value(results, "true_beta_3"),
    true_gamma_1 = scenario_value(results, "true_gamma_1"),
    true_gamma_2 = scenario_value(results, "true_gamma_2"),
    n_replicates = length(unique(results$replicate_index)),
    n_effective_data_sets = length(unique(results$global_data_set_index)),
    n_rows = nrow(results),
    n_skipped = sum(results$skipped),
    n_fit_attempted = nrow(fitted_rows),
    n_successful = nrow(successful_fits),
    convergence_rate = if (nrow(fitted_rows) == 0L) NA_real_ else nrow(successful_fits) / nrow(fitted_rows),
    censoring_proportion = censoring_mean(results),
    successful_fit_censoring_proportion = censoring_mean(successful_fits)
  )

  if (nrow(successful_fits) == 0L) {
    parameter_summary <- data.frame()
  } else {
    parameter_summary <- rbind(
      summary_for(
        successful_fits,
        c("beta_1", "beta_2", "beta_3"),
        c("beta_1_bias", "beta_2_bias", "beta_3_bias"),
        c("beta_1_covered", "beta_2_covered", "beta_3_covered")
      ),
      summary_for(
        successful_fits,
        c("gamma_1", "gamma_2"),
        c("gamma_1_bias", "gamma_2_bias"),
        c("gamma_1_covered", "gamma_2_covered")
      )
    )
    parameter_summary <- data.frame(
      scenario_label = scenario_value(results, "scenario_label"),
      data_n = scenario_value(results, "data_n"),
      true_beta_1 = scenario_value(results, "true_beta_1"),
      true_beta_2 = scenario_value(results, "true_beta_2"),
      true_beta_3 = scenario_value(results, "true_beta_3"),
      parameter_summary,
      check.names = FALSE
    )
  }

  list(
    convergence_summary = convergence_summary,
    parameter_summary = parameter_summary
  )
}

summaries <- lapply(scenario_dirs, summarize_scenario)
convergence_summary <- do.call(rbind, lapply(summaries, `[[`, "convergence_summary"))
parameter_summaries <- lapply(summaries, `[[`, "parameter_summary")
parameter_summaries <- parameter_summaries[vapply(parameter_summaries, nrow, integer(1)) > 0L]
parameter_summary <- if (length(parameter_summaries) == 0L) {
  data.frame()
} else {
  do.call(rbind, parameter_summaries)
}

utils::write.csv(convergence_summary, convergence_summary_file, row.names = FALSE)
utils::write.csv(parameter_summary, parameter_summary_file, row.names = FALSE)

message("wrote ", convergence_summary_file)
message("wrote ", parameter_summary_file)
print(convergence_summary)
print(parameter_summary)
