#!/usr/bin/env Rscript

parse_args <- function(args) {
  values <- list(
    replicate_index = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"),
    data_n = "600",
    n_data_sets = "10",
    output_dir = NA_character_,
    max_data_set_attempts = "100",
    samples = "20000",
    burn_in = "10000",
    thin = "10"
  )

  index <- 1L
  while (index <= length(args)) {
    arg <- args[[index]]
    if (arg %in% c(
      "--replicate-index",
      "--data-n",
      "--n-data-sets",
      "--output-dir",
      "--max-data-set-attempts",
      "--samples",
      "--burn-in",
      "--thin"
    )) {
      key <- gsub("-", "_", sub("^--", "", arg))
      values[[key]] <- args[[index + 1L]]
      index <- index + 2L
    } else {
      stop("unknown argument: ", arg, call. = FALSE)
    }
  }

  values$replicate_index <- as.integer(values$replicate_index)
  values$data_n <- as.integer(values$data_n)
  values$n_data_sets <- as.integer(values$n_data_sets)
  values$max_data_set_attempts <- as.integer(values$max_data_set_attempts)
  values$samples <- as.integer(values$samples)
  values$burn_in <- as.integer(values$burn_in)
  values$thin <- as.integer(values$thin)
  if (is.na(values$output_dir)) {
    values$output_dir <- sprintf("results/fixed-parameter-simulation-n%d", values$data_n)
  }
  values
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

set.seed(20260708 + args$replicate_index)

source_package_r_files <- function() {
  r_files <- list.files("R", pattern = "[.][rR]$", full.names = TRUE)
  for (r_file in r_files) {
    source(r_file)
  }
}

source_package_r_files()

min_group_size <- 60

true_beta <- c(log(1.5), log(1.2), log(2))
true_gamma <- c(1.5, 1.8)
gamma_start_values <- c(1, 1.3)

control <- default_mcmc_control(
  samples = args$samples,
  burn_in = args$burn_in,
  thin = args$thin,
  gamma_mean = c(0.5, 0.5),
  gamma_sd = 1,
  gamma_proposal_sd = 0.5
)

simulate_from_true_parameters <- function(true_beta, true_gamma) {
  simulate_threshold_data(
    n = args$data_n,
    beta = true_beta,
    gamma = true_gamma,
    baseline_hazard = 1,
    biomarker_correlation = 0,
    study_end_range = c(0.4, 3)
  )
}

bias_metric_row <- function(parameter_prefix, value, truth) {
  bias <- value - truth
  absolute_bias <- abs(bias)

  stats::setNames(
    c(
      value,
      bias,
      absolute_bias,
      sqrt(bias^2),
      bias^2
    ),
    c(
      paste0(parameter_prefix, "_value"),
      paste0(parameter_prefix, "_bias"),
      paste0(parameter_prefix, "_absolute_bias"),
      paste0(parameter_prefix, "_rmse"),
      paste0(parameter_prefix, "_squared_bias")
    )
  )
}

empty_bias_metric_row <- function(parameter_prefix) {
  stats::setNames(
    rep(NA_real_, 5),
    c(
      paste0(parameter_prefix, "_value"),
      paste0(parameter_prefix, "_bias"),
      paste0(parameter_prefix, "_absolute_bias"),
      paste0(parameter_prefix, "_rmse"),
      paste0(parameter_prefix, "_squared_bias")
    )
  )
}

coverage_metric_row <- function(parameter_prefix, interval, truth) {
  interval_lower <- interval[1]
  interval_upper <- interval[2]
  covered <- interval_lower <= truth && truth <= interval_upper

  stats::setNames(
    c(
      interval_lower,
      interval_upper,
      covered
    ),
    c(
      paste0(parameter_prefix, "_interval_lower"),
      paste0(parameter_prefix, "_interval_upper"),
      paste0(parameter_prefix, "_covered")
    )
  )
}

empty_coverage_metric_row <- function(parameter_prefix) {
  stats::setNames(
    rep(NA_real_, 3),
    c(
      paste0(parameter_prefix, "_interval_lower"),
      paste0(parameter_prefix, "_interval_upper"),
      paste0(parameter_prefix, "_covered")
    )
  )
}

parameter_metric_row <- function(parameter_prefix, value, interval, truth) {
  c(
    bias_metric_row(parameter_prefix, value, truth),
    coverage_metric_row(parameter_prefix, interval, truth)
  )
}

empty_parameter_metric_row <- function(parameter_prefix) {
  c(
    empty_bias_metric_row(parameter_prefix),
    empty_coverage_metric_row(parameter_prefix)
  )
}

fit_metric_columns <- function(summary, true_beta, true_gamma) {
  c(
    parameter_metric_row("beta_1", summary$beta[1], summary$beta_interval[, 1], true_beta[1]),
    parameter_metric_row("beta_2", summary$beta[2], summary$beta_interval[, 2], true_beta[2]),
    parameter_metric_row("beta_3", summary$beta[3], summary$beta_interval[, 3], true_beta[3]),
    parameter_metric_row("gamma_1", summary$gamma[1], summary$gamma_interval[, 1], true_gamma[1]),
    parameter_metric_row("gamma_2", summary$gamma[2], summary$gamma_interval[, 2], true_gamma[2])
  )
}

empty_fit_metric_columns <- function() {
  c(
    empty_parameter_metric_row("beta_1"),
    empty_parameter_metric_row("beta_2"),
    empty_parameter_metric_row("beta_3"),
    empty_parameter_metric_row("gamma_1"),
    empty_parameter_metric_row("gamma_2")
  )
}

base_result_columns <- function(replicate_index,
                                data_set_index,
                                data_set_attempt,
                                start_index,
                                true_beta,
                                true_gamma,
                                gamma_start,
                                data_summary,
                                skipped,
                                skip_reason,
                                converged,
                                error) {
  data.frame(
    replicate_index = replicate_index,
    data_n = args$data_n,
    data_set_index = data_set_index,
    global_data_set_index = (replicate_index - 1L) * args$n_data_sets + data_set_index,
    data_set_attempt = data_set_attempt,
    start_index = start_index,
    true_beta_1 = true_beta[1],
    true_beta_2 = true_beta[2],
    true_beta_3 = true_beta[3],
    true_gamma_1 = true_gamma[1],
    true_gamma_2 = true_gamma[2],
    gamma_start_1 = gamma_start[1],
    gamma_start_2 = gamma_start[2],
    treated_above = data_summary$group_sizes[["treated_above"]],
    treated_below = data_summary$group_sizes[["treated_below"]],
    control_above = data_summary$group_sizes[["control_above"]],
    control_below = data_summary$group_sizes[["control_below"]],
    skipped = skipped,
    skip_reason = skip_reason,
    converged = converged,
    error = error
  )
}

result_row <- function(base_columns, metric_columns) {
  data.frame(base_columns, as.list(metric_columns), check.names = FALSE)
}

skipped_data_set_row <- function(replicate_index,
                                 data_set_index,
                                 data_set_attempt,
                                 true_beta,
                                 true_gamma,
                                 data_summary,
                                 skip_reason) {
  result_row(
    base_result_columns(
      replicate_index = replicate_index,
      data_set_index = data_set_index,
      data_set_attempt = data_set_attempt,
      start_index = NA_integer_,
      true_beta = true_beta,
      true_gamma = true_gamma,
      gamma_start = c(NA_real_, NA_real_),
      data_summary = data_summary,
      skipped = TRUE,
      skip_reason = skip_reason,
      converged = FALSE,
      error = NA_character_
    ),
    empty_fit_metric_columns()
  )
}

fit_one_data_set <- function(replicate_index, data_set_index, data_set_attempt, data, true_beta, true_gamma, data_summary) {
  start_index <- 1L
  gamma_start <- gamma_start_values

  tryCatch(
    {
      posterior <- suppressWarnings(fit_threshold_model(
        data = data,
        control = control,
        lambda = 0,
        gamma_start = gamma_start
      ))
      summary <- summarize_mcmc(posterior)

      result_row(
        base_result_columns(
          replicate_index = replicate_index,
          data_set_index = data_set_index,
          data_set_attempt = data_set_attempt,
          start_index = start_index,
          true_beta = true_beta,
          true_gamma = true_gamma,
          gamma_start = gamma_start,
          data_summary = data_summary,
          skipped = FALSE,
          skip_reason = NA_character_,
          converged = TRUE,
          error = NA_character_
        ),
        fit_metric_columns(summary, true_beta, true_gamma)
      )
    },
    error = function(err) {
      result_row(
        base_result_columns(
          replicate_index = replicate_index,
          data_set_index = data_set_index,
          data_set_attempt = data_set_attempt,
          start_index = start_index,
          true_beta = true_beta,
          true_gamma = true_gamma,
          gamma_start = gamma_start,
          data_summary = data_summary,
          skipped = FALSE,
          skip_reason = NA_character_,
          converged = FALSE,
          error = conditionMessage(err)
        ),
        empty_fit_metric_columns()
      )
    }
  )
}

simulate_data_set <- function(replicate_index, data_set_index) {
  fit_rows <- list()

  for (data_set_attempt in seq_len(args$max_data_set_attempts)) {
    data <- simulate_from_true_parameters(true_beta, true_gamma)
    data_summary <- summarize_simulated_data(data, true_gamma)

    if (!all(data_summary$group_sizes >= min_group_size)) {
      fit_rows[[data_set_attempt]] <- skipped_data_set_row(
        replicate_index = replicate_index,
        data_set_index = data_set_index,
        data_set_attempt = data_set_attempt,
        true_beta = true_beta,
        true_gamma = true_gamma,
        data_summary = data_summary,
        skip_reason = "insufficient_group_size"
      )
    } else {
      fit_rows[[data_set_attempt]] <- fit_one_data_set(
        replicate_index = replicate_index,
        data_set_index = data_set_index,
        data_set_attempt = data_set_attempt,
        data = data,
        true_beta = true_beta,
        true_gamma = true_gamma,
        data_summary = data_summary
      )
    }

    if (isTRUE(fit_rows[[data_set_attempt]]$converged)) {
      break
    }
  }

  do.call(rbind, fit_rows)
}

study_rows <- vector("list", args$n_data_sets)
for (data_set_index in seq_len(args$n_data_sets)) {
  study_rows[[data_set_index]] <- simulate_data_set(args$replicate_index, data_set_index)
}

result <- do.call(rbind, study_rows)
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(
  args$output_dir,
  sprintf("fixed-parameter-simulation-replicate-%02d.csv", args$replicate_index)
)
utils::write.csv(result, output_file, row.names = FALSE)
print(result)
message("wrote ", output_file)
