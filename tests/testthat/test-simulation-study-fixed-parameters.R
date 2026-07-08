test_that("threshold model records fixed-parameter simulation study results", {
  set.seed(20260708)

  n_data_sets <- 10 # Number of data sets to simulate for the simulation
  min_group_size <- 60 # Appropriate portions for convergence
  max_start_attempts <- 100 # How many times to try fitting each set
  n_biomarkers <- 2 # Number of biomarkers in each set

  true_beta <- c(log(1.5), log(1.2), log(2))
  true_gamma <- c(1.5, 1.8)
  gamma_start_values <- c(1, 1.3)

  # MCMC control data
  control <- default_mcmc_control(
    samples = 20000,
    burn_in = 10000,
    thin = 10,
    gamma_mean = c(0.5, 0.5),
    gamma_sd = 1,
    gamma_proposal_sd = 0.5
  )

  simulate_from_true_parameters <- function(true_beta, true_gamma) {
    simulate_threshold_data(
      n = 600,
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

  fit_metric_columns <- function(beta_values, gamma_values, true_beta, true_gamma) {
    c(
      bias_metric_row("beta_1", beta_values[1], true_beta[1]),
      bias_metric_row("beta_2", beta_values[2], true_beta[2]),
      bias_metric_row("beta_3", beta_values[3], true_beta[3]),
      bias_metric_row("gamma_1", gamma_values[1], true_gamma[1]),
      bias_metric_row("gamma_2", gamma_values[2], true_gamma[2])
    )
  }

  empty_fit_metric_columns <- function() {
    c(
      empty_bias_metric_row("beta_1"),
      empty_bias_metric_row("beta_2"),
      empty_bias_metric_row("beta_3"),
      empty_bias_metric_row("gamma_1"),
      empty_bias_metric_row("gamma_2")
    )
  }

  base_result_columns <- function(data_set_index,
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
      data_set_index = data_set_index,
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

  # If a data set can't get a start after `max_start_attempts` of gamma, it gets a skipped-data-set
  skipped_data_set_row <- function(data_set_index, true_beta, true_gamma, data_summary, skip_reason) {
    result_row(
      base_result_columns(
        data_set_index = data_set_index,
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

  # Helper (wrapper) to make each data set
  simulate_data_set <- function(data_set_index) {
    data <- simulate_from_true_parameters(true_beta, true_gamma)
    data_summary <- summarize_simulated_data(data, true_gamma)

    if (!all(data_summary$group_sizes >= min_group_size)) {
      return(skipped_data_set_row(
        data_set_index = data_set_index,
        true_beta = true_beta,
        true_gamma = true_gamma,
        data_summary = data_summary,
        skip_reason = "insufficient_group_size"
      ))
    }

    fit_one_data_set(
      data_set_index = data_set_index,
      data = data,
      true_beta = true_beta,
      true_gamma = true_gamma,
      data_summary = data_summary
    )
  }

  fit_one_data_set <- function(data_set_index, data, true_beta, true_gamma, data_summary) {
    fit_rows <- list()

    for (start_index in seq_len(max_start_attempts)) {
      gamma_start <- gamma_start_values

      fit_rows[[start_index]] <- tryCatch(
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
              data_set_index = data_set_index,
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
            fit_metric_columns(summary$beta, summary$gamma, true_beta, true_gamma)
          )
        },
        error = function(err) {
          result_row(
            base_result_columns(
              data_set_index = data_set_index,
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

      if (isTRUE(fit_rows[[start_index]]$converged)) {
        break
      }
    }

    do.call(rbind, fit_rows)
  }

  study_rows <- vector("list", n_data_sets)
  for (data_set_index in seq_len(n_data_sets)) {
    study_rows[[data_set_index]] <- simulate_data_set(data_set_index)
  }

  # Combine the sim results, extract the successful sets only
  study_summary <- do.call(rbind, study_rows)
  successful_fits <- study_summary[study_summary$converged, , drop = FALSE]
  fitted_rows <- study_summary[!study_summary$skipped, , drop = FALSE]

  # Extract the data needed to review per-trial values and biases
  trial_results <- successful_fits[, c(
    "data_set_index",
    "true_beta_1",
    "true_beta_2",
    "true_beta_3",
    "true_gamma_1",
    "true_gamma_2",
    "beta_1_value",
    "beta_2_value",
    "beta_3_value",
    "gamma_1_value",
    "gamma_2_value",
    "beta_1_bias",
    "beta_2_bias",
    "beta_3_bias",
    "gamma_1_bias",
    "gamma_2_bias",
    "beta_1_absolute_bias",
    "beta_2_absolute_bias",
    "beta_3_absolute_bias",
    "gamma_1_absolute_bias",
    "gamma_2_absolute_bias",
    "beta_1_rmse",
    "beta_2_rmse",
    "beta_3_rmse",
    "gamma_1_rmse",
    "gamma_2_rmse",
    "beta_1_squared_bias",
    "beta_2_squared_bias",
    "beta_3_squared_bias",
    "gamma_1_squared_bias",
    "gamma_2_squared_bias"
  ), drop = FALSE]

  # Helpful info for the following biases calculations
  n_successful_fits <- nrow(successful_fits)
  n_fitted_rows <- nrow(fitted_rows)
  beta_bias_columns <- c("beta_1_bias", "beta_2_bias", "beta_3_bias")
  gamma_bias_columns <- c("gamma_1_bias", "gamma_2_bias")

  summarize_empirical_bias <- function(fits, parameters, bias_columns) {
    data.frame(
      parameter = parameters,
      # Regular empirical bias
      bias = vapply(fits[, bias_columns, drop = FALSE], mean, numeric(1)),
      # Empirical mean of absolute bias
      absolute_bias = vapply(fits[, bias_columns, drop = FALSE], function(bias) {
        mean(abs(bias))
      }, numeric(1)),
      # RMSE of biases
      rmse = vapply(fits[, bias_columns, drop = FALSE], function(bias) {
        sqrt(mean(bias^2))
      }, numeric(1)),
      # Standard error of biases
      bias_mc_se = vapply(fits[, bias_columns, drop = FALSE], function(bias) {
        if (length(bias) < 2L) {
          return(NA_real_)
        }
        stats::sd(bias) / sqrt(length(bias))
      }, numeric(1))
    )
  }

  empirical_beta_summary <- summarize_empirical_bias(
    successful_fits,
    c("beta_1", "beta_2", "beta_3"),
    beta_bias_columns
  )
  empirical_gamma_summary <- summarize_empirical_bias(
    successful_fits,
    c("gamma_1", "gamma_2"),
    gamma_bias_columns
  )
  convergence_summary <- data.frame(
    n_data_sets = n_data_sets,
    n_skipped = sum(study_summary$skipped),
    n_fit_attempted = n_fitted_rows,
    n_successful = n_successful_fits,
    convergence_rate = n_successful_fits / n_fitted_rows
  )
  print(trial_results)
  print(empirical_beta_summary)
  print(empirical_gamma_summary)
  #print(convergence_summary)

  expect_equal(length(unique(study_summary$data_set_index)), n_data_sets)
  expect_gt(nrow(successful_fits), 0)
  expect_equal(
    nrow(unique(study_summary[, c(
      "true_beta_1",
      "true_beta_2",
      "true_beta_3",
      "true_gamma_1",
      "true_gamma_2"
    ), drop = FALSE])),
    1
  )
  expect_true(all(study_summary$true_beta_1 == true_beta[1]))
  expect_true(all(study_summary$true_beta_2 == true_beta[2]))
  expect_true(all(study_summary$true_beta_3 == true_beta[3]))
  expect_true(all(study_summary$true_gamma_1 == true_gamma[1]))
  expect_true(all(study_summary$true_gamma_2 == true_gamma[2]))
  expect_true(all(fitted_rows$gamma_start_1 == gamma_start_values[1]))
  expect_true(all(fitted_rows$gamma_start_2 == gamma_start_values[2]))
  expect_true(all(table(successful_fits$data_set_index) == 1))
  if (nrow(fitted_rows) > 0) {
    expect_lte(max(fitted_rows$start_index), max_start_attempts)
  }
  expect_true(all(successful_fits$treated_above >= min_group_size))
  expect_true(all(successful_fits$treated_below >= min_group_size))
  expect_true(all(successful_fits$control_above >= min_group_size))
  expect_true(all(successful_fits$control_below >= min_group_size))
  expect_true(all(is.finite(successful_fits$gamma_1_value)))
  expect_true(all(is.finite(successful_fits$gamma_2_value)))
  expect_true(all(is.finite(successful_fits$beta_1_value)))
  expect_true(all(is.finite(successful_fits$beta_2_value)))
  expect_true(all(is.finite(successful_fits$beta_3_value)))
  expect_true(all(is.finite(successful_fits$beta_1_bias)))
  expect_true(all(is.finite(successful_fits$beta_2_bias)))
  expect_true(all(is.finite(successful_fits$beta_3_bias)))
  expect_true(all(is.finite(successful_fits$gamma_1_bias)))
  expect_true(all(is.finite(successful_fits$gamma_2_bias)))
  expect_true(all(is.finite(successful_fits$beta_1_absolute_bias)))
  expect_true(all(is.finite(successful_fits$beta_2_absolute_bias)))
  expect_true(all(is.finite(successful_fits$beta_3_absolute_bias)))
  expect_true(all(is.finite(successful_fits$gamma_1_absolute_bias)))
  expect_true(all(is.finite(successful_fits$gamma_2_absolute_bias)))
  expect_true(all(is.finite(successful_fits$beta_1_rmse)))
  expect_true(all(is.finite(successful_fits$beta_2_rmse)))
  expect_true(all(is.finite(successful_fits$beta_3_rmse)))
  expect_true(all(is.finite(successful_fits$gamma_1_rmse)))
  expect_true(all(is.finite(successful_fits$gamma_2_rmse)))
  expect_true(all(is.finite(successful_fits$beta_1_squared_bias)))
  expect_true(all(is.finite(successful_fits$beta_2_squared_bias)))
  expect_true(all(is.finite(successful_fits$beta_3_squared_bias)))
  expect_true(all(is.finite(successful_fits$gamma_1_squared_bias)))
  expect_true(all(is.finite(successful_fits$gamma_2_squared_bias)))
  expect_true(all(is.finite(empirical_beta_summary$bias)))
  expect_true(all(is.finite(empirical_beta_summary$absolute_bias)))
  expect_true(all(is.finite(empirical_beta_summary$rmse)))
  expect_true(all(is.finite(empirical_beta_summary$bias_mc_se)) || n_successful_fits < 2L)
  expect_true(all(is.finite(empirical_gamma_summary$bias)))
  expect_true(all(is.finite(empirical_gamma_summary$absolute_bias)))
  expect_true(all(is.finite(empirical_gamma_summary$rmse)))
  expect_true(all(is.finite(empirical_gamma_summary$bias_mc_se)) || n_successful_fits < 2L)
  expect_true(all(is.finite(convergence_summary$convergence_rate)))
})
