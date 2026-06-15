test_that("threshold model converges across randomly generated simulation study data sets", {
  set.seed(20260612)

  n_data_sets <- 20
  min_group_size <- 60
  max_start_attempts <- 100
  n_biomarkers <- 2

  control <- default_mcmc_control(
    samples = 5000,
    burn_in = 6000,
    thin = 4,
    gamma_mean = c(0.5, 0.5),
    gamma_sd = 1,
    gamma_proposal_sd = 0.2
  )

  sample_true_gamma <- function() {
    stats::runif(n_biomarkers, min = 1.2, max = 2.2) *
      sample(c(-1, 1), n_biomarkers, replace = TRUE)
  }

  sample_true_beta <- function() {
    log(c(
      stats::runif(1, min = 1.5, max = 3.5),
      stats::runif(1, min = 1.05, max = 1.6),
      stats::runif(1, min = 1.4, max = 2.6)
    ))
  }

  skipped_data_set_row <- function(data_set_index, true_beta, true_gamma, data_summary, skip_reason) {
    data.frame(
      data_set_index = data_set_index,
      start_index = NA_integer_,
      true_beta_1 = true_beta[1],
      true_beta_2 = true_beta[2],
      true_beta_3 = true_beta[3],
      true_gamma_1 = true_gamma[1],
      true_gamma_2 = true_gamma[2],
      gamma_start_1 = NA_real_,
      gamma_start_2 = NA_real_,
      treated_above = data_summary$group_sizes[["treated_above"]],
      treated_below = data_summary$group_sizes[["treated_below"]],
      control_above = data_summary$group_sizes[["control_above"]],
      control_below = data_summary$group_sizes[["control_below"]],
      skipped = TRUE,
      skip_reason = skip_reason,
      converged = FALSE,
      error = NA_character_,
      beta_1 = NA_real_,
      beta_2 = NA_real_,
      beta_3 = NA_real_,
      beta_bias_1 = NA_real_,
      beta_bias_2 = NA_real_,
      beta_bias_3 = NA_real_,
      gamma_1 = NA_real_,
      gamma_2 = NA_real_
    )
  }

  simulate_data_set <- function(data_set_index) {
    true_gamma <- sample_true_gamma()
    true_beta <- sample_true_beta()

    data <- simulate_threshold_data(
      n = 600,
      beta = true_beta,
      gamma = true_gamma,
      baseline_hazard = 1,
      biomarker_correlation = 0,
      study_end_range = c(0.4, 3)
    )
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
      gamma_start <- draw_gamma(
        length(true_gamma),
        mean = control$gamma_mean,
        sd = control$gamma_sd
      )

      fit_rows[[start_index]] <- tryCatch(
        {
          posterior <- suppressWarnings(fit_threshold_model(
            data = data,
            control = control,
            lambda = 0,
            gamma_start = gamma_start
          ))
          summary <- summarize_mcmc(posterior)
          beta_bias <- summary$beta - true_beta

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
            skipped = FALSE,
            skip_reason = NA_character_,
            converged = TRUE,
            error = NA_character_,
            beta_1 = summary$beta[1],
            beta_2 = summary$beta[2],
            beta_3 = summary$beta[3],
            beta_bias_1 = beta_bias[1],
            beta_bias_2 = beta_bias[2],
            beta_bias_3 = beta_bias[3],
            gamma_1 = summary$gamma[1],
            gamma_2 = summary$gamma[2]
          )
        },
        error = function(err) {
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
            skipped = FALSE,
            skip_reason = NA_character_,
            converged = FALSE,
            error = conditionMessage(err),
            beta_1 = NA_real_,
            beta_2 = NA_real_,
            beta_3 = NA_real_,
            beta_bias_1 = NA_real_,
            beta_bias_2 = NA_real_,
            beta_bias_3 = NA_real_,
            gamma_1 = NA_real_,
            gamma_2 = NA_real_
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

  study_summary <- do.call(rbind, study_rows)
  successful_fits <- study_summary[study_summary$converged, , drop = FALSE]
  successful_beta_bias <- successful_fits[, c(
    "data_set_index",
    "true_beta_1",
    "true_beta_2",
    "true_beta_3",
    "beta_1",
    "beta_2",
    "beta_3",
    "beta_bias_1",
    "beta_bias_2",
    "beta_bias_3"
  ), drop = FALSE]
  empirical_beta_bias <- data.frame(
    beta_1 = mean(successful_fits$beta_bias_1),
    beta_2 = mean(successful_fits$beta_bias_2),
    beta_3 = mean(successful_fits$beta_bias_3)
  )
  print(successful_beta_bias)
  print(empirical_beta_bias)

  fitted_rows <- study_summary[!study_summary$skipped, , drop = FALSE]

  expect_equal(length(unique(study_summary$data_set_index)), n_data_sets)
  expect_gt(nrow(successful_fits), 0)
  expect_true(all(table(successful_fits$data_set_index) == 1))
  if (nrow(fitted_rows) > 0) {
    expect_lte(max(fitted_rows$start_index), max_start_attempts)
  }
  expect_true(all(successful_fits$treated_above >= min_group_size))
  expect_true(all(successful_fits$treated_below >= min_group_size))
  expect_true(all(successful_fits$control_above >= min_group_size))
  expect_true(all(successful_fits$control_below >= min_group_size))
  expect_true(all(is.finite(successful_fits$gamma_1)))
  expect_true(all(is.finite(successful_fits$gamma_2)))
  expect_true(all(is.finite(successful_fits$beta_1)))
  expect_true(all(is.finite(successful_fits$beta_2)))
  expect_true(all(is.finite(successful_fits$beta_3)))
  expect_true(all(is.finite(successful_fits$beta_bias_1)))
  expect_true(all(is.finite(successful_fits$beta_bias_2)))
  expect_true(all(is.finite(successful_fits$beta_bias_3)))
  expect_true(all(is.finite(empirical_beta_bias$beta_1)))
  expect_true(all(is.finite(empirical_beta_bias$beta_2)))
  expect_true(all(is.finite(empirical_beta_bias$beta_3)))
})
