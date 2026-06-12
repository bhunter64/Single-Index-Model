test_that("threshold model can be fit repeatedly from different gamma starts", {
  set.seed(20260526)

  true_gamma <- c(1.5, 1.8)
  true_beta <- c(log(3), log(1.2), log(2))

  data <- simulate_threshold_data(
    n = 600,
    beta = true_beta,
    gamma = true_gamma,
    baseline_hazard = 1,
    biomarker_correlation = 0,
    study_end_range = c(0.4, 3)
  )

  #Km test on the simulated data using the true known gamma parameters
  print(summarize_simulated_data(data, true_gamma))
  km_test <- kaplan_meier_threshold_test(
    data = data,
    gamma = true_gamma,
    strata = "treatment_threshold"
  )

  # Print just the Chi-Squared and p-value from the above km test
  print(data.frame(
    km_logrank_chisq = km_test$logrank$chisq,
    km_logrank_p = km_test$p_value
  ))

  # Set up simulation metrics
  control <- default_mcmc_control(
    samples = 20000,
    burn_in = 10000,
    thin = 10,
    gamma_mean = c(0.5, 0.5),
    gamma_sd = 1,
    gamma_proposal_sd = 0.2
  )

  # Correlation helper function
  fit_rows <- list()
  sample_correlation <- function(x, y) {
    if (stats::sd(x) == 0 || stats::sd(y) == 0) {
      return(NA_real_)
    }
    stats::cor(x, y)
  }

  #' Randomly sample starting gamma values until one fit converges.
  #' If it converges, print the fitted values and converged metrics.
  #' If it fails to converge, print NA values and `converged = FALSE`
  max_start_attempts <- 100
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
        gof <- threshold_fit_gof_tests(data, summary$gamma)

        cbind(data.frame(
          start_index = start_index,
          gamma_start_1 = gamma_start[1],
          gamma_start_2 = gamma_start[2],
          converged = TRUE,
          error = NA_character_,
          beta_1 = summary$beta[1],
          beta_2 = summary$beta[2],
          beta_3 = summary$beta[3],
          gamma_1 = summary$gamma[1],
          gamma_2 = summary$gamma[2],
          mean_loglik = mean(posterior$loglik_samples),
          gamma_correlation = sample_correlation(
            posterior$gamma_samples[1, ],
            posterior$gamma_samples[2, ]
          )
        ), gof)
      },
      error = function(err) {
        data.frame(
          start_index = start_index,
          gamma_start_1 = gamma_start[1],
          gamma_start_2 = gamma_start[2],
          converged = FALSE,
          error = conditionMessage(err),
          beta_1 = NA_real_,
          beta_2 = NA_real_,
          beta_3 = NA_real_,
          gamma_1 = NA_real_,
          gamma_2 = NA_real_,
          mean_loglik = NA_real_,
          gamma_correlation = NA_real_,
          gof_converged = FALSE,
          loglik = NA_real_,
          aic = NA_real_,
          concordance = NA_real_,
          likelihood_ratio_p = NA_real_,
          wald_p = NA_real_,
          score_p = NA_real_,
          ph_global_p = NA_real_,
          cox_snell_ks_statistic = NA_real_,
          cox_snell_ks_p = NA_real_
        )
      }
    )

    if (isTRUE(fit_rows[[start_index]]$converged)) {
      break
    }
  }

  # Print all resulting summaries
  fit_summary <- do.call(rbind, fit_rows)
  print(fit_summary)

  # Find all converged fits
  successful_fits <- fit_summary[fit_summary$converged, , drop = FALSE]

  # Logical testing metrics to pass each time
  expect_gt(nrow(successful_fits), 0)
  expect_equal(nrow(successful_fits), 1)
  expect_lte(nrow(fit_summary), max_start_attempts)
  expect_true(all(is.finite(successful_fits$gamma_1)))
  expect_true(all(is.finite(successful_fits$gamma_2)))
  expect_true(all(is.finite(successful_fits$beta_1)))
  expect_true(all(is.finite(successful_fits$beta_2)))
  expect_true(all(is.finite(successful_fits$beta_3)))
  expect_true(all(is.finite(successful_fits$mean_loglik)))
  expect_true(all(successful_fits$gof_converged))
  expect_true(all(is.finite(successful_fits$aic)))
  expect_true(all(is.finite(successful_fits$concordance)))
  expect_true(all(is.finite(successful_fits$cox_snell_ks_p)))
})
