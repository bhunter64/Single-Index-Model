test_that("threshold model can be fit repeatedly from different gamma starts", {
  set.seed(20260526)

  true_gamma <- c(1.5, 1.8)
  true_beta <- c(log(3), log(1.2), log(2))

  data <- simulate_threshold_data(
    n = 160,
    beta = true_beta,
    gamma = true_gamma,
    baseline_hazard = 1,
    biomarker_correlation = 0,
    study_end_range = c(0.4, 3)
  )

  print(summarize_simulated_data(data, true_gamma))

  starting_values <- expand.grid(
    gamma_1 = seq(from = 0, to = 1, by = 0.5),
    gamma_2 = seq(from = 0, to = 1, by = 0.5)
  )
  starting_values <- rbind(
    starting_values,
    data.frame(gamma_1 = true_gamma[1], gamma_2 = true_gamma[2])
  )

  control <- default_mcmc_control(
    samples = 20,
    burn_in = 20,
    thin = 1,
    gamma_mean = c(0, 0),
    gamma_sd = 1,
    gamma_proposal_sd = 0.2
  )

  fit_rows <- vector("list", nrow(starting_values))

  for (start_index in seq_len(nrow(starting_values))) {
    gamma_start <- as.numeric(starting_values[start_index, ])

    fit_rows[[start_index]] <- tryCatch(
      {
        posterior <- fit_threshold_model(
          data = data,
          control = control,
          lambda = 0,
          gamma_start = gamma_start
        )
        summary <- summarize_mcmc(posterior)

        data.frame(
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
          gamma_correlation = stats::cor(
            posterior$gamma_samples[1, ],
            posterior$gamma_samples[2, ]
          )
        )
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
          gamma_correlation = NA_real_
        )
      }
    )
  }

  fit_summary <- do.call(rbind, fit_rows)
  print(fit_summary)

  successful_fits <- fit_summary[fit_summary$converged, , drop = FALSE]

  expect_gt(nrow(successful_fits), 0)
  expect_equal(nrow(fit_summary), nrow(starting_values))
  expect_true(all(is.finite(successful_fits$gamma_1)))
  expect_true(all(is.finite(successful_fits$gamma_2)))
  expect_true(all(is.finite(successful_fits$beta_1)))
  expect_true(all(is.finite(successful_fits$beta_2)))
  expect_true(all(is.finite(successful_fits$beta_3)))
  expect_true(all(is.finite(successful_fits$mean_loglik)))
})
