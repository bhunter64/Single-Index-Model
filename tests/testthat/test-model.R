#' Run the practicum model tests.
#'
#' This is the single executable test entry point for the cleaned direct-gamma
#' workflow. The test simulates a small survival data set with threshold fixed
#' at zero, runs a short MCMC chain, summarizes the posterior, and checks that
#' core dimensions and values are valid.
#'
#' @return Invisibly returns the posterior summary.
run_practicum_tests <- function() {
  set.seed(20260526)

  true_gamma <- c(0.75, 0.25)
  true_beta <- c(log(3), log(1.2), log(2))

  data <- simulate_threshold_data(
    n = 400,
    beta = true_beta,
    gamma = true_gamma,
    baseline_hazard = 1
  )

  control <- default_mcmc_control(
    samples = 500,
    burn_in = 100,
    thin = 1,
    gamma_mean = c(0.7, 0.3),
    gamma_sd = 0.1,
    gamma_proposal_sd = 0.04
  )

  posterior <- fit_threshold_model(
    data = data,
    control = control,
    lambda = 0
  )
  summary <- summarize_mcmc(posterior)

  stopifnot(
    nrow(posterior$gamma_samples) == length(true_gamma),
    length(summary$beta) == length(true_beta),
    all(is.finite(summary$gamma)),
    all(is.finite(summary$beta)),
    identical(summary$threshold, 0),
    identical(summary$score_shift, 1)
  )

  print(summary)
  invisible(summary)
}

run_practicum_tests()
