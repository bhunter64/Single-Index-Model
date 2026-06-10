test_that("default_mcmc_control records sampler settings", {
  control <- default_mcmc_control(
    samples = 12,
    burn_in = 4,
    thin = 3,
    gamma_mean = c(0, 1),
    gamma_sd = c(0.2, 0.3),
    gamma_proposal_sd = c(0.01, 0.02)
  )

  expect_equal(control$samples, 12)
  expect_equal(control$burn_in, 4)
  expect_equal(control$thin, 3)
  expect_equal(control$gamma_mean, c(0, 1))
  expect_equal(control$gamma_sd, c(0.2, 0.3))
  expect_equal(control$gamma_proposal_sd, c(0.01, 0.02))
})

test_that("mcmc_step returns a valid updated sampler state", {
  set.seed(1)
  data <- simulate_threshold_data(
    n = 180,
    beta = c(log(1.4), log(1.2), log(1.3)),
    gamma = c(0.4, 0.6),
    baseline_hazard = 0.5,
    study_end_range = c(0.5, 3)
  )
  sorted <- sort(data$time, decreasing = TRUE, index.return = TRUE)
  x <- data[sorted$ix, c("treatment", "bio_1", "bio_2")]
  y <- survival::Surv(data$time, data$status)[sorted$ix, ]
  gamma <- c(0.4, 0.6)
  beta <- initialize_beta(x, y, gamma)

  state <- mcmc_step(
    x = x,
    y = y,
    beta = beta,
    gamma = gamma,
    control = default_mcmc_control(
      gamma_mean = c(0, 0),
      gamma_sd = 1,
      gamma_proposal_sd = c(0, 0)
    )
  )

  expect_length(state$beta, 3)
  expect_length(state$gamma, 2)
  expect_true(all(is.finite(state$beta)))
  expect_true(all(is.finite(state$gamma)))
  expect_true(is.finite(state$log_likelihood))
})

test_that("fit_threshold_mcmc returns retained beta, gamma, and log-likelihood samples", {
  set.seed(1)
  data <- simulate_threshold_data(
    n = 180,
    beta = c(log(1.4), log(1.2), log(1.3)),
    gamma = c(0.4, 0.6),
    baseline_hazard = 0.5,
    study_end_range = c(0.5, 3)
  )
  sorted <- sort(data$time, decreasing = TRUE, index.return = TRUE)
  x <- data[sorted$ix, c("treatment", "bio_1", "bio_2")]
  y <- survival::Surv(data$time, data$status)[sorted$ix, ]
  control <- default_mcmc_control(
    samples = 5,
    burn_in = 2,
    thin = 2,
    gamma_mean = c(0, 0),
    gamma_sd = 1,
    gamma_proposal_sd = 0
  )

  posterior <- fit_threshold_mcmc(x, y, control = control, gamma_start = c(0.4, 0.6))

  expect_equal(dim(posterior$beta_samples), c(3L, 5L))
  expect_equal(dim(posterior$gamma_samples), c(2L, 5L))
  expect_equal(dim(posterior$loglik_samples), c(1L, 5L))
  expect_true(all(is.finite(posterior$beta_samples)))
  expect_true(all(is.finite(posterior$gamma_samples)))
  expect_true(all(is.finite(posterior$loglik_samples)))
  expect_equal(posterior$threshold, 0)
  expect_equal(posterior$score_shift, 1)
})

test_that("fit_threshold_mcmc rejects invalid inputs", {
  set.seed(3003)
  data <- simulate_threshold_data(
    n = 80,
    beta = c(log(1.4), log(1.2), log(1.3)),
    gamma = c(0.4, 0.6),
    study_end_range = c(0.5, 3)
  )
  sorted <- sort(data$time, decreasing = TRUE, index.return = TRUE)
  x <- data[sorted$ix, c("treatment", "bio_1", "bio_2")]
  y <- survival::Surv(data$time, data$status)[sorted$ix, ]
  control <- default_mcmc_control(samples = 1, burn_in = 0, thin = 1)

  expect_error(fit_threshold_mcmc(x, y, control = control, lambda = -1), "nonnegative")
  expect_error(
    fit_threshold_mcmc(x, y, control = control, gamma_start = 0.1),
    "one coefficient per biomarker"
  )
})

test_that("summarize_mcmc supports mean and median summaries", {
  posterior <- list(
    beta_samples = matrix(c(1, 2, 3, 4, 5, 6), nrow = 3),
    gamma_samples = matrix(c(-1, 1, 0, 2), nrow = 2),
    threshold = 0,
    score_shift = 1
  )

  mean_summary <- summarize_mcmc(posterior, method = "mean")
  median_summary <- summarize_mcmc(posterior, method = "median")

  expect_equal(mean_summary$beta, rowMeans(posterior$beta_samples))
  expect_equal(median_summary$gamma, apply(posterior$gamma_samples, 1, stats::median))
  expect_equal(dim(mean_summary$beta_interval), c(2L, 3L))
  expect_equal(dim(mean_summary$gamma_interval), c(2L, 2L))
  expect_error(summarize_mcmc(posterior, method = "mode"), "mean' or 'median")
})

test_that("metropolis acceptance rejects non-finite ratios", {
  expect_false(singleIndexModel:::accept_metropolis(NA_real_, 0))
  expect_false(singleIndexModel:::accept_metropolis(Inf, Inf))

  withr_seed <- function(seed, expr) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
    on.exit({
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    })
    set.seed(seed)
    force(expr)
  }

  expect_true(withr_seed(1, singleIndexModel:::accept_metropolis(1, 0)))
})
