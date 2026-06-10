test_that("cox_threshold_loglik matches a direct partial-likelihood calculation", {
  x <- data.frame(
    treatment = c(1, 0, 1, 0),
    bio_1 = c(2, -3, 0, 1)
  )
  y <- survival::Surv(time = c(4, 3, 2, 1), event = c(1, 0, 1, 1))
  beta <- c(0.2, -0.1, 0.5)
  gamma <- 1

  design <- make_threshold_design(x, gamma)
  linear_predictor <- as.vector(as.matrix(design) %*% beta)
  expected <- sum(y[, 2] * (linear_predictor - log(cumsum(exp(linear_predictor)))))

  expect_equal(cox_threshold_loglik(x, y, beta, gamma), expected)
})

test_that("cox_threshold_loglik supports named treatment and biomarker columns", {
  x <- data.frame(arm = c(1, 0, 1), marker_a = c(1, -2, 3), marker_b = c(0, 1, -1))
  y <- survival::Surv(time = c(3, 2, 1), event = c(1, 1, 0))

  expect_true(is.finite(cox_threshold_loglik(
    x,
    y,
    beta = c(0.1, 0.2, -0.3),
    gamma = c(0.5, -0.25),
    treatment_col = "arm",
    biomarker_cols = c("marker_a", "marker_b")
  )))
})

test_that("fit_threshold_cox returns a converged Cox model on simulated data", {
  set.seed(2001)
  data <- simulate_threshold_data(
    n = 120,
    beta = c(log(1.5), log(1.1), log(1.4)),
    gamma = c(0.8, -0.4),
    baseline_hazard = 1.2,
    biomarker_correlation = 0.1,
    study_end_range = c(0.5, 3)
  )
  sorted <- sort(data$time, decreasing = TRUE, index.return = TRUE)
  y <- survival::Surv(data$time, data$status)[sorted$ix, ]
  x <- data[sorted$ix, c("treatment", "bio_1", "bio_2")]

  fit <- fit_threshold_cox(x, y, gamma = c(0.8, -0.4))

  expect_s3_class(fit, "coxph")
  expect_true(isTRUE(fit$converged))
  expect_length(stats::coef(fit), 3)
  expect_true(all(is.finite(stats::coef(fit))))
})

test_that("initialize_beta returns the Cox coefficients and validates gamma", {
  set.seed(2002)
  data <- simulate_threshold_data(
    n = 100,
    beta = c(log(1.4), log(1.2), log(1.3)),
    gamma = c(0.7, 0.5),
    study_end_range = c(0.5, 3)
  )
  sorted <- sort(data$time, decreasing = TRUE, index.return = TRUE)
  y <- survival::Surv(data$time, data$status)[sorted$ix, ]
  x <- data[sorted$ix, c("treatment", "bio_1", "bio_2")]

  beta <- initialize_beta(x, y, gamma = c(0.7, 0.5))

  expect_length(beta, 3)
  expect_true(all(is.finite(beta)))
  expect_error(initialize_beta(x, y, gamma = 0.7), "one coefficient per biomarker")
})
