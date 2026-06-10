test_that("simulate_threshold_data returns a complete survival data frame", {
  set.seed(4001)
  data <- simulate_threshold_data(
    n = 50,
    beta = c(log(1.5), log(1.2), log(1.4)),
    gamma = c(0.5, -0.25, 0.1),
    baseline_hazard = 1.1,
    biomarker_correlation = 0.2,
    study_end_range = c(0.4, 2)
  )

  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 50)
  expect_named(data, c("time", "status", "treatment", "bio_1", "bio_2", "bio_3"))
  expect_true(all(data$time > 0))
  expect_true(all(data$status %in% c(0L, 1L)))
  expect_true(all(data$treatment %in% c(0L, 1L)))
  expect_true(all(is.finite(as.matrix(data[, c("bio_1", "bio_2", "bio_3")]))))
})

test_that("simulate_threshold_data validates beta and gamma", {
  expect_error(
    simulate_threshold_data(10, beta = c(1, 2), gamma = c(0.1, 0.2)),
    "beta must contain"
  )
  expect_error(
    simulate_threshold_data(10, beta = c(1, 2, 3), gamma = c(0.1, NA)),
    "finite and numeric"
  )
})

test_that("simulate_threshold_data accepts a one-biomarker model", {
  set.seed(4002)
  data <- simulate_threshold_data(
    n = 20,
    beta = c(0.1, 0.2, 0.3),
    gamma = 0.8,
    biomarker_correlation = 0,
    study_end_range = c(0.4, 2)
  )

  expect_named(data, c("time", "status", "treatment", "bio_1"))
  expect_equal(ncol(data), 4)
})

test_that("summarize_simulated_data reports score, censoring, and subgroup summaries", {
  data <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(1, 0, 1, 0),
    treatment = c(1, 1, 0, 0),
    bio_1 = c(2, -3, 1, -2)
  )

  summary <- summarize_simulated_data(data, gamma = 1)

  expect_named(summary, c(
    "score_quantiles",
    "overall_censoring",
    "group_sizes",
    "group_censoring",
    "score_side_counts",
    "threshold",
    "score_shift"
  ))
  expect_equal(summary$overall_censoring, 0.5)
  expect_equal(sum(summary$group_sizes), nrow(data))
  expect_equal(sum(summary$score_side_counts), nrow(data))
  expect_equal(summary$threshold, 0)
  expect_equal(summary$score_shift, 1)
})

test_that("summarize_simulated_data supports custom column names", {
  data <- data.frame(
    follow_up = c(1, 2, 3, 4),
    event = c(1, 0, 1, 0),
    arm = c(1, 1, 0, 0),
    marker_a = c(2, -3, 1, -2),
    marker_b = c(0, 1, -1, 2)
  )

  summary <- summarize_simulated_data(
    data,
    gamma = c(1, -0.5),
    treatment_col = "arm",
    status_col = "event",
    biomarker_cols = c("marker_a", "marker_b")
  )

  expect_equal(sum(summary$group_sizes), nrow(data))
  expect_true(all(names(summary$group_sizes) %in% c(
    "treated_above",
    "treated_below",
    "control_above",
    "control_below"
  )))
})
