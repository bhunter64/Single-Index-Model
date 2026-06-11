test_that("kaplan_meier_threshold_test returns a log-rank comparison", {
  set.seed(20260610)

  data <- simulate_threshold_data(
    n = 80,
    beta = c(log(1.5), log(1.2), log(1.4)),
    gamma = c(0.8, -0.4),
    baseline_hazard = 0.8,
    biomarker_correlation = 0,
    study_end_range = c(0.5, 2)
  )

  km_test <- kaplan_meier_threshold_test(data, gamma = c(0.8, -0.4))

  expect_s3_class(km_test$km_fit, "survfit")
  expect_s3_class(km_test$logrank, "survdiff")
  expect_true(is.finite(km_test$p_value))
  expect_equal(sum(km_test$group_counts), nrow(data))
})

test_that("threshold_fit_gof_tests returns one row of diagnostics", {
  set.seed(20260610)

  data <- simulate_threshold_data(
    n = 100,
    beta = c(log(1.5), log(1.2), log(1.4)),
    gamma = c(0.8, -0.4),
    baseline_hazard = 0.8,
    biomarker_correlation = 0,
    study_end_range = c(0.5, 2)
  )

  gof <- threshold_fit_gof_tests(data, gamma = c(0.8, -0.4))

  expect_equal(nrow(gof), 1)
  expect_true(gof$gof_converged)
  expect_true(is.finite(gof$aic))
  expect_true(is.finite(gof$concordance))
  expect_true(is.finite(gof$cox_snell_ks_p))
})
