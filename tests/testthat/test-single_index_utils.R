test_that("gamma validation accepts finite numeric vectors", {
  expect_equal(validate_gamma(c(0.1, -0.2), 2), c(0.1, -0.2))
  expect_type(validate_gamma(1:3, 3), "double")
})

test_that("gamma validation rejects invalid dimensions and values", {
  expect_error(validate_gamma(c(1, 2), 3), "one coefficient per biomarker")
  expect_error(validate_gamma(c(1, NA), 2), "finite and numeric")
  expect_error(validate_gamma(c(1, Inf), 2), "finite and numeric")
  expect_error(suppressWarnings(validate_gamma(c("a", "b"), 2)), "finite and numeric")
})

test_that("biomarker_score computes the shifted linear score", {
  biomarkers <- data.frame(bio_1 = c(0, 2, -1), bio_2 = c(1, -1, 3))
  gamma <- c(0.5, -0.25)

  expect_equal(
    biomarker_score(biomarkers, gamma),
    as.vector(1 + as.matrix(biomarkers) %*% gamma)
  )
})

test_that("biomarker_score validates gamma length", {
  biomarkers <- data.frame(bio_1 = 1:3, bio_2 = 4:6)

  expect_error(biomarker_score(biomarkers, 0.2), "one coefficient per biomarker")
})

test_that("resolve_biomarker_columns handles defaults, names, and indices", {
  x <- data.frame(treatment = c(0, 1), bio_1 = c(1, 2), bio_2 = c(3, 4))

  expect_equal(resolve_biomarker_columns(x, "treatment"), c(2L, 3L))
  expect_equal(resolve_biomarker_columns(x, 1, c("bio_2", "bio_1")), c(3L, 2L))
  expect_equal(resolve_biomarker_columns(x, 1, c(3, 2)), c(3L, 2L))
})

test_that("resolve_biomarker_columns rejects missing or invalid columns", {
  x <- data.frame(treatment = c(0, 1), bio_1 = c(1, 2))

  expect_error(resolve_biomarker_columns(x, "treatment", integer()), "at least one")
  expect_error(resolve_biomarker_columns(x, "treatment", "missing"), "unknown column")
  expect_error(resolve_biomarker_columns(x, "treatment", 4), "outside the data")
})

test_that("make_threshold_design builds treatment, subgroup, and interaction columns", {
  x <- data.frame(treatment = c(0, 1, 1, 0), bio_1 = c(-3, -1, 1, 3))

  design <- make_threshold_design(x, gamma = 1, treatment_col = "treatment")

  expect_named(design, c("treatment", "subgroup", "treatment_subgroup"))
  expect_equal(design$treatment, x$treatment)
  expect_equal(design$subgroup, c(0L, 1L, 1L, 1L))
  expect_equal(design$treatment_subgroup, design$treatment * design$subgroup)
})

test_that("gamma prior and proposal helpers have expected dimensions", {
  set.seed(1001)
  gamma <- draw_gamma(3, mean = c(0, 1, 2), sd = 0.1)
  expect_length(gamma, 3)
  expect_true(all(is.finite(gamma)))

  expect_equal(
    log_gamma_prior(c(0, 1), mean = c(0, 1), sd = c(1, 2)),
    sum(stats::dnorm(c(0, 1), mean = c(0, 1), sd = c(1, 2), log = TRUE))
  )

  proposal <- propose_gamma(c(1, 2, 3), proposal_sd = c(0, 0, 0))
  expect_equal(proposal, c(1, 2, 3))
})

test_that("distribution helper lengths are checked", {
  expect_equal(singleIndexModel:::repeat_parameter(2, 3), c(2, 2, 2))
  expect_equal(singleIndexModel:::repeat_parameter(c(1, 2), 2), c(1, 2))
  expect_error(
    singleIndexModel:::repeat_parameter(c(1, 2), 3),
    "length must be 1 or match"
  )
})
