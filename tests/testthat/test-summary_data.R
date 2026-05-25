test_that("summary is factual", {
  dat <-simu_data(400, beta = c(log(3), log(2), log(0.5)), gamma = c(-1.8, -2, 0.9))
  summary_data(dat)
})
