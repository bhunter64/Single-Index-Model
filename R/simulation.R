#' Simulate survival data for the model.
#'
#' @param n Number of observations.
#' @param beta A vector of regression coefficients for treatment,
#'   biomarker subgroup, and their interaction term.
#' @param gamma A vector of biomarker coefficients.
#' @param baseline_hazard Baseline exponential hazard rate.
#' @param biomarker_correlation Pairwise correlation used for simulated
#'   biomarkers.
#' @param study_end_range Two positive values giving the uniform censoring
#'   interval.
#'
#' @return A data frame with columns `time`, `status`, `treatment`, and one column per
#'   biomarker.
#' @export
simulate_threshold_data <- function(n,
                                    beta,
                                    gamma,
                                    baseline_hazard = 1,
                                    biomarker_correlation = 0.5,
                                    study_end_range = c(0.05, 3)) {
  gamma <- validate_gamma(gamma, length(gamma))
  beta <- as.numeric(beta)
  if (length(beta) != 3L) {
    stop("beta must contain treatment, subgroup, and interaction coefficients", call. = FALSE)
  }

  treatment <- stats::rbinom(n, 1, 0.5)
  n_biomarkers <- length(gamma)
  biomarker_mean <- rep(0, n_biomarkers)
  biomarker_covariance <- matrix(biomarker_correlation, nrow = n_biomarkers, ncol = n_biomarkers)
  diag(biomarker_covariance) <- 1
  biomarkers <- mvtnorm::rmvnorm(n, mean = biomarker_mean, sigma = biomarker_covariance)
  colnames(biomarkers) <- paste0("bio_", seq_len(n_biomarkers))

  x <- data.frame(treatment = treatment, biomarkers, check.names = FALSE)
  design <- make_threshold_design(x, gamma)
  hazard <- baseline_hazard * exp(as.matrix(design) %*% beta)

  event_time <- stats::rexp(n, hazard)
  study_end <- stats::runif(n, study_end_range[1], study_end_range[2])
  status <- as.integer(event_time <= study_end)
  observed_time <- pmin(event_time, study_end)

  return(data.frame(time = observed_time, status = status, treatment = treatment, biomarkers, check.names = FALSE))
}

#' Summarize simulated data by fixed threshold subgroup.
#'
#' @param data A simulated threshold-model data frame.
#' @param gamma A vector of biomarker coefficients.
#' @param treatment_col Column name for treatment.
#' @param status_col Column name for event status.
#' @param biomarker_cols Optional column names or indices for biomarkers.
#'
#' @return A list containing biomarker-score quantiles, overall censoring,
#'   subgroup sizes, subgroup censoring proportions, and above/below counts.
#' @export
summarize_simulated_data <- function(data,
                                     gamma,
                                     treatment_col = "treatment",
                                     status_col = "status",
                                     biomarker_cols = NULL) {
  data <- as.data.frame(data)
  if (is.null(biomarker_cols)) {
    biomarker_cols <- setdiff(seq_along(data), match(c("time", status_col, treatment_col), names(data)))
  } else {
    biomarker_cols <- match_columns(data, biomarker_cols)
  }

  score <- biomarker_score(data[, biomarker_cols, drop = FALSE], gamma)
  above_threshold <- score >= 0
  treated <- data[[treatment_col]] == 1
  groups <- list(
    treated_above = treated & above_threshold,
    treated_below = treated & !above_threshold,
    control_above = !treated & above_threshold,
    control_below = !treated & !above_threshold
  )

  group_sizes <- vapply(groups, sum, integer(1))
  censoring_by_group <- vapply(groups, function(idx) {
    if (!any(idx)) {
      return(NA_real_)
    }
    1 - mean(data[[status_col]][idx])
  }, numeric(1))

  return(list(
    score_quantiles = stats::quantile(score),
    overall_censoring = 1 - mean(data[[status_col]]),
    group_sizes = group_sizes,
    group_censoring = censoring_by_group,
    score_side_counts = c(above = sum(above_threshold), below = sum(!above_threshold)),
    threshold = 0,
    score_shift = 1
  ))
}
