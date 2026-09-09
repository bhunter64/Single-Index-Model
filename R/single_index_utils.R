#' Create the fixed-threshold model design matrix.
#'
#' Builds the frame with treatment, biomarker-threshold-subgroup, and treatment-and-subgroup columns
#' used by the Cox threshold model. Biomarker scores are computed directly from
#' `1 + biomarkers %*% gamma` at greater than/equal to and strictly less than zero.
#'
#' @param x A data frame (or matrix) with a treatment column and biomarker columns.
#' @param gamma A vector of biomarker coefficients.
#' @param treatment_col Column name or index for the treatment indicator.
#' @param biomarker_cols Optional column names or indices for biomarkers. When
#'   omitted, EVERY COLUMN EXCEPT `treatment_col` IS TREATED AS A BIOMARKER.
#'
#' @return A data frame with columns `treament`, `subgroup`, and `treatment_subgroup`.
#' @export
make_threshold_design <- function(x,
                                  gamma,
                                  treatment_col = 1,
                                  biomarker_cols = NULL) {
  x <- as.data.frame(x)
  # helper to find biomarker columns
  biomarker_cols <- resolve_biomarker_columns(x, treatment_col, biomarker_cols)

  # helper to calculate 1 + biomarkers %*% gamma so
  score <- biomarker_score(x[, biomarker_cols, drop = FALSE], gamma)

  # calculate over/under the 0-threshold
  subgroup <- as.integer(score >= 0)
  treatment <- x[[treatment_col]]

  return(data.frame(
    treatment = treatment,
    subgroup = subgroup,
    treatment_subgroup = treatment * subgroup
  ))
}

#' Compute the gamma biomarker score for the threshold. Also uses gamma validation.
#'
#' @param biomarkers A data frame (or matrix) of biomarker values.
#' @param gamma A vector of biomarker coefficients.
#'
#' @return A vector containing `1 + biomarkers %*% gamma`.
#' @export
biomarker_score <- function(biomarkers, gamma) {
  biomarkers <- as.matrix(biomarkers)
  # validate gamma length before computation
  gamma <- validate_gamma(gamma, ncol(biomarkers))
  return(as.vector(1 + biomarkers %*% gamma))
}

#' Validate gamma dimensions for calculations and then return it if it passes.
#'
#' @param gamma A vector of biomarker coefficients.
#' @param n_biomarkers Integer number of biomarkers.
#'
#' @return `gamma` as a numeric vector.
#' @export
validate_gamma <- function(gamma, n_biomarkers) {
  gamma <- as.numeric(gamma)
  if (length(gamma) != n_biomarkers) {
    stop("gamma must have one coefficient per biomarker", call. = FALSE)
  }
  if (anyNA(gamma) || any(!is.finite(gamma))) {
    stop("gamma must contain only finite and numeric values", call. = FALSE)
  }
  return(gamma)
}

#' Resolve biomarker columns for model design.
#'
#' @param x A data frame (or matrix) of model covariates.
#' @param treatment_col Column name or index for the treatment indicator.
#' @param biomarker_cols Optional column names or indices for biomarkers.
#'
#' @return Column indices for the biomarker columns.
#' @export
resolve_biomarker_columns <- function(x, treatment_col = 1, biomarker_cols = NULL) {
  x <- as.data.frame(x)
  if (is.null(biomarker_cols)) {
    biomarker_cols <- setdiff(seq_along(x), match_columns(x, treatment_col))
  }

  biomarker_cols <- match_columns(x, biomarker_cols)

  # check if biomarker columns exist
  if (!length(biomarker_cols)) {
    stop("at least one biomarker column is required", call. = FALSE)
  }
  return(biomarker_cols)
}

#' Draw gamma from a normal prior distribution.
#'
#' @param n_biomarkers Number of biomarker coefficients (=gamma components) to draw.
#' @param mean Mean vector (or scalar if identical) for the normal prior.
#' @param sd Standard deviation vector (or scalar if identical) for the normal prior.
#'
#' @return A numeric gamma vector.
#' @export
draw_gamma <- function(n_biomarkers, mean = 0, sd = 0.1) {
  stats::rnorm(n_biomarkers, mean = repeat_parameter(mean, n_biomarkers), sd = repeat_parameter(sd, n_biomarkers))
}

#' Draw gamma from a uniform prior distribution.
#'
#' @param n_biomarkers Number of biomarker coefficients (=gamma components) to draw.
#' @param min Lower bound vector (or scalar if identical) for the uniform prior.
#' @param max Upper bound vector (or scalar if identical) for the uniform prior.
#'
#' @return A numeric gamma vector.
#' @export
draw_gamma_uniform <- function(n_biomarkers, min = -1, max = 1) {
  bounds <- validate_uniform_bounds(min, max, n_biomarkers)
  stats::runif(n_biomarkers, min = bounds$min, max = bounds$max)
}

#' Evaluate the normal gamma prior density.
#'
#' @param gamma A vector of biomarker coefficients.
#' @param mean Mean vector (or scalar if identical) for the normal prior.
#' @param sd Standard deviation vector (or scalar if identical) for the normal prior.
#'
#' @return The log prior density.
#' @export
log_gamma_prior <- function(gamma, mean = 0, sd = 0.1) {
  sum(stats::dnorm(gamma, mean = repeat_parameter(mean, length(gamma)), sd = repeat_parameter(sd, length(gamma)), log = TRUE))
}

#' Evaluate the uniform gamma prior density.
#'
#' Each gamma component has an independent uniform prior. Values outside any
#' component's bounds have log density `-Inf` and are therefore rejected by the
#' Metropolis-Hastings update.
#'
#' @param gamma A vector of biomarker coefficients.
#' @param min Lower bound vector (or scalar if identical) for the uniform prior.
#' @param max Upper bound vector (or scalar if identical) for the uniform prior.
#'
#' @return The log prior density.
#' @export
log_gamma_uniform_prior <- function(gamma, min = -1, max = 1) {
  bounds <- validate_uniform_bounds(min, max, length(gamma))
  sum(stats::dunif(gamma, min = bounds$min, max = bounds$max, log = TRUE))
}

#' Propose an updated gamma vector.
#'
#' @param gamma Current gamma vector.
#' @param proposal_sd Random-walk proposal standard deviation.
#'
#' @return A proposed gamma vector.
#' @export
propose_gamma <- function(gamma, proposal_sd = 0.05) {
  stats::rnorm(length(gamma), mean = gamma, sd = repeat_parameter(proposal_sd, length(gamma)))
}

#' Match requested columns to numeric indices.
#'
#' @param x A data frame or matrix.
#' @param cols Column names or indices.
#'
#' @return Integer column indices.
#' @keywords internal
match_columns <- function(x, cols) {
  if (is.character(cols)) {
    idx <- match(cols, names(x))
    # check for existence of the name
    if (anyNA(idx)) {
      stop("unknown column name in column matching for the model. Please review your data", call. = FALSE)
    }
    return(idx)
  }

  cols <- as.integer(cols)

  # check if indicies are valid
  if (any(cols < 1 | cols > ncol(x))) {
    stop("column index is outside the data", call. = FALSE)
  }
  return(cols)
}

#' Repeat a scalar distribution parameter.
#'
#' @param value Scalar or vector distribution parameter.
#' @param n Required parameter length.
#'
#' @return A vector of length `n`.
#' @keywords internal
repeat_parameter <- function(value, n) {
  if (length(value) == 1L) {
    return(rep(value, n))
  }
  if (length(value) != n) {
    stop("distribution parameter length must be 1 or match gamma length", call. = FALSE)
  }
  value
}

#' Validate component-wise uniform distribution bounds.
#'
#' @param min Lower bound vector or scalar.
#' @param max Upper bound vector or scalar.
#' @param n Required bounds length.
#'
#' @return A list with expanded `min` and `max` vectors.
#' @keywords internal
validate_uniform_bounds <- function(min, max, n) {
  min <- repeat_parameter(min, n)
  max <- repeat_parameter(max, n)

  if (anyNA(min) || anyNA(max) || any(!is.finite(min)) || any(!is.finite(max))) {
    stop("uniform gamma prior bounds must be finite", call. = FALSE)
  }
  if (any(min >= max)) {
    stop("each uniform gamma prior minimum must be less than its maximum", call. = FALSE)
  }

  list(min = min, max = max)
}
