#' Evaluates the Cox partial log likelihood.
#'
#' @param x A data frame of covariates with `treatment` column and one or more biomarker (|treatment|bio1|...|bion|)
#'   columns.
#' @param y A `survival::Surv` response object sorted in descending time order.
#' @param beta A vector of regression coefficients for treatment,
#'   biomarker subgroup, and their interaction term.
#' @param gamma A vector of biomarker coefficients.
#' @param treatment_col Column name or index for the treatment indicator.
#' @param biomarker_cols Optional column names or indices for biomarkers.
#'
#' @return A scalar log partial likelihood.
#' @export
cox_threshold_loglik <- function(x,
                                 y,
                                 beta,
                                 gamma,
                                 treatment_col = 1,
                                 biomarker_cols = NULL) {
  design <- make_threshold_design(
    x = x,
    gamma = gamma,
    treatment_col = treatment_col,
    biomarker_cols = biomarker_cols
  )

  linear_predictor <- as.matrix(design) %*% beta
  risk_cusum <- log(cumsum(exp(linear_predictor)))
  return(as.numeric(sum(y[, 2] * (linear_predictor - risk_cusum))))
}

#' Fit the betas in the Cox model for a given gamma.
#'
#' @param x  A data frame of covariates with `treatment` column and one or more biomarker (|treatment|bio1|...|bion|)
#'   columns.
#' @param y A `survival::Surv` response object.
#' @param gamma Numeric vector of biomarker coefficients.
#' @param treatment_col Column name or index for the treatment indicator.
#' @param biomarker_cols Optional column names or indices for biomarkers.
#'
#' @return A `coxph` fit with an added logical `converged` field. If `coxph`
#'   warns, the warning object is returned with `converged = FALSE`.
#' @export
fit_threshold_cox <- function(x,
                              y,
                              gamma,
                              treatment_col = 1,
                              biomarker_cols = NULL) {
  design <- make_threshold_design(
    x = x,
    gamma = gamma,
    treatment_col = treatment_col,
    biomarker_cols = biomarker_cols
  )

  # standard fit object for cox regression with an additional `converged` argument
  fit <- tryCatch(
    survival::coxph(
      y ~ treatment + subgroup + treatment_subgroup,
      data = design,
      control = survival::coxph.control(iter.max = 1000), singular.ok = FALSE
    ),
    warning = function(w) w
  )

  if (inherits(fit, "warning")) {
    warning(
      "cox regression failed to converge on beta sampling: ",
      conditionMessage(fit),
      call. = FALSE
    )
    fit$converged <- FALSE
  } else {
    fit$converged <- TRUE
  }

  return(fit)
}

#' Initialize beta with the Cox fit.
#'
#' @param x A data frame of covariates with `treatment` column and one or more biomarker (|treatment|bio1|...|bion|)
#'   columns.
#' @param y A `survival::Surv` response object.
#' @param gamma A vector of biomarker coefficients.
#' @param treatment_col Column name or index for the treatment indicator.
#' @param biomarker_cols Optional column names or indices for biomarkers.
#'
#' @return beta coefficients from the Cox model.
#' @export
initialize_beta <- function(x,
                            y,
                            gamma,
                            treatment_col = 1,
                            biomarker_cols = NULL) {
  fit <- fit_threshold_cox(x, y, gamma, treatment_col, biomarker_cols)
  if (!isTRUE(fit$converged)) {
    stop("could not initialize beta because the Cox model did not converge", call. = FALSE)
  }

  return(stats::coef(fit))
}
