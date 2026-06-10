#' Kaplan-Meier threshold-group test for simulated data.
#'
#' @param data A threshold-model data frame with survival, treatment, and
#'   biomarker columns.
#' @param gamma A vector of biomarker coefficients.
#' @param time_col Column name for observed survival time.
#' @param status_col Column name for event status.
#' @param treatment_col Column name for treatment.
#' @param biomarker_cols Optional column names or indices for biomarkers.
#' @param strata Which groups to compare: threshold groups only, or the four
#'   treatment-by-threshold groups.
#'
#' @return A list with the Kaplan-Meier fit, log-rank test, p-value, and group
#'   counts.
#' @export
kaplan_meier_threshold_test <- function(data,
                                        gamma,
                                        time_col = "time",
                                        status_col = "status",
                                        treatment_col = "treatment",
                                        biomarker_cols = NULL,
                                        strata = c("threshold", "treatment_threshold")) {
  strata <- match.arg(strata)
  data <- as.data.frame(data)

  if (!all(c(time_col, status_col, treatment_col) %in% names(data))) {
    stop("time, status, and treatment columns must be present in data", call. = FALSE)
  }

  if (is.null(biomarker_cols)) {
    biomarker_cols <- setdiff(seq_along(data), match(c(time_col, status_col, treatment_col), names(data)))
  } else {
    biomarker_cols <- match_columns(data, biomarker_cols)
  }

  score <- biomarker_score(data[, biomarker_cols, drop = FALSE], gamma)
  threshold_group <- ifelse(score >= 0, "above", "below")

  if (identical(strata, "treatment_threshold")) {
    km_group <- paste0(ifelse(data[[treatment_col]] == 1, "treated", "control"), "_", threshold_group)
  } else {
    km_group <- threshold_group
  }

  km_data <- data.frame(
    time = data[[time_col]],
    status = data[[status_col]],
    group = factor(km_group),
    check.names = FALSE
  )

  surv_formula <- survival::Surv(time, status) ~ group
  km_fit <- survival::survfit(surv_formula, data = km_data)
  logrank <- survival::survdiff(surv_formula, data = km_data)
  degrees_freedom <- length(logrank$n) - 1L
  p_value <- stats::pchisq(logrank$chisq, df = degrees_freedom, lower.tail = FALSE)

  return(list(
    km_fit = km_fit,
    logrank = logrank,
    p_value = p_value,
    group_counts = table(km_data$group)
  ))
}

#' Goodness-of-fit diagnostics for a threshold Cox model.
#'
#' @param data A threshold-model data frame with survival, treatment, and
#'   biomarker columns.
#' @param gamma A vector of biomarker coefficients.
#' @param time_col Column name for observed survival time.
#' @param status_col Column name for event status.
#' @param treatment_col Column name or index for the treatment indicator after
#'   removing survival columns.
#' @param biomarker_cols Optional column names or indices for biomarkers.
#'
#' @return A one-row data frame containing model fit statistics, omnibus Cox
#'   test p-values, a proportional-hazards global p-value, and a Cox-Snell
#'   residual KS test against an exponential(1) distribution.
#' @export
threshold_fit_gof_tests <- function(data,
                                    gamma,
                                    time_col = "time",
                                    status_col = "status",
                                    treatment_col = "treatment",
                                    biomarker_cols = NULL) {
  data <- as.data.frame(data)

  if (!all(c(time_col, status_col) %in% names(data))) {
    stop("time and status columns must be present in data", call. = FALSE)
  }

  x <- data[, setdiff(names(data), c(time_col, status_col)), drop = FALSE]
  y <- survival::Surv(data[[time_col]], data[[status_col]])
  fit <- fit_threshold_cox(
    x = x,
    y = y,
    gamma = gamma,
    treatment_col = treatment_col,
    biomarker_cols = biomarker_cols
  )

  if (!isTRUE(fit$converged)) {
    return(data.frame(
      gof_converged = FALSE,
      loglik = NA_real_,
      aic = NA_real_,
      concordance = NA_real_,
      likelihood_ratio_p = NA_real_,
      wald_p = NA_real_,
      score_p = NA_real_,
      ph_global_p = NA_real_,
      cox_snell_ks_statistic = NA_real_,
      cox_snell_ks_p = NA_real_
    ))
  }

  fit_summary <- summary(fit)
  ph_test <- tryCatch(survival::cox.zph(fit), error = function(err) NULL)
  ph_global_p <- if (is.null(ph_test) || !"GLOBAL" %in% rownames(ph_test$table)) {
    NA_real_
  } else {
    unname(ph_test$table["GLOBAL", "p"])
  }

  cox_snell <- stats::fitted(fit, type = "expected")
  cox_snell_ks <- if (length(unique(cox_snell[is.finite(cox_snell)])) > 1L) {
    suppressWarnings(stats::ks.test(cox_snell, "pexp", rate = 1))
  } else {
    NULL
  }

  return(data.frame(
    gof_converged = TRUE,
    loglik = as.numeric(stats::logLik(fit)),
    aic = stats::AIC(fit),
    concordance = unname(fit_summary$concordance[1]),
    likelihood_ratio_p = unname(fit_summary$logtest["pvalue"]),
    wald_p = unname(fit_summary$waldtest["pvalue"]),
    score_p = unname(fit_summary$sctest["pvalue"]),
    ph_global_p = ph_global_p,
    cox_snell_ks_statistic = if (is.null(cox_snell_ks)) NA_real_ else unname(cox_snell_ks$statistic),
    cox_snell_ks_p = if (is.null(cox_snell_ks)) NA_real_ else cox_snell_ks$p.value
  ))
}
