#' Default MCMC parameters to control sampling and gamma's distribution.
#'
#' @param samples Number of posterior draws.
#' @param burn_in Number of burn-in (throw away) iterations.
#' @param thin Keep every `thin`-th draw (after `burn-in`).
#' @param gamma_mean Mean vector or scalar for the normal gamma prior.
#' @param gamma_sd Standard deviation vector or scalar for the normal gamma
#'   prior.
#' @param gamma_proposal_sd Random-walk proposal standard deviation for gamma.
#'
#' @return A named list of MCMC controls.
#' @export

default_mcmc_control <- function(samples = 5000,
                                 burn_in = 3000,
                                 thin = 10,
                                 gamma_mean = 0,
                                 gamma_sd = 0.1,
                                 gamma_proposal_sd = 0.05) {

  return(list(
    samples = samples,
    burn_in = burn_in,
    thin = thin,
    gamma_mean = gamma_mean,
    gamma_sd = gamma_sd,
    gamma_proposal_sd = gamma_proposal_sd
  ))
}

#' Metropolis-Hastings update.
#'
#' @param x A data frame of covariates with `treatment` column and one or more biomarker (|treatment|bio1|...|bion|)
#'   columns.
#' @param y A `survival::Surv` response object sorted in descending time order.
#' @param beta A 3-dimensional vector of the current regression coefficients.
#' @param gamma An n-dimensional vector of the current biomarker coefficients.
#' @param control A list from `default_mcmc_control()`.
#' @param lambda Nonnegative lasso penalty for gamma.
#' @param treatment_col Column name or index for the treatment indicator.
#' @param biomarker_cols Optional column names or indices for biomarkers.
#'
#' @return A list containing updated `beta`, `gamma`, and `log_likelihood`.
#' @export
mcmc_step <- function(x,
                      y,
                      beta,
                      gamma,
                      control,
                      lambda = 0,
                      treatment_col = 1,
                      biomarker_cols = NULL) {
  x <- as.data.frame(x) # incase x is not passed as a data frame
  previous <- list(beta = beta, gamma = gamma)

  # Loop through gamma selection
  current_log_likelihood <- cox_threshold_loglik(
    x,
    y,
    beta,
    gamma,
    treatment_col = treatment_col,
    biomarker_cols = biomarker_cols
  )
  gamma_proposal_sd <- repeat_parameter(control$gamma_proposal_sd, length(gamma))
  for (gamma_index in seq_along(gamma)) {
    # Get the current likelihood for this step
    current_gamma_logpost <- current_log_likelihood +
      log_gamma_prior(gamma, mean = control$gamma_mean, sd = control$gamma_sd) -
      lambda * sum(abs(gamma))

    # Copy gamma for the new candidate an replace the current index with one step
    candidate_gamma <- gamma
    candidate_gamma[gamma_index] <- stats::rnorm(
      1,
      mean = gamma[gamma_index],
      sd = gamma_proposal_sd[gamma_index]
    )
    candidate_log_likelihood <- cox_threshold_loglik(
      x,
      y,
      beta,
      candidate_gamma,
      treatment_col = treatment_col,
      biomarker_cols = biomarker_cols
    )
    candidate_gamma_logpost <- candidate_log_likelihood +
      log_gamma_prior(candidate_gamma, mean = control$gamma_mean, sd = control$gamma_sd) -
      lambda * sum(abs(candidate_gamma))

    if (accept_metropolis(candidate_gamma_logpost, current_gamma_logpost)) {
      gamma <- candidate_gamma
      current_log_likelihood <- candidate_log_likelihood
    }
  }


  # # gamma proposal step (helper function)
  # candidate_gamma <- propose_gamma(gamma, proposal_sd = control$gamma_proposal_sd)
  #
  # # Formula for logarithm of gamma posterior, candidate selection likelihood
  # candidate_gamma_logpost <- cox_threshold_loglik(
  #   x,
  #   y,
  #   beta,
  #   candidate_gamma,
  #   treatment_col = treatment_col,
  #   biomarker_cols = biomarker_cols
  # ) +
  #   log_gamma_prior(candidate_gamma, mean = control$gamma_mean, sd = control$gamma_sd) -
  #   lambda * sum(abs(candidate_gamma))

  # Criterion for selection (symmetric gamma distribution allows for this formula)
  # if (stats::runif(1) < exp(candidate_gamma_logpost - current_gamma_logpost)) {
  #   gamma <- candidate_gamma
  # }

  fit <- fit_threshold_cox(x, y, gamma, treatment_col, biomarker_cols)

  if (inherits(fit, "warning")) {
    beta <- previous$beta
    gamma <- previous$gamma
  } else {
    beta_loglik <- current_log_likelihood
    beta_candidate <- as.vector(mvtnorm::rmvnorm(1, beta, stats::vcov(fit)))
    if (all(is.finite(beta_candidate))) {
      candidate_loglik <- cox_threshold_loglik(
        x,
        y,
        beta_candidate,
        gamma,
        treatment_col = treatment_col,
        biomarker_cols = biomarker_cols
      )

      if (accept_metropolis(candidate_loglik, beta_loglik)) {
        beta <- beta_candidate
      }
    }
  }


  #log likelihood with new params
  log_likelihood <- cox_threshold_loglik(
    x,
    y,
    beta,
    gamma,
    treatment_col = treatment_col,
    biomarker_cols = biomarker_cols
  )

  return(list(beta = beta, gamma = gamma, log_likelihood = log_likelihood))
}

#' Fit the model using the new mcmc method.
#'
#' @param x A data frame of covariates with treatment column and one or more biomarker (|treatment|bio1|...|bion|)
#'   columns.
#' @param y A `survival::Surv` response object sorted in descending time order.
#' @param control A list from `default_mcmc_control()`.
#' @param lambda Nonnegative lasso penalty for gamma.
#' @param gamma_start Optional initial gamma vector.
#' @param treatment_col Column name or index for the treatment indicator.
#' @param biomarker_cols Optional column names or indices for biomarkers.
#'
#' @return Posterior samples for beta, gamma, and log likelihood.
#' @export
fit_threshold_mcmc <- function(x,
                               y,
                               control = default_mcmc_control(),
                               lambda = 0,
                               gamma_start = NULL,
                               treatment_col = 1,
                               biomarker_cols = NULL) {
  if (lambda < 0) {
    stop("lambda must be nonnegative", call. = FALSE)
  }

  x <- as.data.frame(x)
  biomarker_cols <- resolve_biomarker_columns(x, treatment_col, biomarker_cols)
  n_biomarkers <- length(biomarker_cols)

  gamma <- if (is.null(gamma_start)) {
    draw_gamma(n_biomarkers, mean = control$gamma_mean, sd = control$gamma_sd)
  } else {
    validate_gamma(gamma_start, n_biomarkers)
  }
  beta <- initialize_beta(x, y, gamma, treatment_col, biomarker_cols)
  state <- list(beta = beta, gamma = gamma)

  for (i in seq_len(control$burn_in)) {
    state <- mcmc_step(
      x,
      y,
      beta = state$beta,
      gamma = state$gamma,
      control = control,
      lambda = lambda,
      treatment_col = treatment_col,
      biomarker_cols = biomarker_cols
    )
  }

  retained <- control$samples
  total_iterations <- retained * control$thin
  beta_samples <- matrix(NA_real_, nrow = length(state$beta), ncol = retained)
  gamma_samples <- matrix(NA_real_, nrow = n_biomarkers, ncol = retained)
  loglik_samples <- matrix(NA_real_, nrow = 1, ncol = retained)

  sample_index <- 1
  for (iteration in seq_len(total_iterations)) {
    state <- mcmc_step(
      x,
      y,
      beta = state$beta,
      gamma = state$gamma,
      control = control,
      lambda = lambda,
      treatment_col = treatment_col,
      biomarker_cols = biomarker_cols
    )

    if (iteration %% control$thin == 0) {
      beta_samples[, sample_index] <- state$beta
      gamma_samples[, sample_index] <- state$gamma
      loglik_samples[, sample_index] <- state$log_likelihood
      sample_index <- sample_index + 1
    }
  }

  list(
    beta_samples = beta_samples,
    gamma_samples = gamma_samples,
    loglik_samples = loglik_samples,
    threshold = 0,
    score_shift = 1,
    control = control
  )
}

#' Fit the fixed-threshold model from a survival data frame.
#'
#' @param data Data frame containing survival `time`, `status`, `treatment` (they MUST have these names), and
#'   biomarker columns.
#' @param control A list from `default_mcmc_control()`.
#' @param lambda Nonnegative lasso penalty for gamma.
#' @param gamma_start Optional initial gamma vector.
#' @param time_col Column name for survival time.
#' @param status_col Column name for event status.
#' @param treatment_col Column name for treatment.
#' @param biomarker_cols Optional column names or indices for biomarkers.
#'
#' @return Posterior samples from `fit_threshold_mcmc()`.
#' @export
fit_threshold_model <- function(data,
                                control = default_mcmc_control(),
                                lambda = 0,
                                gamma_start = NULL,
                                time_col = "time",
                                status_col = "status",
                                treatment_col = "treatment",
                                biomarker_cols = NULL) {
  data <- as.data.frame(data)

  # check for basic columns
  required_cols <- c(time_col, status_col, treatment_col)
  if (!all(required_cols %in% names(data))) {
    stop("data must contain time, status, and treatment columns", call. = FALSE)
  }

  # extract biomarker columns
  if (is.null(biomarker_cols)) {
    biomarker_cols <- setdiff(seq_along(data), match(required_cols, names(data)))
  } else {
    biomarker_cols <- match_columns(data, biomarker_cols)
  }

  # make the survival object for mcmc
  y <- survival::Surv(data[[time_col]], data[[status_col]])

  # set up the data for mcmc
  x <- data[, c(treatment_col, names(data)[biomarker_cols]), drop = FALSE]

  # sort y and feed into mcmc for return.
  sorted <- sort(y[, 1], decreasing = TRUE, index.return = TRUE)
  return(fit_threshold_mcmc(
    x = x[sorted$ix, , drop = FALSE],
    y = y[sorted$ix, ],
    control = control,
    lambda = lambda,
    gamma_start = gamma_start,
    treatment_col = 1,
    biomarker_cols = seq_len(length(biomarker_cols)) + 1
  ))
}

#' Summarize posterior samples.
#'
#' @param posterior A posterior sample list from `fit_threshold_mcmc()` (or `from fit_threshold_model()`).
#' @param alpha Error probability for equal-tailed credible intervals.
#' @param method Point estimate method, either `"mean"` or `"median"`.
#'
#' @return A list of posterior point estimates and credible intervals.
#' @export
summarize_mcmc <- function(posterior, alpha = 0.05, method = "mean") {
  if (!method %in% c("mean", "median")) {
    stop("method must be either 'mean' or 'median'", call. = FALSE)
  }

  estimator <- match.fun(method)
  tail_prob <- alpha / 2
  interval_probs <- c(tail_prob, 1 - tail_prob)

  return(list(
    beta = apply(posterior$beta_samples, 1, estimator),
    gamma = apply(posterior$gamma_samples, 1, estimator),
    beta_interval = apply(posterior$beta_samples, 1, stats::quantile, interval_probs),
    gamma_interval = apply(posterior$gamma_samples, 1, stats::quantile, interval_probs),
    threshold = posterior$threshold,
    score_shift = posterior$score_shift
  ))
}

#' Decide whether to accept a Metropolis-Hastings proposal.
#'
#' Non-finite log ratios are rejected, which prevents `NA` acceptance
#' probabilities from crashing the sampler.
#'
#' @param candidate_logpost Candidate log posterior or log likelihood.
#' @param current_logpost Current log posterior or log likelihood.
#'
#' @return A logical scalar.
#' @keywords internal
accept_metropolis <- function(candidate_logpost, current_logpost) {
  log_ratio <- candidate_logpost - current_logpost
  if (!is.finite(log_ratio)) {
    return(FALSE)
  }

  log(stats::runif(1)) < min(0, log_ratio)
}
