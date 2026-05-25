#' Simulate study data based on predetermined coefficients (beta, gamma)
#'
#' @param n The number of observations (subjects).
#' @param beta A vector of effects: 1-Trt & not sensitive. 2-Sensitive & no trt. 1+2+3-Sensitive & trt
#' @param gamma A vector of biomarker effects on individuals sensitvity
#' @returns A data frame of time | status | trt | sensitivity | biomarker 1 | ... | biomarker n |
simu_data <- function (n, beta, gamma, h0 = 1) {
  # Biomarker number implicitly determined by gamma length
  num_of_biomarkers <- length(gamma)

  # Treatment simulation
  x <- stats::rbinom(n, 1, 0.5)

  # MVN Biomarker Simulation
  mean <- rep(0, num_of_biomarkers)
  rho <- 0.5
  sigma <- matrix(rho, ncol=num_of_biomarkers, nrow=num_of_biomarkers)
  diag(sigma) <- 1
  bs <- mvtnorm::rmvnorm(n, mean=mean, sigma=sigma)
  colnames(bs) <- paste0('bio',c(1:num_of_biomarkers))

  # Sensitivity Simulation
  tmpX <- data.frame(trt=x, bs)
  X <- thresholding(tmpX, gamma) #tool to test sensitivity on the given gamma
  sensitivity <- X[,2]
  both <- X[,3]

  # Hazard Calculation and Survival Simulation
  h <- h0 * exp(as.matrix(X) %*% beta)
  stime <- stats::rexp(n, h) # survival time ~ Exp(1/2 * exp(x beta))

  # ——— OPTIONAL WEIBULL ALTERNATIVE———
  # survival time ~ Weibull(4, 1/2 * exp(x beta))
  #a <- 4
  #hb <- exp(-log(h) / a)
  #stime <- rweibull(n, a, hb) # weibull??

  # Simulates end of study period
  endstudy <- stats::runif(n, 0.05, 3)

  # Sets Censoring Indicator
  cstatus <- ifelse(stime > endstudy, 0, 1)
  stime <- ifelse(stime > endstudy, endstudy, stime)

  # Build the final set and return
  dat <- data.frame(cbind(time = stime, status = cstatus, trt = x, sensitivity = sensitivity, both = both, bs))
  return(dat)
}
