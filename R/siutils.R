#' Combine biomarkers with gamma linearly based on the single-index model
#'
#' @param biomarkers A data frame of row(i)=subjects_i, col(j)=biomarker_j.
#' @param gamma A vector of biomarker effects on individuals sensitvity
#' @returns Combined threshold value
#' @examples
cmb_biomarker <- function(biomarkers, gamma) {
  val <- (as.matrix(biomarkers, nrow = nrow(biomarkers), ncol = ncol(biomarkers)) %*% gamma) + 1
  return(val)
}

#' Basic helper to extract Biomarker index
#'
#' @param x A data frame of study data (see format elsewhere)
#' @returns Index where biomarkers are foubd
#' @examples
getBM_col <- function(x) {
  idx <- grep("bio", colnames(x))
  return(idx)
}

#' Combine biomarkers with gamma linearly based on the single-index model
#'
#' @param x A data frame with trt and biomarker values.
#' @param gamma A vector of biomarker effects on individuals sensitvity
#' @returns A data frame: trt | sensitivity | product of trt*sensitivity
#' @examples
thresholding <- function(x, gamma) {
  bm_idx <- getBM_col(x)
  w <- cmb_biomarker(x[,bm_idx], gamma)
  indicator_sensitivity <- ifelse(0 <= w, 1, 0)
  trt_and_sensitive <- x[,1] * indicator_sensitivity
  x_thresholded <- data.frame(trt=x[,1],
                              indicator_sensitivity=indicator_sensitivity,
                              trt_and_sensitive=trt_and_sensitive)
  return(x_thresholded)
}

