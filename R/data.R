#' Real and Anonymus cleaned patient data contatining all information necessary to use the Single-Index-Model
#'
#'
#' @format ## `cleaned_data`
#' A data frame with 569 rows and 53 columns:
#' \describe{
#'   \item{time}{Failure times}
#'   \item{status}{A binary 1 for death and 0 for censoring}
#'   \item{trt}{A binary 1 for "recieved treatment" and 0 for "did not recieve treatment"}
#'   \item{bio1, bio2, ...}{Biomarkers}
#' }
#' @source Private
"cleaned_data"
