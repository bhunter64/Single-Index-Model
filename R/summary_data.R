#' Summarize data
#'
#' @param data A data frame: time | status | trt | sensitivity | biomarker 1 | ... | biomarker n |.
#' @returns A list:
summary_data <- function (data) {
  time <- data$time
  trt <- data$trt
  idx <- grep("b.", colnames(data))
  bs <- data[, idx]
  sensitivity <- data$sensitivity

  grp1 <- data[which(trt==1 & sensitivity == 1), ] # trt and 1+gamma1*z1+gamma2*z2>=0
  grp2 <- data[which(trt==1 & sensitivity == 0), ]  # trt and 1+gamma1*z1+gamma2*z2<0
  grp3 <- data[which(trt==0 & sensitivity == 1), ] # control and 1+gamma1*z1+gamma2*z2>=0
  grp4 <- data[which(trt==0 & sensitivity == 0), ]  # control and 1+gamma1*z1+gamma2*Z2<0

  aboveCount <- nrow(data[which(sensitivity == 1), ])
  belowCount <- nrow(data[which(sensitivity == 0), ])

  above_below_count <- c(aboveCount, belowCount)

  status_count <- c(nrow(grp1), nrow(grp2), nrow(grp3), nrow(grp4))

  censor_prop1 <- 1 - sum(grp1$status) / nrow(grp1)
  censor_prop2 <- 1 - sum(grp2$status) / nrow(grp2)
  censor_prop3 <- 1 - sum(grp3$status) / nrow(grp3)
  censor_prop4 <- 1 - sum(grp4$status) / nrow(grp4)
  cens_status_prop <- c(grp1=censor_prop1, grp2=censor_prop2, grp3=censor_prop3, grp4=censor_prop4) # 20%-35%

  names(status_count) <- c("group1", "group2", "group3", "group4")
  names(cens_status_prop) <- c("group1", "group2", "group3", "group4")

  cens_prop <- 1 - sum(data$status) / nrow(data)

  return(list(cens_prop=cens_prop,status_count=status_count,
              cens_status_prop=cens_status_prop, above_below_count=above_below_count))
}
