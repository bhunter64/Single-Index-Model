test_that("test simulation works for data sets", {
  dat <-simu_data(400, beta = c(log(3), log(2), log(0.5)), gamma = c(-1.8, -2, 0.9))

  par(mfrow = c (1, 2))

  #Define the sensitive and unsensitive groups
  grp1 <- dat[which(dat$sensitivity == 1), ]
  grp2 <- dat[which(dat$sensitivity == 0), ]

  x_1 <- grp1$time
  delta_1 <- grp1$status
  fit_1 <- survival::survfit(formula = Surv(x_1, delta_1) ~ 1)
  plot(fit_1, xlab = "Time", ylab = "Estimated Survival of Sensitive group", mark.time=T, conf.int=F)

  x_2 <- grp2$time
  delta_2 <- grp2$status
  fit_2 <- survival::survfit(formula = Surv(x_2, delta_2) ~ 1)
  plot(fit_2, xlab = "Time", ylab = "Estimated Survival of Unsensitive group", mark.time=T, conf.int=F)

  #COX
  fit_cox<-survival::coxph(survival::Surv(time,status) ~ as.factor(trt)
                 + as.factor(sensitivity) + as.factor(both), data=dat)
  print(summary(fit_cox))

})
