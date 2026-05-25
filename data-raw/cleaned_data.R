## code to prepare `cleaned_data` dataset goes here
cleaned_data <- PA3_bio_base_surv_data
# extract patients with plasma samples
plasma_na <- apply(cleaned_data[,11:ncol(cleaned_data)],1,FUN=function(x){sum(is.na(x))})
cleaned_data <- cleaned_data[plasma_na<50,]
colnames(cleaned_data)[5] <- "status"
cleaned_data$trt <- as.factor(ifelse(cleaned_data$allo1_cd=='D',1,0))
cleaned_data <- cleaned_data %>% relocate(trt, .after = status)
cleaned_data <- cleaned_data[,c(4,5,6, 12:ncol(cleaned_data))]

bmIdx <- 4:(ncol(cleaned_data)-1)

cleaned_data <- cleaned_data %>%
  mutate(across(all_of(bmIdx), ~ {
    num_col <- as.numeric(.x)                  # 1. Force to numbers (turns "." into NA)
    replace_na(num_col, mean(num_col, na.rm = TRUE)) # 2. Replace NAs with the mean
  }))


cleaned_data[,bmIdx] <- log(cleaned_data[,bmIdx])
cleaned_data[,bmIdx] <- scale(cleaned_data[,bmIdx])


usethis::use_data(cleaned_data, overwrite = TRUE)
