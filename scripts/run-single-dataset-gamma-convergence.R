#!/usr/bin/env Rscript

parse_args <- function(args) {
  values <- list(
    data_n = "600",
    output_dir = NA_character_,
    seed = "20260728",
    max_data_set_attempts = "100",
    min_group_size = "60",
    samples = "20000",
    burn_in = "10000",
    thin = "10",
    gamma_proposal_sd = "0.5"
  )

  option_names <- c(
    "--data-n",
    "--output-dir",
    "--seed",
    "--max-data-set-attempts",
    "--min-group-size",
    "--samples",
    "--burn-in",
    "--thin",
    "--gamma-proposal-sd"
  )

  index <- 1L
  while (index <= length(args)) {
    arg <- args[[index]]
    if (arg %in% option_names) {
      if (index == length(args)) {
        stop("missing value for argument: ", arg, call. = FALSE)
      }
      key <- gsub("-", "_", sub("^--", "", arg))
      values[[key]] <- args[[index + 1L]]
      index <- index + 2L
    } else {
      stop("unknown argument: ", arg, call. = FALSE)
    }
  }

  values$data_n <- as.integer(values$data_n)
  values$seed <- as.integer(values$seed)
  values$max_data_set_attempts <- as.integer(values$max_data_set_attempts)
  values$min_group_size <- as.integer(values$min_group_size)
  values$samples <- as.integer(values$samples)
  values$burn_in <- as.integer(values$burn_in)
  values$thin <- as.integer(values$thin)
  values$gamma_proposal_sd <- as.numeric(values$gamma_proposal_sd)

  if (is.na(values$output_dir)) {
    values$output_dir <- sprintf("results/single-dataset-gamma-convergence-n%d", values$data_n)
  }

  values
}

source_package_r_files <- function() {
  r_files <- list.files("R", pattern = "[.][rR]$", full.names = TRUE)
  for (r_file in r_files) {
    source(r_file)
  }
}

simulate_eligible_data_set <- function(args, true_beta, true_gamma) {
  attempts <- vector("list", args$max_data_set_attempts)

  for (attempt in seq_len(args$max_data_set_attempts)) {
    data <- simulate_threshold_data(
      n = args$data_n,
      beta = true_beta,
      gamma = true_gamma,
      baseline_hazard = 1,
      biomarker_correlation = 0,
      study_end_range = c(0.4, 3)
    )
    data_summary <- summarize_simulated_data(data, true_gamma)

    attempts[[attempt]] <- data.frame(
      attempt = attempt,
      censoring_proportion = data_summary$overall_censoring,
      treated_above = data_summary$group_sizes[["treated_above"]],
      treated_below = data_summary$group_sizes[["treated_below"]],
      control_above = data_summary$group_sizes[["control_above"]],
      control_below = data_summary$group_sizes[["control_below"]],
      eligible = all(data_summary$group_sizes >= args$min_group_size)
    )

    if (attempts[[attempt]]$eligible) {
      return(list(
        data = data,
        data_summary = data_summary,
        attempt = attempt,
        attempts = do.call(rbind, attempts[seq_len(attempt)])
      ))
    }
  }

  stop(
    "no eligible data set found after ",
    args$max_data_set_attempts,
    " attempts; try reducing --min-group-size",
    call. = FALSE
  )
}

gamma_trace_data <- function(posterior, true_gamma) {
  gamma_samples <- posterior$gamma_samples
  n_samples <- ncol(gamma_samples)
  n_gamma <- nrow(gamma_samples)
  iteration <- seq_len(n_samples) * posterior$control$thin + posterior$control$burn_in

  do.call(rbind, lapply(seq_len(n_gamma), function(gamma_index) {
    samples <- gamma_samples[gamma_index, ]
    running_mean <- cumsum(samples) / seq_along(samples)

    data.frame(
      retained_sample = seq_len(n_samples),
      iteration = iteration,
      parameter = paste0("gamma_", gamma_index),
      value = samples,
      running_mean = running_mean,
      truth = true_gamma[gamma_index],
      check.names = FALSE
    )
  }))
}

plot_gamma_convergence <- function(trace, output_file, width = 1400, height = 900, res = 140) {
  grDevices::png(output_file, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_gamma_convergence(trace)
}

plot_gamma_convergence_pdf <- function(trace, output_file, width = 10, height = 7) {
  grDevices::pdf(output_file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_gamma_convergence(trace)
}

draw_gamma_convergence <- function(trace) {
  parameters <- unique(trace$parameter)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::par(
    mfrow = c(length(parameters), 1L),
    mar = c(4, 4.2, 2.5, 1),
    oma = c(0, 0, 2, 0)
  )

  for (parameter in parameters) {
    parameter_trace <- trace[trace$parameter == parameter, ]
    y_range <- range(
      parameter_trace$value,
      parameter_trace$running_mean,
      parameter_trace$truth,
      finite = TRUE
    )

    graphics::plot(
      parameter_trace$iteration,
      parameter_trace$value,
      type = "l",
      col = grDevices::adjustcolor("gray35", alpha.f = 0.45),
      lwd = 0.7,
      xlab = "MCMC iteration",
      ylab = parameter,
      ylim = y_range,
      main = paste(parameter, "trace and running mean")
    )
    graphics::lines(
      parameter_trace$iteration,
      parameter_trace$running_mean,
      col = "#0072B2",
      lwd = 2
    )
    graphics::abline(h = parameter_trace$truth[1], col = "#D55E00", lwd = 2, lty = 2)
    graphics::legend(
      "topright",
      legend = c("trace", "running mean", "truth"),
      col = c(grDevices::adjustcolor("gray35", alpha.f = 0.45), "#0072B2", "#D55E00"),
      lty = c(1, 1, 2),
      lwd = c(1, 2, 2),
      bty = "n"
    )
  }

  graphics::mtext("Gamma convergence for one simulated data set", outer = TRUE, font = 2)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
set.seed(args$seed)
source_package_r_files()

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

true_beta <- c(log(1.5), log(1.2), log(2))
true_gamma <- c(1.5, 1.8)
gamma_start <- c(1, 1.3)

control <- default_mcmc_control(
  samples = args$samples,
  burn_in = args$burn_in,
  thin = args$thin,
  gamma_mean = c(0.5, 0.5),
  gamma_sd = 1,
  gamma_proposal_sd = args$gamma_proposal_sd
)

message("simulating one eligible data set")
simulated <- simulate_eligible_data_set(args, true_beta, true_gamma)

message("fitting MCMC")
posterior <- suppressWarnings(fit_threshold_model(
  data = simulated$data,
  control = control,
  lambda = 0,
  gamma_start = gamma_start
))
summary <- summarize_mcmc(posterior)
trace <- gamma_trace_data(posterior, true_gamma)

data_file <- file.path(args$output_dir, "single-dataset.csv")
attempts_file <- file.path(args$output_dir, "data-set-attempts.csv")
summary_file <- file.path(args$output_dir, "single-dataset-summary.csv")
trace_file <- file.path(args$output_dir, "gamma-trace.csv")
plot_file <- file.path(args$output_dir, "gamma-convergence.png")
plot_pdf_file <- file.path(args$output_dir, "gamma-convergence.pdf")
posterior_file <- file.path(args$output_dir, "single-dataset-posterior.rds")

utils::write.csv(simulated$data, data_file, row.names = FALSE)
utils::write.csv(simulated$attempts, attempts_file, row.names = FALSE)
utils::write.csv(trace, trace_file, row.names = FALSE)
saveRDS(posterior, posterior_file)

summary_row <- data.frame(
  data_n = args$data_n,
  seed = args$seed,
  data_set_attempt = simulated$attempt,
  samples = args$samples,
  burn_in = args$burn_in,
  thin = args$thin,
  gamma_proposal_sd = args$gamma_proposal_sd,
  true_gamma_1 = true_gamma[1],
  true_gamma_2 = true_gamma[2],
  gamma_start_1 = gamma_start[1],
  gamma_start_2 = gamma_start[2],
  gamma_mean_1 = summary$gamma[1],
  gamma_mean_2 = summary$gamma[2],
  gamma_1_interval_lower = summary$gamma_interval[1, 1],
  gamma_1_interval_upper = summary$gamma_interval[2, 1],
  gamma_2_interval_lower = summary$gamma_interval[1, 2],
  gamma_2_interval_upper = summary$gamma_interval[2, 2],
  censoring_proportion = simulated$data_summary$overall_censoring,
  treated_above = simulated$data_summary$group_sizes[["treated_above"]],
  treated_below = simulated$data_summary$group_sizes[["treated_below"]],
  control_above = simulated$data_summary$group_sizes[["control_above"]],
  control_below = simulated$data_summary$group_sizes[["control_below"]]
)
utils::write.csv(summary_row, summary_file, row.names = FALSE)

plot_gamma_convergence(trace, plot_file)
plot_gamma_convergence_pdf(trace, plot_pdf_file)

message("wrote data: ", data_file)
message("wrote attempts: ", attempts_file)
message("wrote summary: ", summary_file)
message("wrote trace: ", trace_file)
message("wrote posterior: ", posterior_file)
message("wrote plot: ", plot_file)
message("wrote PDF plot: ", plot_pdf_file)
