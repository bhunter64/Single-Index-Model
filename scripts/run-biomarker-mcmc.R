#!/usr/bin/env Rscript

parse_numeric_vector <- function(value, option, expected_length = NULL) {
  parsed <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (anyNA(parsed) || any(!is.finite(parsed))) {
    stop(option, " must be a comma-separated list of finite numbers", call. = FALSE)
  }
  if (!is.null(expected_length) && length(parsed) != expected_length) {
    stop(option, " must contain exactly ", expected_length, " values", call. = FALSE)
  }
  parsed
}

parse_args <- function(args) {
  absolute_data <- "/data/cleaned_data.rda"
  default_data <- if (file.exists(absolute_data)) absolute_data else "data/cleaned_data.rda"
  values <- list(
    data = default_data,
    data_object = "cleaned_data",
    output_dir = "results/biomarker-mcmc",
    time_col = "time",
    status_col = "status",
    treatment_col = "trt",
    biomarker_cols = "bio30,bio3,bio37",
    biomarker_labels = "HER-2,CA19-9,Axl",
    gamma_start = "-0.6,0.5,0.5",
    beta_start = "auto",
    reference_gamma = "none",
    reference_beta = "none",
    chains = "4",
    workers = Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"),
    samples = "5000",
    burn_in = "5000",
    thin = "5",
    seed = "20260811",
    gamma_prior_mean = "0,0,0",
    gamma_prior_sd = "1,1,1",
    gamma_proposal_sd = "0.05,0.05,0.05",
    lambda = "0"
  )

  allowed <- paste0("--", gsub("_", "-", names(values)))
  index <- 1L
  while (index <= length(args)) {
    option <- args[[index]]
    if (!option %in% allowed) {
      stop("unknown argument: ", option, call. = FALSE)
    }
    if (index == length(args)) {
      stop("missing value for argument: ", option, call. = FALSE)
    }
    key <- gsub("-", "_", sub("^--", "", option))
    values[[key]] <- args[[index + 1L]]
    index <- index + 2L
  }

  values$chains <- as.integer(values$chains)
  values$workers <- as.integer(values$workers)
  values$samples <- as.integer(values$samples)
  values$burn_in <- as.integer(values$burn_in)
  values$thin <- as.integer(values$thin)
  values$seed <- as.integer(values$seed)
  values$lambda <- as.numeric(values$lambda)

  integer_options <- c("chains", "workers", "samples", "burn_in", "thin")
  for (option in integer_options) {
    if (is.na(values[[option]]) || values[[option]] < 1L) {
      if (option == "burn_in" && identical(values[[option]], 0L)) next
      stop("--", gsub("_", "-", option), " must be a positive integer", call. = FALSE)
    }
  }
  if (is.na(values$seed)) stop("--seed must be an integer", call. = FALSE)
  if (!is.finite(values$lambda) || values$lambda < 0) {
    stop("--lambda must be a nonnegative number", call. = FALSE)
  }

  values$biomarker_cols <- trimws(strsplit(values$biomarker_cols, ",", fixed = TRUE)[[1L]])
  values$biomarker_labels <- trimws(strsplit(values$biomarker_labels, ",", fixed = TRUE)[[1L]])
  n_biomarkers <- length(values$biomarker_cols)
  if (length(values$biomarker_labels) != n_biomarkers) {
    stop("--biomarker-labels must have one label per biomarker", call. = FALSE)
  }
  values$gamma_start <- parse_numeric_vector(values$gamma_start, "--gamma-start", n_biomarkers)
  values$gamma_prior_mean <- parse_numeric_vector(values$gamma_prior_mean, "--gamma-prior-mean")
  values$gamma_prior_sd <- parse_numeric_vector(values$gamma_prior_sd, "--gamma-prior-sd")
  values$gamma_proposal_sd <- parse_numeric_vector(values$gamma_proposal_sd, "--gamma-proposal-sd")
  for (field in c("gamma_prior_mean", "gamma_prior_sd", "gamma_proposal_sd")) {
    if (!length(values[[field]]) %in% c(1L, n_biomarkers)) {
      stop("--", gsub("_", "-", field), " must have one or ", n_biomarkers, " values", call. = FALSE)
    }
  }
  if (any(values$gamma_prior_sd <= 0) || any(values$gamma_proposal_sd < 0)) {
    stop("prior SDs must be positive and proposal SDs must be nonnegative", call. = FALSE)
  }

  values$beta_start <- if (identical(tolower(values$beta_start), "auto")) {
    NULL
  } else {
    parse_numeric_vector(values$beta_start, "--beta-start", 3L)
  }
  values$reference_gamma <- if (tolower(values$reference_gamma) %in% c("none", "na")) {
    NULL
  } else {
    parse_numeric_vector(values$reference_gamma, "--reference-gamma", n_biomarkers)
  }
  values$reference_beta <- if (tolower(values$reference_beta) %in% c("none", "na")) {
    NULL
  } else {
    parse_numeric_vector(values$reference_beta, "--reference-beta", 3L)
  }
  values
}

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1L]]) else "scripts/run-biomarker-mcmc.R"
project_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

source_package_r_files <- function(root) {
  r_files <- list.files(file.path(root, "R"), pattern = "[.][rR]$", full.names = TRUE)
  for (r_file in r_files) source(r_file)
}

resolve_path <- function(path, root) {
  if (grepl("^/", path)) path else file.path(root, path)
}

load_analysis_data <- function(path, object_name) {
  extension <- tolower(tools::file_ext(path))
  if (!file.exists(path)) stop("data file does not exist: ", path, call. = FALSE)
  if (extension %in% c("rda", "rdata")) {
    environment <- new.env(parent = emptyenv())
    loaded <- load(path, envir = environment)
    if (object_name %in% loaded) return(environment[[object_name]])
    data_frames <- loaded[vapply(loaded, function(name) is.data.frame(environment[[name]]), logical(1))]
    if (length(data_frames) == 1L) return(environment[[data_frames]])
    stop("could not uniquely select a data frame; use --data-object", call. = FALSE)
  }
  if (extension == "rds") return(readRDS(path))
  if (extension == "csv") return(utils::read.csv(path, check.names = FALSE))
  stop("supported input formats are .rda, .RData, .rds, and .csv", call. = FALSE)
}

prepare_data <- function(data, args) {
  required <- c(args$time_col, args$status_col, args$treatment_col, args$biomarker_cols)
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) {
    stop("missing required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  analysis <- data[, required, drop = FALSE]
  names(analysis)[seq_len(3L)] <- c("time", "status", "treatment")

  convert_numeric <- function(value, name) {
    if (is.factor(value)) value <- as.character(value)
    converted <- suppressWarnings(as.numeric(value))
    if (anyNA(converted) && !anyNA(value)) stop(name, " could not be converted to numeric", call. = FALSE)
    converted
  }
  for (column in names(analysis)) analysis[[column]] <- convert_numeric(analysis[[column]], column)
  complete <- stats::complete.cases(analysis) & apply(is.finite(as.matrix(analysis)), 1L, all)
  removed <- sum(!complete)
  analysis <- analysis[complete, , drop = FALSE]
  if (!nrow(analysis)) stop("no complete finite observations remain", call. = FALSE)
  if (any(analysis$time <= 0)) stop("survival times must be positive", call. = FALSE)
  if (!all(analysis$status %in% c(0, 1))) stop("status must contain only 0 and 1", call. = FALSE)
  if (!all(analysis$treatment %in% c(0, 1))) stop("treatment must contain only 0 and 1", call. = FALSE)
  list(data = analysis, removed = removed)
}

initial_beta <- function(data, gamma_start, biomarker_cols) {
  y <- survival::Surv(data$time, data$status)
  sorted <- sort(data$time, decreasing = TRUE, index.return = TRUE)
  x <- data[sorted$ix, c("treatment", biomarker_cols), drop = FALSE]
  initialize_beta(
    x = x,
    y = y[sorted$ix, ],
    gamma = gamma_start,
    treatment_col = "treatment",
    biomarker_cols = biomarker_cols
  )
}

split_rhat <- function(chains) {
  chains <- as.matrix(chains)
  if (ncol(chains) < 2L || nrow(chains) < 4L) return(NA_real_)
  half <- floor(nrow(chains) / 2L)
  split <- do.call(cbind, lapply(seq_len(ncol(chains)), function(index) {
    cbind(chains[seq_len(half), index], chains[seq.int(nrow(chains) - half + 1L, nrow(chains)), index])
  }))
  within <- mean(apply(split, 2L, stats::var))
  if (!is.finite(within) || within <= 0) return(NA_real_)
  between <- half * stats::var(colMeans(split))
  variance <- (half - 1) / half * within + between / half
  sqrt(variance / within)
}

effective_sample_size <- function(chains) {
  chains <- as.matrix(chains)
  n <- nrow(chains)
  m <- ncol(chains)
  if (n < 3L) return(NA_real_)
  max_lag <- min(n - 1L, 1000L)
  autocorrelations <- vapply(seq_len(m), function(index) {
    as.numeric(stats::acf(chains[, index], lag.max = max_lag, plot = FALSE, na.action = na.pass)$acf)[-1L]
  }, numeric(max_lag))
  if (m == 1L) autocorrelations <- matrix(autocorrelations, ncol = 1L)
  if (any(!is.finite(autocorrelations))) return(NA_real_)
  rho <- rowMeans(autocorrelations)
  pair_count <- floor(length(rho) / 2L)
  if (pair_count) {
    pair_sums <- rho[2L * seq_len(pair_count) - 1L] + rho[2L * seq_len(pair_count)]
    first_nonpositive <- match(TRUE, pair_sums <= 0)
    if (!is.na(first_nonpositive)) pair_sums <- pair_sums[seq_len(first_nonpositive - 1L)]
    correlation_time <- 1 + 2 * sum(pair_sums)
  } else {
    correlation_time <- 1
  }
  min(m * n, m * n / max(correlation_time, 1))
}

posterior_arrays <- function(posteriors, beta_names, gamma_names) {
  parameter_names <- c(beta_names, gamma_names)
  samples <- lapply(posteriors, function(posterior) {
    values <- rbind(posterior$beta_samples, posterior$gamma_samples)
    rownames(values) <- parameter_names
    t(values)
  })
  list(names = parameter_names, chains = samples)
}

make_summary <- function(arrays, starts, references, labels) {
  do.call(rbind, lapply(seq_along(arrays$names), function(index) {
    chain_matrix <- do.call(cbind, lapply(arrays$chains, function(chain) chain[, index]))
    values <- as.vector(chain_matrix)
    interval <- stats::quantile(values, c(0.025, 0.975), names = FALSE)
    reference <- references[index]
    data.frame(
      parameter = arrays$names[index],
      label = labels[index],
      initial_value = starts[index],
      reference_value = reference,
      estimate = mean(values),
      posterior_sd = stats::sd(values),
      empirical_bias = if (is.na(reference)) NA_real_ else mean(values) - reference,
      interval_95_lower = interval[1],
      interval_95_upper = interval[2],
      interval_95_width = diff(interval),
      interval_covers_reference = if (is.na(reference)) NA else interval[1] <= reference && reference <= interval[2],
      split_rhat = split_rhat(chain_matrix),
      effective_sample_size = effective_sample_size(chain_matrix),
      stringsAsFactors = FALSE
    )
  }))
}

make_trace_data <- function(arrays, control) {
  rows <- vector("list", length(arrays$names) * length(arrays$chains))
  row_index <- 1L
  for (parameter_index in seq_along(arrays$names)) {
    for (chain_index in seq_along(arrays$chains)) {
      values <- arrays$chains[[chain_index]][, parameter_index]
      rows[[row_index]] <- data.frame(
        chain = chain_index,
        retained_sample = seq_along(values),
        iteration = control$burn_in + seq_along(values) * control$thin,
        parameter = arrays$names[parameter_index],
        value = values,
        running_mean = cumsum(values) / seq_along(values)
      )
      row_index <- row_index + 1L
    }
  }
  do.call(rbind, rows)
}

draw_traces <- function(trace, parameter_names, references) {
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(3L, 2L), mar = c(3.2, 3.8, 2.2, 1), oma = c(0, 0, 2, 0))
  colors <- rep(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00"),
    length.out = length(unique(trace$chain)))
  for (parameter_index in seq_along(parameter_names)) {
    subset <- trace[trace$parameter == parameter_names[parameter_index], , drop = FALSE]
    limits <- range(subset$value, references[parameter_index], finite = TRUE)
    graphics::plot(NA, xlim = range(subset$iteration), ylim = limits,
      xlab = "MCMC iteration", ylab = "value", main = parameter_names[parameter_index])
    for (chain in unique(subset$chain)) {
      chain_data <- subset[subset$chain == chain, ]
      graphics::lines(chain_data$iteration, chain_data$value,
        col = grDevices::adjustcolor(colors[chain], alpha.f = 0.55), lwd = 0.6)
    }
    if (is.finite(references[parameter_index])) {
      graphics::abline(h = references[parameter_index], lty = 2L, lwd = 2L, col = "black")
    }
  }
  graphics::mtext("MCMC convergence traces", outer = TRUE, font = 2)
}

draw_acf <- function(arrays) {
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(3L, 2L), mar = c(3.2, 3.8, 2.2, 1), oma = c(0, 0, 2, 0))
  colors <- rep(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00"),
    length.out = length(arrays$chains))
  lag_max <- min(100L, nrow(arrays$chains[[1L]]) - 1L)
  for (parameter_index in seq_along(arrays$names)) {
    graphics::plot(NA, xlim = c(0, lag_max), ylim = c(-0.2, 1),
      xlab = "retained-sample lag", ylab = "autocorrelation", main = arrays$names[parameter_index])
    graphics::abline(h = 0, col = "gray70")
    for (chain_index in seq_along(arrays$chains)) {
      acf_values <- stats::acf(arrays$chains[[chain_index]][, parameter_index], lag.max = lag_max, plot = FALSE)
      graphics::lines(as.numeric(acf_values$lag), as.numeric(acf_values$acf), col = colors[chain_index], lwd = 1.2)
    }
  }
  graphics::mtext("Within-chain autocorrelation", outer = TRUE, font = 2)
}

save_plot <- function(file, format, draw, width = 12, height = 12) {
  if (format == "png") grDevices::png(file, width = width, height = height, units = "in", res = 150)
  if (format == "pdf") grDevices::pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw()
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
source_package_r_files(project_root)
if (!requireNamespace("survival", quietly = TRUE) || !requireNamespace("mvtnorm", quietly = TRUE)) {
  stop("the survival and mvtnorm R packages must be installed", call. = FALSE)
}

data_path <- resolve_path(args$data, project_root)
output_dir <- resolve_path(args$output_dir, project_root)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
loaded <- load_analysis_data(data_path, args$data_object)
prepared <- prepare_data(loaded, args)
analysis_data <- prepared$data

beta_names <- c("beta_treatment", "beta_subgroup", "beta_treatment_x_subgroup")
gamma_names <- paste0("gamma_", args$biomarker_cols)
beta_labels <- c("Treatment", "Biomarker-defined subgroup", "Treatment x subgroup")
gamma_labels <- paste0(args$biomarker_labels, " (", args$biomarker_cols, ")")

beta_start <- args$beta_start
if (is.null(beta_start)) {
  beta_start <- initial_beta(analysis_data, args$gamma_start, args$biomarker_cols)
  message("Cox-derived beta start: ", paste(format(beta_start, digits = 7), collapse = ", "))
}
names(beta_start) <- beta_names
names(args$gamma_start) <- gamma_names

control <- default_mcmc_control(
  samples = args$samples,
  burn_in = args$burn_in,
  thin = args$thin,
  gamma_mean = args$gamma_prior_mean,
  gamma_sd = args$gamma_prior_sd,
  gamma_proposal_sd = args$gamma_proposal_sd
)

fit_chain <- function(chain_index) {
  set.seed(args$seed + chain_index - 1L)
  message("starting chain ", chain_index)
  posterior <- suppressWarnings(fit_threshold_model(
    data = analysis_data,
    control = control,
    lambda = args$lambda,
    gamma_start = args$gamma_start,
    time_col = "time",
    status_col = "status",
    treatment_col = "treatment",
    biomarker_cols = args$biomarker_cols,
    beta_start = beta_start
  ))
  message("finished chain ", chain_index)
  posterior
}

workers <- min(args$workers, args$chains)
if (.Platform$OS.type == "unix" && workers > 1L) {
  posteriors <- parallel::mclapply(seq_len(args$chains), fit_chain, mc.cores = workers, mc.preschedule = FALSE)
} else {
  posteriors <- lapply(seq_len(args$chains), fit_chain)
}
failed <- vapply(posteriors, inherits, logical(1), what = "try-error")
if (any(failed)) stop("one or more MCMC chains failed", call. = FALSE)

arrays <- posterior_arrays(posteriors, beta_names, gamma_names)
starts <- c(beta_start, args$gamma_start)
references <- c(
  if (is.null(args$reference_beta)) rep(NA_real_, 3L) else args$reference_beta,
  if (is.null(args$reference_gamma)) rep(NA_real_, length(args$gamma_start)) else args$reference_gamma
)
labels <- c(beta_labels, gamma_labels)
summary <- make_summary(arrays, starts, references, labels)
trace <- make_trace_data(arrays, control)

utils::write.csv(summary, file.path(output_dir, "parameter-summary.csv"), row.names = FALSE, na = "")
trace_connection <- gzfile(file.path(output_dir, "posterior-traces.csv.gz"), open = "wt")
utils::write.csv(trace, trace_connection, row.names = FALSE)
close(trace_connection)
saveRDS(posteriors, file.path(output_dir, "posterior-chains.rds"))

settings <- data.frame(
  data_file = normalizePath(data_path),
  data_object = args$data_object,
  observations_used = nrow(analysis_data),
  incomplete_rows_removed = prepared$removed,
  events = sum(analysis_data$status),
  treated = sum(analysis_data$treatment == 1),
  chains = args$chains,
  workers = workers,
  retained_samples_per_chain = args$samples,
  burn_in = args$burn_in,
  thin = args$thin,
  seed = args$seed,
  lambda = args$lambda,
  stringsAsFactors = FALSE
)
utils::write.csv(settings, file.path(output_dir, "run-settings.csv"), row.names = FALSE)

start_values <- data.frame(
  parameter = c(beta_names, gamma_names),
  label = labels,
  initial_value = starts
)
utils::write.csv(start_values, file.path(output_dir, "initial-values.csv"), row.names = FALSE)

for (format in c("png", "pdf")) {
  save_plot(file.path(output_dir, paste0("convergence-traces.", format)), format,
    function() draw_traces(trace, arrays$names, references))
  save_plot(file.path(output_dir, paste0("convergence-autocorrelation.", format)), format,
    function() draw_acf(arrays))
}

notes <- c(
  "parameter-summary.csv contains posterior means and equal-tailed 95% Bayesian credible intervals.",
  "empirical_bias and interval_covers_reference are blank unless --reference-beta and/or --reference-gamma are supplied.",
  "Starting values are algorithm inputs, not true values, and therefore are not used as bias/coverage references.",
  "For a simulation study, aggregate empirical bias and coverage over independent generated datasets, not only MCMC chains.",
  "Inspect split_rhat (ideally near 1), effective_sample_size, trace plots, and autocorrelation plots before interpreting estimates."
)
writeLines(notes, file.path(output_dir, "README-results.txt"))
capture.output(utils::sessionInfo(), file = file.path(output_dir, "session-info.txt"))

print(summary)
message("results written to ", normalizePath(output_dir))
