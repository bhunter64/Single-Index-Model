#!/usr/bin/env Rscript

parse_args <- function(args) {
  values <- list(
    data = if (file.exists("/data/cleaned_data.rda")) "/data/cleaned_data.rda" else "data/cleaned_data.rda",
    data_object = "cleaned_data",
    three_summary = NA_character_,
    pairwise_summary = NA_character_,
    output_dir = "results/km-from-mcmc",
    time_col = "time",
    status_col = "status",
    treatment_col = "trt"
  )
  allowed <- paste0("--", gsub("_", "-", names(values)))

  index <- 1L
  while (index <= length(args)) {
    option <- args[[index]]
    if (identical(option, "--help")) {
      cat(
        "Usage: Rscript scripts/run-km-from-mcmc-results.R [options]\n",
        "  --three-summary PATH     Three-biomarker parameter-summary.csv\n",
        "  --pairwise-summary PATH  Pairwise pairwise-parameter-summary.csv\n",
        "  --data PATH              .rda, .rds, or .csv survival data\n",
        "  --data-object NAME       Object in an .rda file (default cleaned_data)\n",
        "  --output-dir PATH        Output directory\n",
        "  --time-col NAME          Survival-time column (default time)\n",
        "  --status-col NAME        Event column (default status)\n",
        "  --treatment-col NAME     Binary treatment column (default trt)\n",
        sep = ""
      )
      quit(status = 0L)
    }
    if (!option %in% allowed) stop("unknown argument: ", option, call. = FALSE)
    if (index == length(args)) stop("missing value for argument: ", option, call. = FALSE)
    key <- gsub("-", "_", sub("^--", "", option))
    values[[key]] <- args[[index + 1L]]
    index <- index + 2L
  }

  if (is.na(values$three_summary) && is.na(values$pairwise_summary)) {
    stop("provide --three-summary, --pairwise-summary, or both", call. = FALSE)
  }
  values
}

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1L]]) else "scripts/run-km-from-mcmc-results.R"
project_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

resolve_path <- function(path) {
  if (is.na(path) || grepl("^/", path)) path else file.path(project_root, path)
}

load_analysis_data <- function(path, object_name) {
  if (!file.exists(path)) stop("data file does not exist: ", path, call. = FALSE)
  extension <- tolower(tools::file_ext(path))
  if (extension %in% c("rda", "rdata")) {
    environment <- new.env(parent = emptyenv())
    loaded <- load(path, envir = environment)
    if (object_name %in% loaded) return(environment[[object_name]])
    frames <- loaded[vapply(loaded, function(name) is.data.frame(environment[[name]]), logical(1))]
    if (length(frames) == 1L) return(environment[[frames]])
    stop("could not uniquely select a data frame; use --data-object", call. = FALSE)
  }
  if (extension == "rds") return(readRDS(path))
  if (extension == "csv") return(utils::read.csv(path, check.names = FALSE))
  stop("supported data formats are .rda, .RData, .rds, and .csv", call. = FALSE)
}

read_summary <- function(path) {
  if (!file.exists(path)) stop("summary file does not exist: ", path, call. = FALSE)
  summary <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("parameter", "label", "estimate")
  missing <- setdiff(required, names(summary))
  if (length(missing)) stop("summary is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  summary
}

model_specification <- function(name, summary) {
  gamma_rows <- grepl("^gamma_", summary$parameter)
  if (!any(gamma_rows)) stop("no gamma rows found for model ", name, call. = FALSE)
  gamma <- summary[gamma_rows, , drop = FALSE]
  biomarkers <- sub("^gamma_", "", gamma$parameter)
  estimates <- as.numeric(gamma$estimate)
  if (anyNA(estimates) || any(!is.finite(estimates))) {
    stop("gamma estimates are not finite for model ", name, call. = FALSE)
  }
  list(
    name = name,
    biomarkers = biomarkers,
    labels = gamma$label,
    gamma = stats::setNames(estimates, biomarkers)
  )
}

load_model_specifications <- function(three_path, pairwise_path) {
  specifications <- list()
  if (!is.na(three_path)) {
    summary <- read_summary(three_path)
    specification <- model_specification("HER2_CA19-9_Axl", summary)
    specifications[[specification$name]] <- specification
  }
  if (!is.na(pairwise_path)) {
    summary <- read_summary(pairwise_path)
    if (!"biomarker_pair" %in% names(summary)) {
      stop("pairwise summary must contain biomarker_pair", call. = FALSE)
    }
    for (name in unique(summary$biomarker_pair)) {
      specification <- model_specification(name, summary[summary$biomarker_pair == name, , drop = FALSE])
      specifications[[specification$name]] <- specification
    }
  }
  specifications
}

as_binary_numeric <- function(value, name) {
  if (is.factor(value)) value <- as.character(value)
  converted <- suppressWarnings(as.numeric(value))
  if (anyNA(converted) || !all(converted %in% c(0, 1))) {
    stop(name, " must contain only numeric 0 and 1 values", call. = FALSE)
  }
  converted
}

prepare_model_data <- function(data, specification, args) {
  required <- c(args$time_col, args$status_col, args$treatment_col, specification$biomarkers)
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("data is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  selected <- data[, required, drop = FALSE]
  names(selected)[seq_len(3L)] <- c("time", "status", "treatment")

  if (is.factor(selected$time)) selected$time <- as.character(selected$time)
  selected$time <- suppressWarnings(as.numeric(selected$time))
  selected$status <- as_binary_numeric(selected$status, "status")
  selected$treatment <- as_binary_numeric(selected$treatment, "treatment")
  for (biomarker in specification$biomarkers) {
    value <- selected[[biomarker]]
    if (is.factor(value)) value <- as.character(value)
    selected[[biomarker]] <- suppressWarnings(as.numeric(value))
  }

  finite <- stats::complete.cases(selected) & apply(is.finite(as.matrix(selected)), 1L, all)
  selected <- selected[finite, , drop = FALSE]
  if (!nrow(selected)) stop("no finite complete observations for model ", specification$name, call. = FALSE)
  if (any(selected$time <= 0)) stop("survival times must be positive", call. = FALSE)

  selected$score <- as.vector(1 + as.matrix(selected[, specification$biomarkers, drop = FALSE]) %*% specification$gamma)
  selected$subgroup_numeric <- as.integer(selected$score >= 0)
  selected$subgroup <- factor(
    selected$subgroup_numeric,
    levels = c(0, 1),
    labels = c("Below threshold", "Above threshold")
  )
  selected$treatment_label <- factor(
    selected$treatment,
    levels = c(0, 1),
    labels = c("Control", "Treated")
  )
  selected$km_group <- interaction(selected$treatment_label, selected$subgroup, sep = " / ", drop = TRUE)
  selected
}

safe_logrank <- function(data, group, comparison) {
  group <- droplevels(factor(group))
  if (nlevels(group) < 2L) {
    return(data.frame(
      comparison = comparison,
      n = nrow(data),
      groups = nlevels(group),
      chisq = NA_real_,
      df = NA_integer_,
      p_value = NA_real_
    ))
  }
  fit <- survival::survdiff(survival::Surv(data$time, data$status) ~ group)
  degrees_freedom <- length(fit$n) - 1L
  data.frame(
    comparison = comparison,
    n = nrow(data),
    groups = length(fit$n),
    chisq = unname(fit$chisq),
    df = degrees_freedom,
    p_value = stats::pchisq(fit$chisq, degrees_freedom, lower.tail = FALSE)
  )
}

interaction_results <- function(data) {
  combination_counts <- table(data$treatment, data$subgroup_numeric)
  if (
    length(unique(data$treatment)) < 2L ||
      length(unique(data$subgroup_numeric)) < 2L ||
      any(combination_counts == 0L)
  ) {
    return(data.frame(
      interaction_estimable = FALSE,
      beta3_estimate = NA_real_,
      beta3_standard_error = NA_real_,
      beta3_hazard_ratio = NA_real_,
      beta3_hr_95_lower = NA_real_,
      beta3_hr_95_upper = NA_real_,
      beta3_wald_z = NA_real_,
      beta3_wald_p = NA_real_,
      interaction_likelihood_ratio_p = NA_real_,
      proportional_hazards_global_p = NA_real_
    ))
  }
  response <- survival::Surv(data$time, data$status)
  additive <- survival::coxph(response ~ treatment + subgroup_numeric, data = data, ties = "efron")
  interaction <- survival::coxph(response ~ treatment * subgroup_numeric, data = data, ties = "efron")
  coefficient_name <- "treatment:subgroup_numeric"
  coefficient_table <- summary(interaction)$coefficients
  confidence <- stats::confint(interaction, coefficient_name)
  comparison <- stats::anova(additive, interaction, test = "LRT")
  p_column <- grep("P|Pr", names(comparison), value = TRUE)
  likelihood_ratio_p <- if (length(p_column)) as.numeric(comparison[2L, p_column[1L]]) else NA_real_
  proportional_hazards <- tryCatch(survival::cox.zph(interaction), error = function(error) NULL)
  global_ph_p <- if (
    is.null(proportional_hazards) ||
      !"GLOBAL" %in% rownames(proportional_hazards$table)
  ) NA_real_ else unname(proportional_hazards$table["GLOBAL", "p"])

  estimate <- unname(coefficient_table[coefficient_name, "coef"])
  standard_error <- unname(coefficient_table[coefficient_name, "se(coef)"])
  data.frame(
    interaction_estimable = TRUE,
    beta3_estimate = estimate,
    beta3_standard_error = standard_error,
    beta3_hazard_ratio = exp(estimate),
    beta3_hr_95_lower = exp(confidence[1L]),
    beta3_hr_95_upper = exp(confidence[2L]),
    beta3_wald_z = estimate / standard_error,
    beta3_wald_p = unname(coefficient_table[coefficient_name, "Pr(>|z|)"]),
    interaction_likelihood_ratio_p = likelihood_ratio_p,
    proportional_hazards_global_p = global_ph_p
  )
}

group_summary <- function(data) {
  groups <- split(data, data$km_group, drop = TRUE)
  do.call(rbind, lapply(names(groups), function(name) {
    group <- groups[[name]]
    data.frame(
      group = name,
      n = nrow(group),
      events = sum(group$status),
      censored = sum(group$status == 0),
      median_time = unname(summary(survival::survfit(survival::Surv(time, status) ~ 1, data = group))$table["median"]),
      stringsAsFactors = FALSE
    )
  }))
}

format_p <- function(value) {
  if (!is.finite(value)) "p unavailable" else if (value < 0.001) "p < 0.001" else paste0("p = ", formatC(value, digits = 3, format = "f"))
}

draw_km_with_risk_table <- function(data, group, title, subtitle = NULL) {
  if (!nrow(data)) {
    layout(1)
    graphics::par(mar = c(4, 4, 4, 2))
    graphics::plot.new()
    graphics::title(main = title)
    graphics::text(0.5, 0.55, "No observations in this fitted-score subgroup", cex = 1.2)
    if (!is.null(subtitle)) graphics::text(0.5, 0.43, subtitle, cex = 0.9)
    return(invisible(NULL))
  }
  group <- droplevels(factor(group))
  plot_data <- data.frame(time = data$time, status = data$status, group = group)
  fit <- survival::survfit(survival::Surv(time, status) ~ group, data = plot_data)
  group_names <- levels(group)
  colors <- rep(c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), length.out = length(group_names))
  line_types <- seq_along(group_names)
  maximum_time <- max(plot_data$time)
  time_points <- pretty(c(0, maximum_time), n = 6L)
  time_points <- time_points[time_points >= 0 & time_points <= maximum_time]

  layout(matrix(c(1, 2), nrow = 2L), heights = c(3.8, 1.35))
  graphics::par(mar = c(3.2, 4.5, 4.2, 1.5))
  graphics::plot(
    fit,
    col = colors,
    lty = line_types,
    lwd = 2.2,
    conf.int = FALSE,
    mark.time = TRUE,
    xlab = "Time",
    ylab = "Survival probability",
    main = title,
    xlim = c(0, maximum_time),
    ylim = c(0, 1)
  )
  if (!is.null(subtitle)) graphics::mtext(subtitle, side = 3, line = 0.4, cex = 0.85)
  graphics::legend("bottomleft", legend = group_names, col = colors, lty = line_types, lwd = 2.2, bty = "n", cex = 0.85)

  graphics::par(mar = c(3.5, 11, 0.5, 1.5))
  graphics::plot.new()
  graphics::plot.window(xlim = c(-0.28 * maximum_time, maximum_time), ylim = c(0.5, length(group_names) + 1.15), xaxs = "i")
  graphics::axis(1, at = time_points)
  graphics::mtext("Time", side = 1, line = 2.2)
  graphics::text(-0.27 * maximum_time, length(group_names) + 0.8, "Number at risk", adj = c(0, 0.5), xpd = NA, font = 2)
  for (group_index in seq_along(group_names)) {
    group_data <- plot_data[plot_data$group == group_names[group_index], , drop = FALSE]
    risks <- vapply(time_points, function(time) sum(group_data$time >= time), integer(1))
    y <- length(group_names) - group_index + 1L
    graphics::text(-0.02 * maximum_time, y, group_names[group_index], adj = c(1, 0.5), xpd = NA, col = colors[group_index], cex = 0.78)
    graphics::text(time_points, y, labels = risks, cex = 0.78)
  }
  layout(1)
}

save_plot <- function(path, format, draw) {
  if (format == "png") grDevices::png(path, width = 12, height = 8.5, units = "in", res = 180)
  if (format == "pdf") grDevices::pdf(path, width = 12, height = 8.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw()
}

run_model <- function(data, specification, args, output_root) {
  message("Analyzing ", specification$name)
  model_data <- prepare_model_data(data, specification, args)
  model_dir <- file.path(output_root, specification$name)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

  four_group <- safe_logrank(model_data, model_data$km_group, "Four treatment-by-subgroup groups")
  below <- model_data[model_data$subgroup_numeric == 0L, , drop = FALSE]
  above <- model_data[model_data$subgroup_numeric == 1L, , drop = FALSE]
  tests <- rbind(
    four_group,
    safe_logrank(below, below$treatment_label, "Treatment within below-threshold subgroup"),
    safe_logrank(above, above$treatment_label, "Treatment within above-threshold subgroup")
  )
  interaction <- interaction_results(model_data)
  counts <- group_summary(model_data)
  gamma <- data.frame(
    biomarker = specification$biomarkers,
    label = specification$labels,
    gamma_estimate = unname(specification$gamma),
    score_definition = "1 + biomarkers %*% gamma; above if score >= 0",
    stringsAsFactors = FALSE
  )

  utils::write.csv(model_data, file.path(model_dir, "classified-data.csv"), row.names = FALSE)
  utils::write.csv(counts, file.path(model_dir, "group-summary.csv"), row.names = FALSE)
  utils::write.csv(tests, file.path(model_dir, "logrank-tests.csv"), row.names = FALSE)
  utils::write.csv(interaction, file.path(model_dir, "cox-interaction-test.csv"), row.names = FALSE)
  utils::write.csv(gamma, file.path(model_dir, "gamma-score-definition.csv"), row.names = FALSE)

  interaction_text <- if (isTRUE(interaction$interaction_estimable)) {
    paste0("beta-3 interaction Wald ", format_p(interaction$beta3_wald_p))
  } else {
    "beta-3 interaction not estimable"
  }
  four_group_subtitle <- paste0(
    "Four-group log-rank ", format_p(four_group$p_value),
    "; ", interaction_text
  )
  below_test <- tests$p_value[tests$comparison == "Treatment within below-threshold subgroup"]
  above_test <- tests$p_value[tests$comparison == "Treatment within above-threshold subgroup"]
  for (format in c("png", "pdf")) {
    save_plot(file.path(model_dir, paste0("km-four-groups.", format)), format, function() {
      draw_km_with_risk_table(model_data, model_data$km_group, paste0(specification$name, ": treatment and biomarker subgroup"), four_group_subtitle)
    })
    save_plot(file.path(model_dir, paste0("km-treatment-below-threshold.", format)), format, function() {
      draw_km_with_risk_table(below, below$treatment_label, paste0(specification$name, ": below threshold"), paste0("Treatment log-rank ", format_p(below_test)))
    })
    save_plot(file.path(model_dir, paste0("km-treatment-above-threshold.", format)), format, function() {
      draw_km_with_risk_table(above, above$treatment_label, paste0(specification$name, ": above threshold"), paste0("Treatment log-rank ", format_p(above_test)))
    })
  }

  list(
    logrank = data.frame(model = specification$name, tests, check.names = FALSE),
    interaction = data.frame(model = specification$name, n = nrow(model_data), interaction, check.names = FALSE),
    groups = data.frame(model = specification$name, counts, check.names = FALSE),
    gamma = data.frame(model = specification$name, gamma, check.names = FALSE)
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!requireNamespace("survival", quietly = TRUE)) stop("the survival R package must be installed", call. = FALSE)
args$data <- resolve_path(args$data)
args$three_summary <- resolve_path(args$three_summary)
args$pairwise_summary <- resolve_path(args$pairwise_summary)
args$output_dir <- resolve_path(args$output_dir)

specifications <- load_model_specifications(args$three_summary, args$pairwise_summary)
data <- load_analysis_data(args$data, args$data_object)
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
results <- lapply(specifications, function(specification) run_model(data, specification, args, args$output_dir))

utils::write.csv(do.call(rbind, lapply(results, `[[`, "logrank")), file.path(args$output_dir, "all-logrank-tests.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, lapply(results, `[[`, "interaction")), file.path(args$output_dir, "all-cox-interaction-tests.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, lapply(results, `[[`, "groups")), file.path(args$output_dir, "all-group-summaries.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, lapply(results, `[[`, "gamma")), file.path(args$output_dir, "all-gamma-score-definitions.csv"), row.names = FALSE)

notes <- c(
  "Subgroups use the fitted posterior-mean gamma values and the package score definition: 1 + biomarkers %*% gamma.",
  "Scores greater than or equal to zero are classified as above threshold.",
  "Kaplan-Meier curves and log-rank tests describe survival differences but do not directly test beta 3.",
  "The beta-3 Wald and likelihood-ratio p-values come from Cox models with treatment, subgroup, and treatment-by-subgroup interaction.",
  "Because gamma and subgroup membership were estimated from these same observations, all p-values are exploratory/post-selection and may be anti-conservative.",
  "A confirmatory analysis requires a prespecified score evaluated in independent validation data or a method that propagates subgroup-estimation uncertainty."
)
writeLines(notes, file.path(args$output_dir, "README-KM-results.txt"))
capture.output(utils::sessionInfo(), file = file.path(args$output_dir, "session-info.txt"))

message("KM analysis complete: ", normalizePath(args$output_dir))
