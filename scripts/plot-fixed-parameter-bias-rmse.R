#!/usr/bin/env Rscript

# Plot mean absolute bias and RMSE for the fixed-parameter simulation.

suppressPackageStartupMessages(library(ggplot2))

sample_sizes <- c(600, 800, 1000, 1200, 1400)

# Values transcribed from cov_probs_and_biases.tex and Table appFixedParameter.
simulation_error <- data.frame(
  sample_size = rep(sample_sizes, times = 10),
  parameter = rep(
    rep(c("beta[1]", "beta[2]", "beta[3]", "gamma[1]", "gamma[2]"), each = 5),
    times = 2
  ),
  metric = rep(c("Mean absolute bias", "RMSE"), each = 25),
  value = c(
    0.1246964, 0.1165518, 0.09935666, 0.08897241, 0.08252119,
    0.1160485, 0.1041444, 0.08871557, 0.07805910, 0.07054184,
    0.1531722, 0.13899059, 0.12450412, 0.10620848, 0.09941074,
    0.1757917, 0.1503266, 0.1190034, 0.1016028, 0.09259884,
    0.1936907, 0.1629776, 0.1219236, 0.1103291, 0.09452585,
    0.1576896, 0.1461301, 0.1258421, 0.11069432, 0.10368373,
    0.1444895, 0.1305509, 0.1110494, 0.09900556, 0.08826997,
    0.1919513, 0.1723994, 0.1562753, 0.13264493, 0.12415340,
    0.2369967, 0.2091948, 0.1600997, 0.1414929, 0.1265932,
    0.2576818, 0.2212743, 0.1661922, 0.1552442, 0.1331465
  )
)

simulation_error$parameter <- factor(
  simulation_error$parameter,
  levels = c("beta[1]", "beta[2]", "beta[3]", "gamma[1]", "gamma[2]")
)
simulation_error$metric <- factor(
  simulation_error$metric,
  levels = c("Mean absolute bias", "RMSE")
)

figure <- ggplot(
  simulation_error,
  aes(
    x = sample_size,
    y = value,
    colour = metric,
    linetype = metric,
    shape = metric,
    group = metric
  )
) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.8, stroke = 0.45) +
  facet_wrap(~parameter, ncol = 3, labeller = label_parsed, axes = "all_x") +
  scale_x_continuous(breaks = c(600, 1000, 1400)) +
  scale_y_continuous(
    limits = c(0.06, 0.27),
    breaks = seq(0.08, 0.24, by = 0.04),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_colour_manual(
    values = c("Mean absolute bias" = "#0072B2", "RMSE" = "#D55E00")
  ) +
  scale_linetype_manual(
    values = c("Mean absolute bias" = "solid", "RMSE" = "dashed")
  ) +
  scale_shape_manual(
    values = c("Mean absolute bias" = 16, "RMSE" = 15)
  ) +
  labs(
    x = "Sample size, n",
    y = "Estimation error",
    colour = NULL,
    linetype = NULL,
    shape = NULL
  ) +
  theme_bw(base_family = "serif", base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    panel.border = element_rect(linewidth = 0.35),
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8, colour = "black"),
    legend.position = "bottom",
    legend.box.margin = margin(t = -2),
    legend.key.width = grid::unit(1.6, "lines"),
    plot.margin = margin(4, 5, 3, 3)
  )

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
project_root <- normalizePath(file.path(dirname(script_path), ".."))
output_dir <- file.path(project_root, "results", "fixed-parameter-simulation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pdf_path <- file.path(output_dir, "simulation_bias_rmse.pdf")
png_path <- file.path(output_dir, "simulation_bias_rmse.png")

ggsave(pdf_path, figure, width = 7.2, height = 4.8, units = "in", device = "pdf")
ggsave(png_path, figure, width = 7.2, height = 4.8, units = "in", dpi = 300)

writeLines(c(pdf_path, png_path))
