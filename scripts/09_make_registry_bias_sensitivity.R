#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  script_path <- sub("^--file=", "", file_arg[1])
  if (!grepl("^[A-Za-z]:|^/", script_path)) {
    root <- if (.Platform$OS.type == "windows") utils::shortPathName(getwd()) else getwd()
    script_path <- file.path(root, script_path)
  }
  SCRIPT_DIR <- dirname(script_path)
  PAPER_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = TRUE)
} else {
  PAPER_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
EXEC_DIR <- if (.Platform$OS.type == "windows") utils::shortPathName(PAPER_DIR) else PAPER_DIR
TABLE_DIR <- file.path(EXEC_DIR, "result_tables")
FIG_DIR <- file.path(EXEC_DIR, "figures")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(20260609L)
N_SIM <- as.integer(Sys.getenv("WIDER_I2_REGISTRY_BIAS_N_SIM", unset = "5000"))
if (!is.finite(N_SIM) || N_SIM < 500L) {
  stop("WIDER_I2_REGISTRY_BIAS_N_SIM must be at least 500.", call. = FALSE)
}

BIAS_LEVELS <- seq(0, 0.50, by = 0.025)
PRIMARY_UNIT <- "surname_region_terminal"
PRIMARY_CONTROL <- "non_I2_excl_R1a"
PRIMARY_METRICS <- c(
  "own_family_social_xsurname_xlocation",
  "own_family_cloud_xsurname_xlocation"
)

metric_labels <- c(
  own_family_social_xsurname_xlocation = "Shared Slava or region",
  own_family_cloud_xsurname_xlocation = "Coordinate cloud"
)

scenario_labels <- c(
  random_registry_noise = "Non-directional registry noise",
  directed_adversarial_bias = "Coordinated adverse bias"
)

theme_wider <- function(base_size = 10.5) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "#172033"),
      plot.title = element_text(face = "bold", size = base_size + 3, margin = margin(b = 6)),
      plot.subtitle = element_text(size = base_size, color = "#4B5563", margin = margin(b = 10)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#243447"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", color = "#172033"),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.margin = margin(12, 16, 12, 16)
    )
}

wrap_lab <- function(x, width = 95) {
  paste(strwrap(x, width = width), collapse = "\n")
}

save_both <- function(plot, basename, width = 8.4, height = 6.2) {
  png_path <- file.path(FIG_DIR, paste0(basename, ".png"))
  pdf_path <- file.path(FIG_DIR, paste0(basename, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 320, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, bg = "white")
  c(png = png_path, pdf = pdf_path)
}

read_required <- function(name) {
  path <- file.path(TABLE_DIR, name)
  if (!file.exists(path)) stop("Missing required table: ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE)
}

did_from_counts <- function(y_tr, n_tr, y_tf, n_tf, y_cr, n_cr, y_cf, n_cf) {
  (y_tr / n_tr - y_tf / n_tf) - (y_cr / n_cr - y_cf / n_cf)
}

simulate_counts <- function(row, bias_level, scenario) {
  bias_level <- as.numeric(bias_level[1])
  scenario <- as.character(scenario[1])

  y_tr <- as.integer(row$target_relic_yes)
  n_tr <- as.integer(row$target_relic_n)
  y_tf <- as.integer(row$target_founder_yes)
  n_tf <- as.integer(row$target_founder_n)
  y_cr <- as.integer(row$control_relic_yes)
  n_cr <- as.integer(row$control_relic_n)
  y_cf <- as.integer(row$control_founder_yes)
  n_cf <- as.integer(row$control_founder_n)

  if (scenario == "random_registry_noise") {
    y_tr <- rbinom(1, y_tr, 1 - bias_level) + rbinom(1, n_tr - y_tr, bias_level)
    y_tf <- rbinom(1, y_tf, 1 - bias_level) + rbinom(1, n_tf - y_tf, bias_level)
    y_cr <- rbinom(1, y_cr, 1 - bias_level) + rbinom(1, n_cr - y_cr, bias_level)
    y_cf <- rbinom(1, y_cf, 1 - bias_level) + rbinom(1, n_cf - y_cf, bias_level)
  } else if (scenario == "directed_adversarial_bias") {
    y_tr <- rbinom(1, y_tr, 1 - bias_level)
    y_tf <- y_tf + rbinom(1, n_tf - y_tf, bias_level)
    y_cr <- y_cr + rbinom(1, n_cr - y_cr, bias_level)
    y_cf <- rbinom(1, y_cf, 1 - bias_level)
  } else {
    stop("Unknown scenario: ", scenario, call. = FALSE)
  }

  did_from_counts(y_tr, n_tr, y_tf, n_tf, y_cr, n_cr, y_cf, n_cf)
}

crossing_level <- function(curve, value_col = "median_did", threshold = 0) {
  hit <- curve |>
    filter(.data[[value_col]] <= threshold) |>
    arrange(bias_level) |>
    slice_head(n = 1)
  if (nrow(hit) == 0) return(NA_real_)
  hit$bias_level[1]
}

results <- read_required("own_family_reciprocal_results.csv") |>
  filter(
    unit_type == PRIMARY_UNIT,
    control_set == PRIMARY_CONTROL,
    metric %in% PRIMARY_METRICS
  ) |>
  arrange(match(metric, PRIMARY_METRICS))

if (nrow(results) != length(PRIMARY_METRICS)) {
  stop("Could not find both primary wider-I2 metrics in own_family_reciprocal_results.csv.", call. = FALSE)
}

observed <- results |>
  transmute(
    unit_type,
    control_set,
    metric,
    metric_label = metric_labels[metric],
    observed_did = did_lift,
    target_relic = paste0(target_relic_yes, "/", target_relic_n),
    target_founder = paste0(target_founder_yes, "/", target_founder_n),
    control_relic = paste0(control_relic_yes, "/", control_relic_n),
    control_founder = paste0(control_founder_yes, "/", control_founder_n)
  )

draws <- bind_rows(lapply(seq_len(nrow(results)), function(i) {
  row <- results[i, ]
  bind_rows(lapply(names(scenario_labels), function(scenario) {
    bind_rows(lapply(BIAS_LEVELS, function(bias_level) {
      tibble(
        metric = row$metric,
        metric_label = metric_labels[row$metric],
        scenario = scenario,
        scenario_label = scenario_labels[scenario],
        bias_level = bias_level,
        sim_id = seq_len(N_SIM),
        did_lift = replicate(N_SIM, simulate_counts(row, bias_level, scenario))
      )
    }))
  }))
}))

curves <- draws |>
  group_by(metric, metric_label, scenario, scenario_label, bias_level) |>
  summarise(
    n_sim = n(),
    median_did = median(did_lift, na.rm = TRUE),
    p05 = quantile(did_lift, 0.05, na.rm = TRUE),
    p95 = quantile(did_lift, 0.95, na.rm = TRUE),
    prob_positive = mean(did_lift > 0, na.rm = TRUE),
    .groups = "drop"
  )

tipping <- curves |>
  group_by(metric, metric_label, scenario, scenario_label) |>
  group_modify(~ tibble(
    median_zero_bias = crossing_level(.x, "median_did", 0),
    p05_zero_bias = crossing_level(.x, "p05", 0),
    did_at_10pct = .x$median_did[which.min(abs(.x$bias_level - 0.10))],
    did_at_15pct = .x$median_did[which.min(abs(.x$bias_level - 0.15))],
    did_at_20pct = .x$median_did[which.min(abs(.x$bias_level - 0.20))],
    positive_prob_at_10pct = .x$prob_positive[which.min(abs(.x$bias_level - 0.10))],
    positive_prob_at_15pct = .x$prob_positive[which.min(abs(.x$bias_level - 0.15))],
    positive_prob_at_20pct = .x$prob_positive[which.min(abs(.x$bias_level - 0.20))]
  )) |>
  ungroup() |>
  mutate(
    median_zero_label = ifelse(
      is.na(median_zero_bias),
      "not reached by 50%",
      percent(median_zero_bias, accuracy = 0.1)
    )
  )

readr::write_csv(observed, file.path(TABLE_DIR, "registry_bias_sensitivity_observed.csv"))
readr::write_csv(curves, file.path(TABLE_DIR, "registry_bias_sensitivity_curves.csv"))
readr::write_csv(tipping, file.path(TABLE_DIR, "registry_bias_sensitivity_tipping.csv"))

curves_plot <- curves |>
  mutate(
    metric_label = factor(metric_label, levels = metric_labels),
    scenario_label = factor(scenario_label, levels = scenario_labels)
  )

scenario_colors <- c(
  "Non-directional registry noise" = "#64748B",
  "Coordinated adverse bias" = "#9B2226"
)

p_curves <- ggplot(
  curves_plot,
  aes(x = bias_level, y = median_did, color = scenario_label, fill = scenario_label)
) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#9B2226", alpha = 0.07) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#172033", linewidth = 0.45) +
  geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.13, color = NA) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 1.35) +
  facet_wrap(~ metric_label, ncol = 2) +
  scale_x_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 0.50, 0.10)) +
  scale_y_continuous(labels = label_number(accuracy = 0.05)) +
  coord_cartesian(ylim = c(-0.36, 0.42)) +
  scale_color_manual(values = scenario_colors) +
  scale_fill_manual(values = scenario_colors) +
  labs(
    title = "Registry-bias stress test",
    subtitle = wrap_lab("Random noise flips outcomes symmetrically; coordinated adverse bias recodes all four DID cells against the I2 result.", 96),
    x = "Simulated registry-bias level",
    y = "Median DID"
  ) +
  theme_wider() +
  theme(panel.grid.major.y = element_line(color = "#E5E7EB", linewidth = 0.35))

tipping_plot <- tipping |>
  mutate(
    metric_label = factor(metric_label, levels = metric_labels),
    scenario_label = factor(scenario_label, levels = rev(scenario_labels)),
    plot_bias = ifelse(is.na(median_zero_bias), 0.50, median_zero_bias),
    threshold_label = ifelse(
      is.na(median_zero_bias) | plot_bias >= 0.50,
      "\u226550%",
      percent(median_zero_bias, accuracy = 0.1)
    ),
    label_hjust = ifelse(plot_bias >= 0.49, 1.12, -0.18)
  )

p_tip <- ggplot(tipping_plot, aes(x = plot_bias, y = scenario_label, color = scenario_label)) +
  geom_segment(aes(x = 0, xend = plot_bias, yend = scenario_label), linewidth = 0.9, alpha = 0.55) +
  geom_point(size = 3.1) +
  geom_text(aes(label = threshold_label, hjust = label_hjust), size = 3.1, fontface = "bold", color = "#172033") +
  facet_wrap(~ metric_label, ncol = 2) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.56), breaks = seq(0, 0.50, 0.10)) +
  scale_color_manual(values = scenario_colors, guide = "none") +
  labs(
    title = "Bias level needed to erase the median DID",
    subtitle = "Random noise reaches zero only at the edge of the 50% noise envelope; directed bias is a stronger adversarial test.",
    x = "Median-zero tipping point",
    y = NULL
  ) +
  theme_wider(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text = element_blank()
  )

figure <- p_curves / p_tip +
  patchwork::plot_layout(heights = c(1.25, 0.75)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 11, color = "#172033"))
  )

save_both(figure, "Figure_7_registry_bias_sensitivity", width = 8.4, height = 6.4)

manifest_path <- file.path(FIG_DIR, "figure_manifest.csv")
existing_manifest <- if (file.exists(manifest_path)) {
  readr::read_csv(manifest_path, show_col_types = FALSE)
} else {
  tibble(figure = character(), png = character(), pdf = character())
}

new_manifest <- tibble(
  figure = "Figure_7_registry_bias_sensitivity",
  png = file.path("figures", "Figure_7_registry_bias_sensitivity.png"),
  pdf = file.path("figures", "Figure_7_registry_bias_sensitivity.pdf")
)

manifest <- existing_manifest |>
  filter(figure != "Figure_7_registry_bias_sensitivity") |>
  bind_rows(new_manifest)

readr::write_csv(manifest, manifest_path)

message("Wrote registry-bias sensitivity tables and Figure 7 to: ", PAPER_DIR)
