#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(scales)
})

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(script_path)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else getwd()
}

PAPER_DIR <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
EXEC_DIR <- if (.Platform$OS.type == "windows") utils::shortPathName(PAPER_DIR) else PAPER_DIR
TABLE_DIR <- file.path(EXEC_DIR, "result_tables")
FIG_DIR <- file.path(EXEC_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

read_required <- function(name) {
  path <- file.path(TABLE_DIR, name)
  if (!file.exists(path)) stop("Missing required table: ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE)
}

save_both <- function(plot, basename, width = 7.2, height = 5.0) {
  png_path <- file.path(FIG_DIR, paste0(basename, ".png"))
  pdf_path <- file.path(FIG_DIR, paste0(basename, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 320, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, bg = "white")
  c(png = png_path, pdf = pdf_path)
}

metric_labels <- c(
  own_family_social_xsurname_xlocation = "Shared Slava or region\nsame surname/location excluded",
  own_family_cloud_xsurname_xlocation = "Coordinate cloud\nsame surname/location excluded"
)

family_palette <- c(
  I2 = "#005F73",
  J2 = "#CA6702",
  G = "#7B7F35",
  R1a = "#8D99AE",
  E = "#6D6875",
  R1b = "#9B2226",
  I1 = "#495057"
)

theme_wider <- function(base_size = 10.5) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "#172033"),
      plot.title = element_text(face = "bold", size = base_size + 3, margin = margin(b = 6)),
      plot.subtitle = element_text(size = base_size, color = "#4B5563", margin = margin(b = 10)),
      plot.caption = element_text(size = base_size - 2, color = "#6B7280", hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#243447"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", color = "#172033"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.margin = margin(12, 16, 12, 16)
    )
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p))))
}

wrap_lab <- function(x, width = 95) {
  paste(strwrap(x, width = width), collapse = "\n")
}

evidence_summary <- read_required("reciprocal_i2_evidence_summary.csv")
family_summary <- read_required("own_family_reciprocal_family_summary.csv")
own_results <- read_required("own_family_reciprocal_results.csv")
polarity_obs <- read_required("exact_region_polarity_obs.csv")
polarity_best <- read_required("signal_class_polarity_summary.csv")

# Figure 1: conceptual design -------------------------------------------------
box_df <- tibble::tibble(
  step = c(
    "Lineage families",
    "Own relic field",
    "Downweighted links",
    "Reciprocal contrast",
    "Ranked result"
  ),
  x = c(0.5, 2.1, 3.7, 5.3, 6.9),
  y = 1,
  fill = c("#E0F2FE", "#E9F5DB", "#FFF1D6", "#EDE9FE", "#D9F99D"),
  text = c(
    "I2, J2, G,\nR1a, R1b, E, I1",
    "Each family scored\nagainst its own older\nfamily field",
    "Same surname and\nsame location removed",
    "Relic branches vs\nfounder branches",
    "I2 ranks first;\nJ2 follows;\nother families lower"
  )
)

arrow_df <- tibble::tibble(
  x = box_df$x[-nrow(box_df)] + 0.63,
  xend = box_df$x[-1] - 0.63,
  y = 1,
  yend = 1
)

fig1 <- ggplot() +
  geom_rect(
    data = box_df,
    aes(xmin = x - 0.62, xmax = x + 0.62, ymin = y - 0.34, ymax = y + 0.34, fill = fill),
    color = "#334155",
    linewidth = 0.55
  ) +
  geom_segment(
    data = arrow_df,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.16, "in")),
    linewidth = 0.65,
    color = "#334155"
  ) +
  geom_text(data = box_df, aes(x = x, y = y + 0.11, label = step), fontface = "bold", size = 3.0) +
  geom_text(data = box_df, aes(x = x, y = y - 0.08, label = text), size = 2.65, lineheight = 0.95) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(-0.15, 7.55), ylim = c(0.45, 1.55), expand = FALSE) +
  labs(
    title = wrap_lab("The reciprocal own-family test asks whether relic branches retain deeper family neighborhoods", 78),
    subtitle = wrap_lab("The design makes I2 compete against other lineage families on the same logic, while removing direct surname and locality carry.", 90)
  ) +
  theme_void(base_size = 10.5) +
  theme(
    text = element_text(color = "#172033"),
    plot.title = element_text(face = "bold", size = 13, margin = margin(b = 6)),
    plot.subtitle = element_text(size = 10, color = "#4B5563", margin = margin(b = 10)),
    plot.caption = element_text(size = 8.5, color = "#6B7280", hjust = 0),
    plot.margin = margin(18, 22, 18, 22)
  )

save_both(fig1, "Figure_1_own_family_test_design", width = 8.4, height = 3.45)

# Figure 2: lineage-family ranking ------------------------------------------
family_plot_data <- family_summary |>
  filter(
    unit_type == "cluster",
    metric %in% names(metric_labels)
  ) |>
  mutate(
    metric_label = factor(metric_labels[metric], levels = metric_labels),
    reservoir_family = factor(reservoir_family, levels = rev(c("I2", "J2", "G", "R1a", "E", "R1b", "I1"))),
    value_label = ifelse(abs(relic_minus_founder) < 0.0005, "0.000", sprintf("%+.3f", relic_minus_founder)),
    family_fill = ifelse(reservoir_family %in% names(family_palette), as.character(reservoir_family), "Other")
  )

fig2 <- ggplot(family_plot_data, aes(x = relic_minus_founder, y = reservoir_family, fill = reservoir_family)) +
  geom_vline(xintercept = 0, color = "#94A3B8", linewidth = 0.7) +
  geom_col(width = 0.66, color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = value_label),
    hjust = ifelse(family_plot_data$relic_minus_founder >= 0, -0.12, 1.12),
    size = 3.0,
    color = "#172033"
  ) +
  facet_wrap(~ metric_label, ncol = 2) +
  scale_fill_manual(values = family_palette, guide = "none") +
  scale_x_continuous(labels = label_number(accuracy = 0.01), limits = c(-0.36, 0.34)) +
  labs(
    title = wrap_lab("I2 is the top-ranked family when each retained lineage is judged against its own relic field", 74),
    subtitle = wrap_lab("Positive values mean relic-class branches retain more own-family structure than founder-class branches after same surname and same location are excluded.", 86),
    x = "Relic minus founder retention rate",
    y = NULL
  ) +
  theme_wider()

save_both(fig2, "Figure_2_lineage_family_ranking", width = 7.6, height = 5.2)

# Figure 3: anatomy of the primary DID ---------------------------------------
primary_rows <- own_results |>
  filter(
    unit_type == "surname_region_terminal",
    control_set == "non_I2_excl_R1a",
    metric %in% names(metric_labels)
  ) |>
  select(unit_type, control_set, metric, did_lift,
         target_relic_rate, target_founder_rate, control_relic_rate, control_founder_rate,
         target_relic_yes, target_relic_n, target_founder_yes, target_founder_n,
         control_relic_yes, control_relic_n, control_founder_yes, control_founder_n,
         target_region_branch_p, target_region_size_branch_p, class_target_region_p) |>
  pivot_longer(
    cols = ends_with("_rate"),
    names_to = "series",
    values_to = "rate"
  ) |>
  mutate(
    family = ifelse(grepl("^target", series), "I2 target", "Non-I2 controls"),
    branch = ifelse(grepl("relic", series), "Relic class", "Founder class"),
    metric_label = factor(metric_labels[metric], levels = metric_labels),
    branch = factor(branch, levels = c("Founder class", "Relic class")),
    family = factor(family, levels = c("Non-I2 controls", "I2 target"))
  )

fig3 <- ggplot(primary_rows, aes(x = branch, y = rate, fill = family)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = percent(rate, accuracy = 0.1)),
    position = position_dodge(width = 0.72),
    vjust = -0.35,
    size = 3.0
  ) +
  facet_wrap(~ metric_label, ncol = 2) +
  scale_fill_manual(values = c("I2 target" = "#005F73", "Non-I2 controls" = "#64748B")) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.48), expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = wrap_lab("The signal is a contrast between relic retention and founder carry", 78),
    subtitle = wrap_lab("The I2 target gains sharply from founder to relic class; the non-I2 controls do not show the same shift.", 86),
    x = NULL,
    y = "Retention rate",
    fill = NULL
  ) +
  theme_wider()

save_both(fig3, "Figure_3_primary_did_anatomy", width = 7.4, height = 5.2)

# Figure 4: comparative evidence summary -------------------------------------
comparison_summary <- evidence_summary |>
  mutate(
    evidence_short = case_when(
      rank == 1 ~ "Reciprocal own-family social\ncross-surname/location",
      rank == 2 ~ "Reciprocal own-family coordinate\ncross-surname/location",
      rank == 3 ~ "Mutual I2 cloud\ncross-surname/location",
      TRUE ~ "Previous mutual I2\nsocial-spatial baseline"
    ),
    evidence_short = factor(evidence_short, levels = rev(evidence_short)),
    p_label = paste0(
      "base ", fmt_p(base_target_p),
      ifelse(!is.na(hard_target_slavafreq_p), paste0("; hard ", fmt_p(hard_target_slavafreq_p)), "")
    ),
    highlight = rank <= 2
  )

fig4 <- ggplot(comparison_summary, aes(x = did_lift, y = evidence_short, fill = highlight)) +
  geom_col(width = 0.66, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("DID = %.3f", did_lift)), hjust = -0.10, size = 3.1, fontface = "bold") +
  geom_text(aes(label = p_label), hjust = -0.06, nudge_y = -0.18, size = 2.45, color = "#4B5563") +
  scale_fill_manual(values = c("TRUE" = "#005F73", "FALSE" = "#94A3B8"), guide = "none") +
  scale_x_continuous(limits = c(0, 0.46), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = wrap_lab("The reciprocal own-family design exceeds earlier mutual-I2 co-retention baselines", 78),
    subtitle = wrap_lab("The leading result preserves the reservoir logic while adding symmetric family controls and stricter surname/location exclusions.", 86),
    x = "Difference-in-differences lift",
    y = NULL
  ) +
  theme_wider()

save_both(fig4, "Figure_4_comparative_evidence_summary", width = 8.4, height = 4.9)

# Figure 5: robustness heatmap -----------------------------------------------
robust <- evidence_summary |>
  filter(rank <= 3) |>
  select(rank, evidence,
         min_leave_one_region_did,
         min_leave_one_slava_did,
         min_leave_one_surname_did,
         min_leave_one_i2_subgroup_did) |>
  mutate(
    evidence_short = case_when(
      rank == 1 ~ "Reciprocal social",
      rank == 2 ~ "Reciprocal coordinate",
      TRUE ~ "Mutual I2 cross-location"
    )
  ) |>
  pivot_longer(
    cols = starts_with("min_leave_one"),
    names_to = "stress_test",
    values_to = "min_did"
  ) |>
  mutate(
    stress_test = recode(
      stress_test,
      min_leave_one_region_did = "Drop one region",
      min_leave_one_slava_did = "Drop one Slava",
      min_leave_one_surname_did = "Drop one surname",
      min_leave_one_i2_subgroup_did = "Drop one I2 subgroup"
    ),
    evidence_short = factor(evidence_short, levels = rev(c("Reciprocal social", "Reciprocal coordinate", "Mutual I2 cross-location"))),
    stress_test = factor(stress_test, levels = c("Drop one region", "Drop one Slava", "Drop one surname", "Drop one I2 subgroup"))
  )

fig5 <- ggplot(robust, aes(x = stress_test, y = evidence_short, fill = min_did)) +
  geom_tile(color = "white", linewidth = 1.1) +
  geom_text(aes(label = sprintf("%.3f", min_did)), size = 3.2, fontface = "bold", color = "#102A43") +
  scale_fill_gradientn(
    colors = c("#FEE2E2", "#FEF3C7", "#BFDBFE", "#0F766E"),
    values = rescale(c(0.15, 0.20, 0.27, 0.34)),
    limits = c(0.15, 0.34),
    name = "Minimum DID"
  ) +
  labs(
    title = wrap_lab("The leading result survives leave-one deletion of regions, Slavas, surnames, and I2 subgroups", 60),
    subtitle = wrap_lab("Cells show the weakest retained DID after deleting one category at a time.", 86),
    x = NULL,
    y = NULL
  ) +
  theme_wider() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

save_both(fig5, "Figure_5_adversarial_robustness", width = 7.4, height = 4.6)

# Figure 6: polarity bridge ---------------------------------------------------
polarity_plot <- polarity_obs |>
  filter(unit_type %in% c("surname_region_terminal", "location_terminal", "cluster")) |>
  select(
    unit_type,
    class_polarity_contrast,
    ph908_relic_core_share,
    ph908_founder_core_share,
    r1a_relic_core_share,
    r1a_founder_core_share
  ) |>
  pivot_longer(
    cols = ends_with("_core_share"),
    names_to = "group_class",
    values_to = "core_share"
  ) |>
  mutate(
    family = ifelse(grepl("^ph908", group_class), "PH908", "R1a"),
    branch_class = ifelse(grepl("relic", group_class), "Relic class", "Founder class"),
    unit_label = recode(
      unit_type,
      surname_region_terminal = "Surname-region-terminal",
      location_terminal = "Location-terminal",
      cluster = "Cluster"
    ),
    unit_label = factor(unit_label, levels = c("Surname-region-terminal", "Location-terminal", "Cluster")),
    branch_class = factor(branch_class, levels = c("Founder class", "Relic class")),
    family = factor(family, levels = c("PH908", "R1a"))
  )

polarity_p <- polarity_best |>
  filter(unit_type %in% c("surname_region_terminal", "location_terminal", "cluster")) |>
  mutate(
    unit_label = recode(
      unit_type,
      surname_region_terminal = "Surname-region-terminal",
      location_terminal = "Location-terminal",
      cluster = "Cluster"
    ),
    label = paste0("polarity DID ", sprintf("%.3f", observed), "\n", fmt_p(max_p)),
    y = 0.54
  )

fig6 <- ggplot(polarity_plot, aes(x = branch_class, y = core_share, fill = family)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = percent(core_share, accuracy = 0.1)),
    position = position_dodge(width = 0.72),
    vjust = -0.35,
    size = 2.8
  ) +
  geom_text(
    data = polarity_p,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.8,
    color = "#334155",
    fontface = "bold"
  ) +
  facet_wrap(~ unit_label, ncol = 3) +
  scale_fill_manual(values = c("PH908" = "#005F73", "R1a" = "#8D99AE")) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.62), expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = wrap_lab("The polarity bridge separates PH908 relic retention from R1a founder carry", 72),
    subtitle = wrap_lab("PH908 points toward relic-class core retention, while R1a does not reproduce the same relic polarity.", 88),
    x = NULL,
    y = "Core-field share",
    fill = NULL
  ) +
  theme_wider()

save_both(fig6, "Figure_6_polarity_bridge", width = 7.8, height = 4.9)

manifest <- tibble::tibble(
  figure = c(
    "Figure_1_own_family_test_design",
    "Figure_2_lineage_family_ranking",
    "Figure_3_primary_did_anatomy",
    "Figure_4_comparative_evidence_summary",
    "Figure_5_adversarial_robustness",
    "Figure_6_polarity_bridge"
  ),
  png = file.path("figures", paste0(figure, ".png")),
  pdf = file.path("figures", paste0(figure, ".pdf"))
)

readr::write_csv(manifest, file.path(FIG_DIR, "figure_manifest.csv"))
message("Wrote wider-I2 figures to: ", FIG_DIR)
