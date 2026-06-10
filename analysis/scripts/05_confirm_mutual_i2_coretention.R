
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])
  script_path <- chartr("\\", "/", script_path)
} else {
  script_path <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

SCRIPT_DIR <- if (file.exists(script_path)) dirname(script_path) else script_path
PROJECT_DIR_LONG <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), winslash = "/", mustWork = TRUE)
PROJECT_DIR <- if (.Platform$OS.type == "windows") utils::shortPathName(PROJECT_DIR_LONG) else PROJECT_DIR_LONG
OUT_TABLES <- file.path(PROJECT_DIR, "analysis", "outputs", "tables")
OUT_LOGS <- file.path(PROJECT_DIR, "analysis", "outputs", "logs")
dir.create(OUT_TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_LOGS, recursive = TRUE, showWarnings = FALSE)

set.seed(20260605L + 35L)
N_PERM <- as.integer(Sys.getenv("WIDER_I2_N_PERM", unset = "5000"))
if (!is.finite(N_PERM) || N_PERM < 1000L) {
  stop("WIDER_I2_N_PERM must be at least 1000 for confirmation.", call. = FALSE)
}

I2_GROUPS <- c(
  "PH908_primary",
  "PH908_map_nonprimary",
  "Y3120_deep_mid_nonPH908",
  "Y3120_shallow_nonPH908",
  "I2_unmapped_nonPH908"
)

is_i2_unit <- function(d) {
  d$phylo_group %in% I2_GROUPS
}

empirical_p_greater <- function(null_vals, obs) {
  null_vals <- null_vals[is.finite(null_vals)]
  (1 + sum(null_vals >= obs)) / (length(null_vals) + 1)
}

ci95 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  as.numeric(stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
}

did_detail <- function(target, relic, yes) {
  n_tr <- sum(target & relic)
  n_tf <- sum(target & !relic)
  n_cr <- sum(!target & relic)
  n_cf <- sum(!target & !relic)
  if (min(n_tr, n_tf, n_cr, n_cf) < 5) return(NULL)
  y_tr <- sum(yes[target & relic], na.rm = TRUE)
  y_tf <- sum(yes[target & !relic], na.rm = TRUE)
  y_cr <- sum(yes[!target & relic], na.rm = TRUE)
  y_cf <- sum(yes[!target & !relic], na.rm = TRUE)
  r_tr <- y_tr / n_tr
  r_tf <- y_tf / n_tf
  r_cr <- y_cr / n_cr
  r_cf <- y_cf / n_cf
  tibble(
    target_relic_rate = r_tr,
    target_founder_rate = r_tf,
    control_relic_rate = r_cr,
    control_founder_rate = r_cf,
    target_relic_n = n_tr,
    target_founder_n = n_tf,
    control_relic_n = n_cr,
    control_founder_n = n_cf,
    target_relic_yes = y_tr,
    target_founder_yes = y_tf,
    control_relic_yes = y_cr,
    control_founder_yes = y_cf,
    target_relic_minus_founder = r_tr - r_tf,
    control_relic_minus_founder = r_cr - r_cf,
    did_lift = (r_tr - r_tf) - (r_cr - r_cf)
  )
}

did_value <- function(target, relic, yes) {
  st <- did_detail(target, relic, yes)
  if (is.null(st)) return(NA_real_)
  st$did_lift
}

make_region_size_bin <- function(region) {
  counts <- table(region)
  n <- as.numeric(counts[region])
  qs <- unique(stats::quantile(n, probs = c(0, 0.33, 0.66, 1), na.rm = TRUE))
  if (length(qs) < 2) return(rep("bin_1", length(region)))
  out <- paste0("bin_", cut(n, breaks = qs, include.lowest = TRUE, labels = FALSE))
  out[is.na(out)] <- "bin_1"
  out
}

permute_logical_by_strata <- function(x, split_idx) {
  out <- x
  for (idx in split_idx) {
    if (length(idx) > 1) out[idx] <- sample(out[idx], length(idx), replace = FALSE)
  }
  out
}

null_target <- function(target, relic, yes, strata, n_perm) {
  split_idx <- split(seq_along(strata), as.character(strata))
  vals <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    vals[i] <- did_value(permute_logical_by_strata(target, split_idx), relic, yes)
  }
  vals
}

null_class <- function(target, relic, yes, strata, n_perm) {
  split_idx <- split(seq_along(strata), as.character(strata))
  vals <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    vals[i] <- did_value(target, permute_logical_by_strata(relic, split_idx), yes)
  }
  vals
}

prepare_pair <- function(d, control_set) {
  target_filter <- is_i2_unit(d)
  control_filter <- if (control_set == "non_I2_excl_R1a") {
    d$broad_family == "non_I2" & d$phylo_group != "R1a"
  } else if (control_set == "all_non_I2") {
    d$broad_family == "non_I2"
  } else {
    stop("Unknown control set.", call. = FALSE)
  }

  d |>
    mutate(
      .is_target_candidate = target_filter,
      .is_control_candidate = control_filter,
      .metric_yes = other_i2_cloud_social_spatial
    ) |>
    filter(
      (.is_target_candidate | .is_control_candidate),
      branch_class %in% c("Relic", "Founder"),
      !is.na(.metric_yes)
    ) |>
    mutate(is_target = .is_target_candidate) |>
    select(-.is_target_candidate, -.is_control_candidate)
}

observed_row <- function(d, unit_type, control_set, n_perm) {
  pair <- prepare_pair(d, control_set) |>
    filter(!is.na(Region), Region != "")
  target <- as.logical(pair$is_target)
  relic <- pair$branch_class == "Relic"
  yes <- as.logical(pair$.metric_yes)
  region <- as.character(pair$Region)
  region_size <- make_region_size_bin(region)
  st <- did_detail(target, relic, yes)
  if (is.null(st)) return(NULL)

  n1 <- null_target(target, relic, yes, paste(region, relic, sep = " | "), n_perm)
  n2 <- null_target(target, relic, yes, paste(region, region_size, relic, sep = " | "), n_perm)
  n3 <- null_class(target, relic, yes, paste(target, region, sep = " | "), n_perm)
  q1 <- ci95(n1)
  q2 <- ci95(n2)
  q3 <- ci95(n3)

  st |>
    mutate(
      variant = "high_resolution_upstream_excluded",
      unit_type = unit_type,
      target_mode = "I2_all_mutual",
      control_set = control_set,
      metric = "other_i2_cloud_social_spatial",
      n_perm = n_perm,
      target_region_branch_perm_p = empirical_p_greater(n1, did_lift),
      target_region_size_branch_perm_p = empirical_p_greater(n2, did_lift),
      class_family_region_perm_p = empirical_p_greater(n3, did_lift),
      target_region_branch_ci_low = q1[1],
      target_region_branch_ci_high = q1[2],
      target_region_size_branch_ci_low = q2[1],
      target_region_size_branch_ci_high = q2[2],
      class_family_region_ci_low = q3[1],
      class_family_region_ci_high = q3[2]
    )
}

leave_one_region_rows <- function(d, unit_type, control_set) {
  pair <- prepare_pair(d, control_set) |>
    filter(!is.na(Region), Region != "")
  bind_rows(lapply(sort(unique(pair$Region)), function(reg) {
    x <- pair |> filter(Region != reg)
    st <- did_detail(x$is_target, x$branch_class == "Relic", as.logical(x$.metric_yes))
    if (is.null(st)) return(NULL)
    st |> mutate(unit_type = unit_type, control_set = control_set, omitted_region = reg)
  }))
}

leave_one_target_rows <- function(d, unit_type, control_set, factor_name) {
  pair <- prepare_pair(d, control_set) |>
    filter(!is.na(Region), Region != "")
  target_levels <- pair |>
    filter(is_target) |>
    distinct(.drop_level = .data[[factor_name]]) |>
    filter(!is.na(.drop_level), .drop_level != "") |>
    pull(.drop_level)
  bind_rows(lapply(sort(target_levels), function(level) {
    x <- pair |> filter(!(is_target & .data[[factor_name]] == level))
    st <- did_detail(x$is_target, x$branch_class == "Relic", as.logical(x$.metric_yes))
    if (is.null(st)) return(NULL)
    st |> mutate(
      unit_type = unit_type,
      control_set = control_set,
      factor_name = factor_name,
      omitted_level = level
    )
  }))
}

units_path <- file.path(OUT_TABLES, "mutual_i2_coretention_units.csv")
if (!file.exists(units_path)) {
  stop("Run 04_build_mutual_i2_units.R before this step.", call. = FALSE)
}
units <- readr::read_csv(units_path, show_col_types = FALSE) |>
  filter(variant == "high_resolution_upstream_excluded")

unit_types <- c("cluster", "surname_region_terminal", "location_terminal")
control_sets <- c("non_I2_excl_R1a", "all_non_I2")
target_factors <- c("phylo_group", "terminal_snp_norm", "Slava", "Region")

confirmed <- bind_rows(lapply(unit_types, function(ut) {
  d <- units |> filter(unit_type == ut)
  bind_rows(lapply(control_sets, function(cs) {
    cat("Confirming ", ut, " / ", cs, "\n", sep = "")
    observed_row(d, ut, cs, N_PERM)
  }))
})) |>
  select(variant, unit_type, target_mode, control_set, metric, n_perm, everything()) |>
  arrange(target_region_branch_perm_p, target_region_size_branch_perm_p, class_family_region_perm_p)

region_delete <- bind_rows(lapply(unit_types, function(ut) {
  d <- units |> filter(unit_type == ut)
  bind_rows(lapply(control_sets, function(cs) leave_one_region_rows(d, ut, cs)))
}))

target_delete <- bind_rows(lapply(unit_types, function(ut) {
  d <- units |> filter(unit_type == ut)
  bind_rows(lapply(control_sets, function(cs) {
    bind_rows(lapply(target_factors, function(f) leave_one_target_rows(d, ut, cs, f)))
  }))
}))

robust <- confirmed |>
  left_join(
    region_delete |>
      group_by(unit_type, control_set) |>
      summarise(
        min_leave_one_region_did = min(did_lift, na.rm = TRUE),
        weakest_region = omitted_region[which.min(did_lift)],
        .groups = "drop"
      ),
    by = c("unit_type", "control_set")
  ) |>
  left_join(
    target_delete |>
      group_by(unit_type, control_set, factor_name) |>
      summarise(
        min_leave_one_target_factor_did = min(did_lift, na.rm = TRUE),
        weakest_target_level = omitted_level[which.min(did_lift)],
        .groups = "drop"
      ) |>
      tidyr::pivot_wider(
        names_from = factor_name,
        values_from = c(min_leave_one_target_factor_did, weakest_target_level),
        names_glue = "{.value}_{factor_name}"
      ),
    by = c("unit_type", "control_set")
  ) |>
  arrange(target_region_branch_perm_p, target_region_size_branch_perm_p, class_family_region_perm_p)

group_detail <- units |>
  filter(
    unit_type %in% unit_types,
    phylo_group %in% I2_GROUPS,
    branch_class %in% c("Relic", "Founder")
  ) |>
  group_by(unit_type, phylo_group, branch_class) |>
  summarise(
    n = n(),
    yes = sum(other_i2_cloud_social_spatial, na.rm = TRUE),
    rate = yes / n,
    terminal_snps = n_distinct(terminal_snp_norm),
    slavas = n_distinct(Slava),
    regions = n_distinct(Region),
    .groups = "drop"
  )

readr::write_csv(confirmed, file.path(OUT_TABLES, "mutual_i2_coretention_confirmation_5000.csv"))
readr::write_csv(region_delete, file.path(OUT_TABLES, "mutual_i2_coretention_confirmation_leave_region.csv"))
readr::write_csv(target_delete, file.path(OUT_TABLES, "mutual_i2_coretention_confirmation_leave_target.csv"))
readr::write_csv(robust, file.path(OUT_TABLES, "mutual_i2_coretention_confirmation_robust.csv"))
readr::write_csv(group_detail, file.path(OUT_TABLES, "mutual_i2_coretention_confirmation_group_detail.csv"))

best_lines <- robust |>
  mutate(line = paste0(
    "- ", unit_type, " / ", control_set,
    ": target relic=", target_relic_yes, "/", target_relic_n,
    " (", sprintf("%.3f", target_relic_rate), ")",
    ", target founder=", target_founder_yes, "/", target_founder_n,
    " (", sprintf("%.3f", target_founder_rate), ")",
    ", control relic=", control_relic_yes, "/", control_relic_n,
    " (", sprintf("%.3f", control_relic_rate), ")",
    ", control founder=", control_founder_yes, "/", control_founder_n,
    " (", sprintf("%.3f", control_founder_rate), ")",
    ", DID=", sprintf("%.3f", did_lift),
    ", target-region-branch p=", sprintf("%.5f", target_region_branch_perm_p),
    ", target-region-size-branch p=", sprintf("%.5f", target_region_size_branch_perm_p),
    ", class-family-region p=", sprintf("%.5f", class_family_region_perm_p),
    ", min leave-one-region DID=", sprintf("%.3f", min_leave_one_region_did),
    ", min leave-one-phylo DID=", sprintf("%.3f", min_leave_one_target_factor_did_phylo_group),
    ", min leave-one-SNP DID=", sprintf("%.3f", min_leave_one_target_factor_did_terminal_snp_norm),
    ", min leave-one-Slava DID=", sprintf("%.3f", min_leave_one_target_factor_did_Slava),
    ", min leave-one-target-region DID=", sprintf("%.3f", min_leave_one_target_factor_did_Region)
  )) |>
  pull(line)

writeLines(
  c(
    paste0("Completed 35_mutual_i2_coretention_confirmation.R at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("N_PERM=", N_PERM),
    paste0("Rows=", nrow(confirmed))
  ),
  file.path(OUT_LOGS, "35_mutual_i2_coretention_confirmation.log")
)
