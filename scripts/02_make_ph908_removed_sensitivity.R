#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
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
DATA_DIR <- file.path(EXEC_DIR, "data")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(20260609L + 2L)
N_PERM <- as.integer(Sys.getenv("PH908_SENSITIVITY_N_PERM", unset = "5000"))
if (!is.finite(N_PERM) || N_PERM < 100L) {
  stop("PH908_SENSITIVITY_N_PERM must be at least 100.", call. = FALSE)
}

UNIT_TYPES <- c("surname_region_terminal")
CONTROL_SETS <- c("non_I2_excl_R1a")
METRICS <- c(
  "own_family_social_xsurname_xlocation",
  "own_family_cloud_xsurname_xlocation"
)
I2_GROUPS <- c(
  "PH908_primary",
  "PH908_map_nonprimary",
  "Y3120_deep_mid_nonPH908",
  "Y3120_shallow_nonPH908",
  "I2_unmapped_nonPH908"
)
PH908_GROUPS <- c("PH908_primary", "PH908_map_nonprimary")
FAMILIES <- c("I2", "R1a", "R1b", "E", "J2", "G", "I1")

deg2rad <- function(x) x * pi / 180

haversine_km <- function(lat1, lon1, lat2, lon2) {
  r <- 6371.0088
  p1 <- deg2rad(lat1)
  p2 <- deg2rad(lat2)
  dp <- deg2rad(lat2 - lat1)
  dl <- deg2rad(lon2 - lon1)
  a <- sin(dp / 2)^2 + cos(p1) * cos(p2) * sin(dl / 2)^2
  2 * r * atan2(sqrt(a), sqrt(1 - a))
}

family_label <- function(phylo_group) {
  ifelse(phylo_group %in% I2_GROUPS, "I2", as.character(phylo_group))
}

subgroup_label <- function(family, phylo_group, terminal_snp_norm) {
  ifelse(family == "I2", as.character(phylo_group), as.character(terminal_snp_norm))
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

permute_by_strata <- function(x, strata) {
  out <- x
  split_idx <- split(seq_along(strata), as.character(strata))
  for (idx in split_idx) {
    if (length(idx) > 1) out[idx] <- sample(out[idx], length(idx), replace = FALSE)
  }
  out
}

empirical_p_greater <- function(null_vals, obs) {
  null_vals <- null_vals[is.finite(null_vals)]
  (1 + sum(null_vals >= obs)) / (length(null_vals) + 1)
}

control_mask <- function(d, control_set, target) {
  if (control_set == "non_I2_excl_R1a") {
    d$broad_family == "non_I2" & d$phylo_group != "R1a" & !target
  } else if (control_set == "all_non_I2") {
    d$broad_family == "non_I2" & !target
  } else {
    stop("Unknown control set.", call. = FALSE)
  }
}

add_own_family_fields <- function(d) {
  d <- d |>
    mutate(
      reservoir_family = family_label(phylo_group),
      reservoir_subgroup = subgroup_label(reservoir_family, phylo_group, terminal_snp_norm),
      slava_region = paste(slava_id, region_id, sep = " | ")
    )

  refs <- d |>
    filter(
      reservoir_family %in% FAMILIES,
      branch_class == "Relic",
      is.finite(lat),
      is.finite(long)
    )

  rows <- lapply(seq_len(nrow(d)), function(i) {
    row <- d[i, ]
    ref <- refs |>
      filter(
        reservoir_family == row$reservoir_family,
        reservoir_subgroup != row$reservoir_subgroup,
        surname_id != row$surname_id,
        location_id != row$location_id
      )
    if (nrow(ref) < 2 || !is.finite(row$lat) || !is.finite(row$long)) {
      return(tibble(
        own_family_social_xsurname_xlocation = NA,
        own_family_cloud_xsurname_xlocation = NA,
        own_family_ref_n = nrow(ref)
      ))
    }
    dist <- haversine_km(row$lat, row$long, ref$lat, ref$long)
    same_cell <- ref$slava_region == row$slava_region
    near <- dist <= 30
    tibble(
      own_family_social_xsurname_xlocation = any(same_cell, na.rm = TRUE),
      own_family_cloud_xsurname_xlocation = any(same_cell, na.rm = TRUE) & any(near, na.rm = TRUE),
      own_family_ref_n = nrow(ref)
    )
  })

  bind_cols(d, bind_rows(rows))
}

prepare_pair <- function(d, control_set, metric_col) {
  target <- d$reservoir_family == "I2"
  control <- control_mask(d, control_set, target)
  d |>
    mutate(is_target = target, metric_yes = .data[[metric_col]]) |>
    filter(
      (is_target | control),
      branch_class %in% c("Relic", "Founder"),
      !is.na(metric_yes),
      !is.na(region_id),
      region_id != ""
    )
}

score_case <- function(d, unit_type, control_set, metric_col, n_perm) {
  pair <- prepare_pair(d, control_set, metric_col)
  target <- pair$is_target
  relic <- pair$branch_class == "Relic"
  yes <- as.logical(pair$metric_yes)
  st <- did_detail(target, relic, yes)
  if (is.null(st)) return(NULL)

  region_size <- make_region_size_bin(pair$region_id)
  strata_target <- paste(pair$region_id, pair$branch_class, sep = " | ")
  strata_target_size <- paste(pair$region_id, region_size, pair$branch_class, sep = " | ")
  strata_class <- paste(target, pair$region_id, sep = " | ")
  null_target <- numeric(n_perm)
  null_target_size <- numeric(n_perm)
  null_class <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    null_target[i] <- did_value(permute_by_strata(target, strata_target), relic, yes)
    null_target_size[i] <- did_value(permute_by_strata(target, strata_target_size), relic, yes)
    null_class[i] <- did_value(target, permute_by_strata(relic, strata_class), yes)
  }

  st |>
    mutate(
      unit_type = unit_type,
      control_set = control_set,
      metric = metric_col,
      n_perm = n_perm,
      ph908_removed = TRUE,
      ph908_groups_removed = paste(PH908_GROUPS, collapse = ";"),
      target_region_branch_p = empirical_p_greater(null_target, did_lift),
      target_region_size_branch_p = empirical_p_greater(null_target_size, did_lift),
      class_target_region_p = empirical_p_greater(null_class, did_lift)
    )
}

source_path <- file.path(DATA_DIR, "ph908_sensitivity_unit_source.csv")
if (!file.exists(source_path)) {
  stop("Missing required PH908 sensitivity source table: ", source_path, call. = FALSE)
}

source_units <- readr::read_csv(source_path, show_col_types = FALSE) |>
  mutate(
    lat = as.numeric(lat),
    long = as.numeric(long),
    unit_type = as.character(unit_type),
    phylo_group = as.character(phylo_group),
    broad_family = as.character(broad_family),
    branch_class = as.character(branch_class),
    terminal_snp_norm = as.character(terminal_snp_norm),
    surname_id = as.character(surname_id),
    slava_id = as.character(slava_id),
    region_id = as.character(region_id),
    location_id = as.character(location_id)
  ) |>
  filter(unit_type %in% UNIT_TYPES)

ph908_removed_units <- source_units |>
  filter(!phylo_group %in% PH908_GROUPS)

scored_units <- bind_rows(lapply(UNIT_TYPES, function(ut) {
  ph908_removed_units |>
    filter(unit_type == ut) |>
    add_own_family_fields()
}))

results <- bind_rows(lapply(UNIT_TYPES, function(ut) {
  d <- scored_units |> filter(unit_type == ut)
  bind_rows(lapply(CONTROL_SETS, function(cs) {
    bind_rows(lapply(METRICS, function(metric_col) {
      score_case(d, ut, cs, metric_col, N_PERM)
    }))
  }))
})) |>
  arrange(unit_type, control_set, metric)

unit_counts <- scored_units |>
  mutate(is_i2 = reservoir_family == "I2") |>
  group_by(unit_type, is_i2, branch_class) |>
  summarise(n = n(), .groups = "drop") |>
  arrange(unit_type, desc(is_i2), branch_class)

readr::write_csv(results, file.path(TABLE_DIR, "ph908_removed_sensitivity.csv"))
readr::write_csv(unit_counts, file.path(TABLE_DIR, "ph908_removed_unit_counts.csv"))

message("Wrote PH908-removed sensitivity tables to: ", TABLE_DIR)
