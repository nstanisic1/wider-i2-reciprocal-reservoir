
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringi)
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

set.seed(20260606L + 57L)
N_PERM <- as.integer(Sys.getenv("WIDER_I2_N_PERM", unset = "500"))
if (!is.finite(N_PERM) || N_PERM < 100L) {
  stop("WIDER_I2_N_PERM must be at least 100.", call. = FALSE)
}

UNIT_TYPES <- c("cluster", "surname_region_terminal", "location_terminal")
CONTROL_SETS <- c("non_I2_excl_R1a", "all_non_I2")
I2_GROUPS <- c(
  "PH908_primary",
  "PH908_map_nonprimary",
  "Y3120_deep_mid_nonPH908",
  "Y3120_shallow_nonPH908",
  "I2_unmapped_nonPH908"
)
DEFAULT_METRICS <- c(
  "other_i2_social_xsurname",
  "other_i2_cloud_xsurname",
  "other_i2_strict_xsurname",
  "other_i2_social_xsurname_xlocation",
  "other_i2_cloud_xsurname_xlocation",
  "other_i2_strict_xsurname_xlocation"
)
metric_env <- Sys.getenv("WIDER_I2_METRICS", unset = "")
METRICS <- if (nzchar(metric_env)) {
  trimws(strsplit(metric_env, ",", fixed = TRUE)[[1]])
} else {
  DEFAULT_METRICS
}

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

to_ascii_latin <- function(x) {
  x <- ifelse(is.na(x), "", x)
  x <- stringi::stri_trans_general(x, "Any-Latin; Latin-ASCII")
  tolower(x)
}

clean_primary_surname <- function(x) {
  y <- to_ascii_latin(x)
  y <- sub("\\s*\\(.*\\)\\s*$", "", y)
  y <- gsub("\"", "", y)
  y <- gsub("[^a-z]+", "", y)
  y
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

ci95 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  as.numeric(stats::quantile(x, c(0.025, 0.975), na.rm = TRUE))
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

add_cross_surname_fields <- function(d) {
  d <- d |>
    mutate(
      surname_ascii = clean_primary_surname(Surname),
      slava_region = paste(Slava, Region, sep = " | ")
    )

  refs <- d |>
    filter(
      phylo_group %in% I2_GROUPS,
      branch_class == "Relic",
      is.finite(lat),
      is.finite(long)
    )

  rows <- lapply(seq_len(nrow(d)), function(i) {
    row <- d[i, ]
    ref <- refs |>
      filter(
        phylo_group != row$phylo_group,
        surname_ascii != row$surname_ascii
      )
    if (nrow(ref) < 5 || !is.finite(row$lat) || !is.finite(row$long)) {
      return(tibble(
        other_i2_social_xsurname = NA,
        other_i2_cloud_xsurname = NA,
        other_i2_strict_xsurname = NA,
        other_i2_social_xsurname_xlocation = NA,
        other_i2_cloud_xsurname_xlocation = NA,
        other_i2_strict_xsurname_xlocation = NA,
        min_km_to_other_i2_xsurname = NA_real_
      ))
    }
    dist <- haversine_km(row$lat, row$long, ref$lat, ref$long)
    same_cell <- ref$slava_region == row$slava_region
    diff_location <- ref$Location != row$Location
    near <- dist <= 30
    tibble(
      other_i2_social_xsurname = any(same_cell, na.rm = TRUE),
      other_i2_cloud_xsurname = any(same_cell, na.rm = TRUE) & any(near, na.rm = TRUE),
      other_i2_strict_xsurname = any(same_cell & near, na.rm = TRUE),
      other_i2_social_xsurname_xlocation = any(same_cell & diff_location, na.rm = TRUE),
      other_i2_cloud_xsurname_xlocation = any(same_cell & diff_location, na.rm = TRUE) & any(near, na.rm = TRUE),
      other_i2_strict_xsurname_xlocation = any(same_cell & diff_location & near, na.rm = TRUE),
      min_km_to_other_i2_xsurname = min(dist, na.rm = TRUE)
    )
  })

  bind_cols(d, bind_rows(rows))
}

prepare_pair <- function(d, control_set, metric_col) {
  target <- d$phylo_group %in% I2_GROUPS
  control <- control_mask(d, control_set, target)
  d |>
    mutate(
      is_target = target,
      metric_yes = .data[[metric_col]]
    ) |>
    filter(
      (is_target | control),
      branch_class %in% c("Relic", "Founder"),
      !is.na(metric_yes),
      !is.na(Region),
      Region != ""
    )
}

score_case <- function(d, unit_type, control_set, metric_col, n_perm) {
  pair <- prepare_pair(d, control_set, metric_col)
  target <- pair$is_target
  relic <- pair$branch_class == "Relic"
  yes <- as.logical(pair$metric_yes)
  st <- did_detail(target, relic, yes)
  if (is.null(st)) return(NULL)

  region_size <- make_region_size_bin(pair$Region)
  strata_target <- paste(pair$Region, pair$branch_class, sep = " | ")
  strata_target_size <- paste(pair$Region, region_size, pair$branch_class, sep = " | ")
  strata_class <- paste(target, pair$Region, sep = " | ")
  null_target <- numeric(n_perm)
  null_target_size <- numeric(n_perm)
  null_class <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    null_target[i] <- did_value(permute_by_strata(target, strata_target), relic, yes)
    null_target_size[i] <- did_value(permute_by_strata(target, strata_target_size), relic, yes)
    null_class[i] <- did_value(target, permute_by_strata(relic, strata_class), yes)
  }
  q1 <- ci95(null_target)
  q2 <- ci95(null_target_size)
  q3 <- ci95(null_class)

  st |>
    mutate(
      unit_type = unit_type,
      control_set = control_set,
      metric = metric_col,
      n_perm = n_perm,
      target_region_branch_p = empirical_p_greater(null_target, did_lift),
      target_region_size_branch_p = empirical_p_greater(null_target_size, did_lift),
      class_target_region_p = empirical_p_greater(null_class, did_lift),
      target_region_branch_ci_low = q1[1],
      target_region_branch_ci_high = q1[2],
      target_region_size_branch_ci_low = q2[1],
      target_region_size_branch_ci_high = q2[2],
      class_target_region_ci_low = q3[1],
      class_target_region_ci_high = q3[2]
    )
}

leave_one <- function(d, unit_type, control_set, metric_col, factor_name) {
  pair <- prepare_pair(d, control_set, metric_col)
  levels <- pair |>
    filter(is_target) |>
    distinct(.drop = .data[[factor_name]]) |>
    filter(!is.na(.drop), .drop != "") |>
    pull(.drop)
  bind_rows(lapply(sort(levels), function(level) {
    x <- pair |> filter(!(is_target & .data[[factor_name]] == level))
    st <- did_detail(x$is_target, x$branch_class == "Relic", as.logical(x$metric_yes))
    if (is.null(st)) return(NULL)
    st |>
      mutate(
        unit_type = unit_type,
        control_set = control_set,
        metric = metric_col,
        factor_name = factor_name,
        omitted_level = level
      )
  }))
}

units <- readr::read_csv(file.path(OUT_TABLES, "mutual_i2_coretention_units.csv"), show_col_types = FALSE) |>
  filter(variant == "high_resolution_upstream_excluded", unit_type %in% UNIT_TYPES)

scored_units <- bind_rows(lapply(UNIT_TYPES, function(ut) {
  cat("Scoring cross-surname fields for ", ut, "\n", sep = "")
  units |> filter(unit_type == ut) |> add_cross_surname_fields()
}))

cat("Running cross-surname/location reciprocal I2 tests with ", N_PERM, " permutations per case.\n", sep = "")
results <- bind_rows(lapply(UNIT_TYPES, function(ut) {
  d <- scored_units |> filter(unit_type == ut)
  bind_rows(lapply(CONTROL_SETS, function(cs) {
    bind_rows(lapply(METRICS, function(metric_col) {
      cat("Case ", ut, " / ", cs, " / ", metric_col, "\n", sep = "")
      score_case(d, ut, cs, metric_col, N_PERM)
    }))
  }))
})) |>
  arrange(target_region_branch_p, target_region_size_branch_p, class_target_region_p)

primary_metric <- "other_i2_cloud_xsurname_xlocation"
robust <- bind_rows(lapply(UNIT_TYPES, function(ut) {
  d <- scored_units |> filter(unit_type == ut)
  bind_rows(lapply(CONTROL_SETS, function(cs) {
    bind_rows(
      leave_one(d, ut, cs, primary_metric, "Region"),
      leave_one(d, ut, cs, primary_metric, "phylo_group"),
      leave_one(d, ut, cs, primary_metric, "terminal_snp_norm"),
      leave_one(d, ut, cs, primary_metric, "Surname"),
      leave_one(d, ut, cs, primary_metric, "Slava")
    )
  }))
}))

robust_summary <- robust |>
  group_by(unit_type, control_set, metric, factor_name) |>
  summarise(
    n_omissions = n(),
    min_did = min(did_lift, na.rm = TRUE),
    median_did = median(did_lift, na.rm = TRUE),
    max_did = max(did_lift, na.rm = TRUE),
    worst_omitted = omitted_level[which.min(did_lift)],
    .groups = "drop"
  )

group_detail <- scored_units |>
  filter(
    unit_type %in% UNIT_TYPES,
    phylo_group %in% I2_GROUPS,
    branch_class %in% c("Relic", "Founder")
  ) |>
  group_by(unit_type, phylo_group, branch_class) |>
  summarise(
    n = n(),
    cross_surname_social = sum(other_i2_social_xsurname, na.rm = TRUE),
    cross_surname_cloud = sum(other_i2_cloud_xsurname, na.rm = TRUE),
    cross_surname_location_cloud = sum(other_i2_cloud_xsurname_xlocation, na.rm = TRUE),
    cross_surname_location_strict = sum(other_i2_strict_xsurname_xlocation, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(scored_units, file.path(OUT_TABLES, "cross_surname_location_i2_units.csv"))
readr::write_csv(results, file.path(OUT_TABLES, "cross_surname_location_i2_results.csv"))
readr::write_csv(robust, file.path(OUT_TABLES, "cross_surname_location_i2_leave_one.csv"))
readr::write_csv(robust_summary, file.path(OUT_TABLES, "cross_surname_location_i2_robust_summary.csv"))
readr::write_csv(group_detail, file.path(OUT_TABLES, "cross_surname_location_i2_group_detail.csv"))

best_lines <- results |>
  slice_head(n = 24) |>
  mutate(line = paste0(
    "- ", unit_type, " / ", control_set, " / ", metric,
    ": DID=", sprintf("%.3f", did_lift),
    ", target p=", sprintf("%.4f", target_region_branch_p),
    ", target-size p=", sprintf("%.4f", target_region_size_branch_p),
    ", class p=", sprintf("%.4f", class_target_region_p),
    ", target relic/founder=", target_relic_yes, "/", target_relic_n,
    " vs ", target_founder_yes, "/", target_founder_n,
    ", control relic/founder=", control_relic_yes, "/", control_relic_n,
    " vs ", control_founder_yes, "/", control_founder_n
  )) |>
  pull(line)

robust_lines <- robust_summary |>
  filter(metric == primary_metric, control_set == "non_I2_excl_R1a") |>
  mutate(line = paste0(
    "- ", unit_type, " / leave-one-", factor_name,
    ": min DID=", sprintf("%.3f", min_did),
    ", worst omitted=", worst_omitted
  )) |>
  pull(line)

group_lines <- group_detail |>
  filter(unit_type == "cluster") |>
  mutate(line = paste0(
    "- ", phylo_group, " / ", branch_class,
    ": cloud x-surname-location=", cross_surname_location_cloud, "/", n,
    ", strict=", cross_surname_location_strict, "/", n
  )) |>
  pull(line)

writeLines(
  c(
    paste0("Completed 57_cross_surname_location_i2_coretention.R at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("N_PERM=", N_PERM),
    paste0("Rows=", nrow(results))
  ),
  file.path(OUT_LOGS, "57_cross_surname_location_i2_coretention.log")
)

cat("Cross-surname/location reciprocal I2 test complete.\n")
