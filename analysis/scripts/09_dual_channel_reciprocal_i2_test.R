
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
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

set.seed(20260606L + 54L)
N_PERM <- as.integer(Sys.getenv("WIDER_I2_N_PERM", unset = "500"))
UNIT_TYPES <- c("cluster", "surname_region_terminal", "location_terminal")
CONTROL_SETS <- c("non_I2_excl_R1a", "all_non_I2")
I2_GROUPS <- c(
  "PH908_primary",
  "PH908_map_nonprimary",
  "Y3120_deep_mid_nonPH908",
  "Y3120_shallow_nonPH908",
  "I2_unmapped_nonPH908"
)

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

freq_bin <- function(n, breaks, labels) {
  cut(n, breaks = breaks, labels = labels, include.lowest = TRUE, right = TRUE)
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

joint_p <- function(null_geo, null_social, obs_geo, obs_social) {
  ok <- is.finite(null_geo) & is.finite(null_social)
  (1 + sum(null_geo[ok] >= obs_geo & null_social[ok] >= obs_social)) / (1 + sum(ok))
}

prepare_pair <- function(d, control_set) {
  target <- d$phylo_group %in% I2_GROUPS
  control <- if (control_set == "non_I2_excl_R1a") {
    d$broad_family == "non_I2" & d$phylo_group != "R1a"
  } else {
    d$broad_family == "non_I2"
  }
  d |>
    mutate(is_target = target) |>
    filter(
      (is_target | control),
      branch_class %in% c("Relic", "Founder"),
      !is.na(other_i2_near_30km),
      !is.na(other_i2_cloud_social_spatial),
      !is.na(Region),
      Region != ""
    )
}

score_case <- function(d, unit_type, control_set, n_perm) {
  pair <- prepare_pair(d, control_set)
  target <- pair$is_target
  relic <- pair$branch_class == "Relic"
  geo <- as.logical(pair$other_i2_near_30km)
  social <- as.logical(pair$other_i2_cloud_social_spatial)

  st_geo <- did_detail(target, relic, geo)
  st_social <- did_detail(target, relic, social)
  if (is.null(st_geo) || is.null(st_social)) return(NULL)
  obs_geo <- st_geo$did_lift
  obs_social <- st_social$did_lift

  strata <- list(
    target_slavafreq = paste(pair$Region, pair$branch_class, pair$surname_ic, pair$slava_freq_bin, sep = " | "),
    target_cellfreq = paste(pair$Region, pair$branch_class, pair$surname_ic, pair$slava_region_freq_bin, sep = " | "),
    class_slavafreq = paste(target, pair$Region, pair$surname_ic, pair$slava_freq_bin, sep = " | "),
    class_cellfreq = paste(target, pair$Region, pair$surname_ic, pair$slava_region_freq_bin, sep = " | ")
  )

  null <- lapply(strata, function(x) list(geo = numeric(n_perm), social = numeric(n_perm)))
  for (i in seq_len(n_perm)) {
    t1 <- permute_by_strata(target, strata$target_slavafreq)
    t2 <- permute_by_strata(target, strata$target_cellfreq)
    r1 <- permute_by_strata(relic, strata$class_slavafreq)
    r2 <- permute_by_strata(relic, strata$class_cellfreq)

    null$target_slavafreq$geo[i] <- did_value(t1, relic, geo)
    null$target_slavafreq$social[i] <- did_value(t1, relic, social)
    null$target_cellfreq$geo[i] <- did_value(t2, relic, geo)
    null$target_cellfreq$social[i] <- did_value(t2, relic, social)
    null$class_slavafreq$geo[i] <- did_value(target, r1, geo)
    null$class_slavafreq$social[i] <- did_value(target, r1, social)
    null$class_cellfreq$geo[i] <- did_value(target, r2, geo)
    null$class_cellfreq$social[i] <- did_value(target, r2, social)
  }

  tibble(
    unit_type = unit_type,
    control_set = control_set,
    n_perm = n_perm,
    geo_did = obs_geo,
    social_spatial_did = obs_social,
    geo_target_relic_rate = st_geo$target_relic_rate,
    geo_target_founder_rate = st_geo$target_founder_rate,
    social_target_relic_rate = st_social$target_relic_rate,
    social_target_founder_rate = st_social$target_founder_rate,
    target_relic_n = st_geo$target_relic_n,
    target_founder_n = st_geo$target_founder_n,
    target_slavafreq_joint_p = joint_p(null$target_slavafreq$geo, null$target_slavafreq$social, obs_geo, obs_social),
    target_cellfreq_joint_p = joint_p(null$target_cellfreq$geo, null$target_cellfreq$social, obs_geo, obs_social),
    class_slavafreq_joint_p = joint_p(null$class_slavafreq$geo, null$class_slavafreq$social, obs_geo, obs_social),
    class_cellfreq_joint_p = joint_p(null$class_cellfreq$geo, null$class_cellfreq$social, obs_geo, obs_social),
    target_slavafreq_geo_p = empirical_p_greater(null$target_slavafreq$geo, obs_geo),
    target_slavafreq_social_p = empirical_p_greater(null$target_slavafreq$social, obs_social),
    target_cellfreq_geo_p = empirical_p_greater(null$target_cellfreq$geo, obs_geo),
    target_cellfreq_social_p = empirical_p_greater(null$target_cellfreq$social, obs_social),
    class_slavafreq_geo_p = empirical_p_greater(null$class_slavafreq$geo, obs_geo),
    class_slavafreq_social_p = empirical_p_greater(null$class_slavafreq$social, obs_social),
    class_cellfreq_geo_p = empirical_p_greater(null$class_cellfreq$geo, obs_geo),
    class_cellfreq_social_p = empirical_p_greater(null$class_cellfreq$social, obs_social)
  )
}

units0 <- readr::read_csv(file.path(OUT_TABLES, "mutual_i2_coretention_units.csv"), show_col_types = FALSE) |>
  filter(variant == "high_resolution_upstream_excluded", unit_type %in% UNIT_TYPES) |>
  mutate(
    surname_ascii = clean_primary_surname(Surname),
    surname_ic = grepl("ic$", surname_ascii),
    slava_region = paste(Slava, Region, sep = " | ")
  )

slava_counts <- units0 |>
  group_by(unit_type, Slava) |>
  summarise(slava_n = n(), .groups = "drop")
cell_counts <- units0 |>
  group_by(unit_type, slava_region) |>
  summarise(slava_region_n = n(), .groups = "drop")

units <- units0 |>
  left_join(slava_counts, by = c("unit_type", "Slava")) |>
  left_join(cell_counts, by = c("unit_type", "slava_region")) |>
  mutate(
    slava_freq_bin = as.character(freq_bin(
      slava_n,
      breaks = c(-Inf, 10, 30, 100, Inf),
      labels = c("01_10", "11_30", "31_100", "gt100")
    )),
    slava_region_freq_bin = as.character(freq_bin(
      slava_region_n,
      breaks = c(-Inf, 1, 3, 10, Inf),
      labels = c("1", "2_3", "4_10", "gt10")
    ))
  )

cat("Running dual-channel reciprocal I2 joint test with ", N_PERM, " permutations per case.\n", sep = "")
results <- bind_rows(lapply(UNIT_TYPES, function(ut) {
  d <- units |> filter(unit_type == ut)
  bind_rows(lapply(CONTROL_SETS, function(cs) {
    cat("Case ", ut, " / ", cs, "\n", sep = "")
    score_case(d, ut, cs, N_PERM)
  }))
})) |>
  arrange(target_cellfreq_joint_p, target_slavafreq_joint_p)

readr::write_csv(results, file.path(OUT_TABLES, "dual_channel_reciprocal_i2_joint_results.csv"))

lines <- results |>
  mutate(line = paste0(
    "- ", unit_type, " / ", control_set,
    ": geo DID=", sprintf("%.3f", geo_did),
    ", social DID=", sprintf("%.3f", social_spatial_did),
    ", target Slava-freq joint p=", sprintf("%.4f", target_slavafreq_joint_p),
    ", target cell-freq joint p=", sprintf("%.4f", target_cellfreq_joint_p),
    ", class Slava-freq joint p=", sprintf("%.4f", class_slavafreq_joint_p),
    ", class cell-freq joint p=", sprintf("%.4f", class_cellfreq_joint_p)
  )) |>
  pull(line)

writeLines(
  c(
    paste0("Completed 54_dual_channel_reciprocal_i2_joint_test.R at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Rows=", nrow(results)),
    paste0("N_PERM=", N_PERM)
  ),
  file.path(OUT_LOGS, "54_dual_channel_reciprocal_i2_joint_test.log")
)

cat("Dual-channel reciprocal I2 joint test complete.\n")
