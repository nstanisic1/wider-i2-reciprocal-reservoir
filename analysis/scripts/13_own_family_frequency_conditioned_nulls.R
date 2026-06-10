
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

set.seed(20260606L + 59L)
N_PERM <- as.integer(Sys.getenv("WIDER_I2_N_PERM", unset = "500"))
if (!is.finite(N_PERM) || N_PERM < 100L) {
  stop("WIDER_I2_N_PERM must be at least 100.", call. = FALSE)
}

parse_env_list <- function(name, default) {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) return(default)
  trimws(strsplit(val, ",", fixed = TRUE)[[1]])
}

UNIT_TYPES <- parse_env_list("WIDER_I2_UNIT_TYPES", c("cluster", "surname_region_terminal", "location_terminal"))
CONTROL_SETS <- parse_env_list("WIDER_I2_CONTROL_SETS", c("non_I2_excl_R1a", "all_non_I2"))
METRICS <- parse_env_list(
  "WIDER_I2_METRICS",
  c("own_family_social_xsurname_xlocation", "own_family_cloud_xsurname_xlocation")
)
RUN_BRANCHSIZE_NULL <- Sys.getenv("WIDER_I2_BRANCHSIZE_NULL", unset = "0") == "1"

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

bin_by_breaks <- function(x, breaks, labels) {
  out <- cut(x, breaks = breaks, labels = labels, include.lowest = TRUE, right = TRUE)
  as.character(out)
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

prepare_pair <- function(d, control_set, metric_col) {
  target <- d$reservoir_family == "I2"
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
  obs <- st$did_lift

  strata <- list(
    target_slavafreq = paste(pair$Region, pair$branch_class, pair$surname_ic, pair$slava_freq_bin, pair$ref_n_bin, sep = " | "),
    target_cellfreq = paste(pair$Region, pair$branch_class, pair$surname_ic, pair$slava_region_freq_bin, pair$ref_n_bin, sep = " | "),
    class_slavafreq = paste(target, pair$Region, pair$surname_ic, pair$slava_freq_bin, pair$ref_n_bin, sep = " | "),
    class_cellfreq = paste(target, pair$Region, pair$surname_ic, pair$slava_region_freq_bin, pair$ref_n_bin, sep = " | ")
  )
  if (RUN_BRANCHSIZE_NULL) {
    strata$target_branchsize_cellfreq <- paste(pair$Region, pair$branch_class, pair$surname_ic, pair$slava_region_freq_bin, pair$ref_n_bin, pair$branch_size_bin, sep = " | ")
    strata$class_branchsize_cellfreq <- paste(target, pair$Region, pair$surname_ic, pair$slava_region_freq_bin, pair$ref_n_bin, pair$branch_size_bin, sep = " | ")
  }

  null <- lapply(strata, function(x) numeric(n_perm))
  for (i in seq_len(n_perm)) {
    null$target_slavafreq[i] <- did_value(permute_by_strata(target, strata$target_slavafreq), relic, yes)
    null$target_cellfreq[i] <- did_value(permute_by_strata(target, strata$target_cellfreq), relic, yes)
    null$class_slavafreq[i] <- did_value(target, permute_by_strata(relic, strata$class_slavafreq), yes)
    null$class_cellfreq[i] <- did_value(target, permute_by_strata(relic, strata$class_cellfreq), yes)
    if (RUN_BRANCHSIZE_NULL) {
      null$target_branchsize_cellfreq[i] <- did_value(permute_by_strata(target, strata$target_branchsize_cellfreq), relic, yes)
      null$class_branchsize_cellfreq[i] <- did_value(target, permute_by_strata(relic, strata$class_branchsize_cellfreq), yes)
    }
  }

  st |>
    mutate(
      unit_type = unit_type,
      control_set = control_set,
      metric = metric_col,
      n_perm = n_perm,
      target_slavafreq_p = empirical_p_greater(null$target_slavafreq, obs),
      target_cellfreq_p = empirical_p_greater(null$target_cellfreq, obs),
      target_branchsize_cellfreq_p = if (RUN_BRANCHSIZE_NULL) empirical_p_greater(null$target_branchsize_cellfreq, obs) else NA_real_,
      class_slavafreq_p = empirical_p_greater(null$class_slavafreq, obs),
      class_cellfreq_p = empirical_p_greater(null$class_cellfreq, obs),
      class_branchsize_cellfreq_p = if (RUN_BRANCHSIZE_NULL) empirical_p_greater(null$class_branchsize_cellfreq, obs) else NA_real_
    )
}

units_path <- file.path(OUT_TABLES, "own_family_reciprocal_units.csv")
if (!file.exists(units_path)) {
  stop("Run 12_own_family_reciprocal_reservoir_test.R before this step.", call. = FALSE)
}

units0 <- readr::read_csv(units_path, show_col_types = FALSE) |>
  filter(unit_type %in% UNIT_TYPES)

if (!"surname_ascii" %in% names(units0)) {
  units0 <- units0 |> mutate(surname_ascii = clean_primary_surname(Surname))
}

units0 <- units0 |>
  mutate(
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
    slava_freq_bin = bin_by_breaks(slava_n, c(-Inf, 10, 30, 100, Inf), c("01_10", "11_30", "31_100", "gt100")),
    slava_region_freq_bin = bin_by_breaks(slava_region_n, c(-Inf, 1, 3, 10, Inf), c("1", "2_3", "4_10", "gt10")),
    ref_n_bin = bin_by_breaks(own_family_ref_n, c(-Inf, 5, 15, 40, 100, Inf), c("000_005", "006_015", "016_040", "041_100", "gt100")),
    branch_size_bin = bin_by_breaks(branch_size, c(-Inf, 1, 3, 10, 30, 100, Inf), c("1", "2_3", "4_10", "11_30", "31_100", "gt100"))
  )

cat("Running own-family hard-conditioned nulls with ", N_PERM, " permutations per case.\n", sep = "")
results <- bind_rows(lapply(UNIT_TYPES, function(ut) {
  d <- units |> filter(unit_type == ut)
  bind_rows(lapply(CONTROL_SETS, function(cs) {
    bind_rows(lapply(METRICS, function(metric_col) {
      cat("Case ", ut, " / ", cs, " / ", metric_col, "\n", sep = "")
      score_case(d, ut, cs, metric_col, N_PERM)
    }))
  }))
})) |>
  arrange(target_cellfreq_p, target_slavafreq_p, class_cellfreq_p, desc(did_lift))

readr::write_csv(results, file.path(OUT_TABLES, "own_family_frequency_conditioned_results.csv"))

lines <- results |>
  mutate(line = paste0(
    "- ", unit_type, " / ", control_set, " / ", metric,
    ": DID=", sprintf("%.3f", did_lift),
    ", target Slava-freq p=", sprintf("%.4f", target_slavafreq_p),
    ", target cell-freq p=", sprintf("%.4f", target_cellfreq_p),
    ifelse(RUN_BRANCHSIZE_NULL, paste0(", target branch-size+cell p=", sprintf("%.4f", target_branchsize_cellfreq_p)), ""),
    ", class Slava-freq p=", sprintf("%.4f", class_slavafreq_p),
    ", class cell-freq p=", sprintf("%.4f", class_cellfreq_p),
    ifelse(RUN_BRANCHSIZE_NULL, paste0(", class branch-size+cell p=", sprintf("%.4f", class_branchsize_cellfreq_p)), "")
  )) |>
  pull(line)

writeLines(
  c(
    paste0("Completed 59_own_family_frequency_conditioned_nulls.R at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("N_PERM=", N_PERM),
    paste0("Rows=", nrow(results))
  ),
  file.path(OUT_LOGS, "59_own_family_frequency_conditioned_nulls.log")
)

cat("Own-family frequency-conditioned nulls complete.\n")
