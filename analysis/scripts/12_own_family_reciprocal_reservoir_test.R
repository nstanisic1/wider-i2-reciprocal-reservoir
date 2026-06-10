
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

set.seed(20260606L + 58L)
N_PERM <- as.integer(Sys.getenv("WIDER_I2_N_PERM", unset = "100"))
if (!is.finite(N_PERM) || N_PERM < 50L) {
  stop("WIDER_I2_N_PERM must be at least 50.", call. = FALSE)
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
FAMILIES <- c("I2", "R1a", "R1b", "E", "J2", "G", "I1")
METRICS <- c(
  "own_family_social_xsurname_xlocation",
  "own_family_near30",
  "own_family_cloud_xsurname_xlocation",
  "own_family_strict_xsurname_xlocation"
)

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
      surname_ascii = clean_primary_surname(Surname),
      slava_region = paste(Slava, Region, sep = " | ")
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
        surname_ascii != row$surname_ascii,
        Location != row$Location
      )
    if (nrow(ref) < 2 || !is.finite(row$lat) || !is.finite(row$long)) {
      return(tibble(
        own_family_social_xsurname_xlocation = NA,
        own_family_near30 = NA,
        own_family_cloud_xsurname_xlocation = NA,
        own_family_strict_xsurname_xlocation = NA,
        min_km_to_own_family_relic = NA_real_,
        own_family_ref_n = nrow(ref)
      ))
    }
    dist <- haversine_km(row$lat, row$long, ref$lat, ref$long)
    same_cell <- ref$slava_region == row$slava_region
    near <- dist <= 30
    tibble(
      own_family_social_xsurname_xlocation = any(same_cell, na.rm = TRUE),
      own_family_near30 = any(near, na.rm = TRUE),
      own_family_cloud_xsurname_xlocation = any(same_cell, na.rm = TRUE) & any(near, na.rm = TRUE),
      own_family_strict_xsurname_xlocation = any(same_cell & near, na.rm = TRUE),
      min_km_to_own_family_relic = min(dist, na.rm = TRUE),
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

  st |>
    mutate(
      unit_type = unit_type,
      control_set = control_set,
      metric = metric_col,
      n_perm = n_perm,
      target_region_branch_p = empirical_p_greater(null_target, did_lift),
      target_region_size_branch_p = empirical_p_greater(null_target_size, did_lift),
      class_target_region_p = empirical_p_greater(null_class, did_lift)
    )
}

family_observed <- function(d, unit_type, metric_col) {
  d |>
    filter(
      reservoir_family %in% FAMILIES,
      branch_class %in% c("Relic", "Founder"),
      !is.na(.data[[metric_col]])
    ) |>
    group_by(reservoir_family, branch_class) |>
    summarise(
      n = n(),
      yes = sum(.data[[metric_col]], na.rm = TRUE),
      rate = yes / n,
      .groups = "drop"
    ) |>
    select(reservoir_family, branch_class, n, yes, rate) |>
    pivot_wider(
      names_from = branch_class,
      values_from = c(n, yes, rate),
      names_sep = "_"
    ) |>
    mutate(
      relic_minus_founder = rate_Relic - rate_Founder,
      unit_type = unit_type,
      metric = metric_col
    ) |>
    arrange(desc(relic_minus_founder))
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
  cat("Scoring own-family fields for ", ut, "\n", sep = "")
  units |> filter(unit_type == ut) |> add_own_family_fields()
}))

cat("Running own-family reciprocal reservoir test with ", N_PERM, " permutations per case.\n", sep = "")
results <- bind_rows(lapply(UNIT_TYPES, function(ut) {
  d <- scored_units |> filter(unit_type == ut)
  bind_rows(lapply(CONTROL_SETS, function(cs) {
    bind_rows(lapply(METRICS, function(metric_col) {
      cat("Case ", ut, " / ", cs, " / ", metric_col, "\n", sep = "")
      score_case(d, ut, cs, metric_col, N_PERM)
    }))
  }))
})) |>
  arrange(target_region_branch_p, target_region_size_branch_p, class_target_region_p, desc(did_lift))

family_summary <- bind_rows(lapply(UNIT_TYPES, function(ut) {
  d <- scored_units |> filter(unit_type == ut)
  bind_rows(lapply(METRICS, function(metric_col) family_observed(d, ut, metric_col)))
})) |>
  group_by(unit_type, metric) |>
  mutate(
    family_rank = rank(-relic_minus_founder, ties.method = "min"),
    exact_rank_p = family_rank / n()
  ) |>
  ungroup()

ROBUST_METRICS <- intersect(
  METRICS,
  c("own_family_social_xsurname_xlocation", "own_family_cloud_xsurname_xlocation")
)
primary_metric <- "own_family_cloud_xsurname_xlocation"
robust <- bind_rows(lapply(UNIT_TYPES, function(ut) {
  d <- scored_units |> filter(unit_type == ut)
  bind_rows(lapply(CONTROL_SETS, function(cs) {
    bind_rows(lapply(ROBUST_METRICS, function(metric_col) {
      bind_rows(
        leave_one(d, ut, cs, metric_col, "Region"),
        leave_one(d, ut, cs, metric_col, "reservoir_subgroup"),
        leave_one(d, ut, cs, metric_col, "terminal_snp_norm"),
        leave_one(d, ut, cs, metric_col, "Surname"),
        leave_one(d, ut, cs, metric_col, "Slava")
      )
    }))
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

readr::write_csv(scored_units, file.path(OUT_TABLES, "own_family_reciprocal_units.csv"))
readr::write_csv(results, file.path(OUT_TABLES, "own_family_reciprocal_results.csv"))
readr::write_csv(family_summary, file.path(OUT_TABLES, "own_family_reciprocal_family_summary.csv"))
readr::write_csv(robust, file.path(OUT_TABLES, "own_family_reciprocal_leave_one.csv"))
readr::write_csv(robust_summary, file.path(OUT_TABLES, "own_family_reciprocal_robust_summary.csv"))

best_lines <- results |>
  slice_head(n = 24) |>
  mutate(line = paste0(
    "- ", unit_type, " / ", control_set, " / ", metric,
    ": DID=", sprintf("%.3f", did_lift),
    ", target p=", sprintf("%.4f", target_region_branch_p),
    ", target-size p=", sprintf("%.4f", target_region_size_branch_p),
    ", class p=", sprintf("%.4f", class_target_region_p),
    ", I2 relic/founder=", target_relic_yes, "/", target_relic_n,
    " vs ", target_founder_yes, "/", target_founder_n,
    ", control relic/founder=", control_relic_yes, "/", control_relic_n,
    " vs ", control_founder_yes, "/", control_founder_n
  )) |>
  pull(line)

family_lines <- family_summary |>
  filter(unit_type == "cluster", metric == primary_metric) |>
  arrange(family_rank) |>
  mutate(line = paste0(
    "- rank ", family_rank, ": ", reservoir_family,
    ", relic-founder=", sprintf("%.3f", relic_minus_founder),
    ", relic=", yes_Relic, "/", n_Relic,
    ", founder=", yes_Founder, "/", n_Founder,
    ", exact rank p=", sprintf("%.4f", exact_rank_p)
  )) |>
  pull(line)

robust_lines <- robust_summary |>
  filter(metric %in% ROBUST_METRICS, control_set == "non_I2_excl_R1a") |>
  mutate(line = paste0(
    "- ", unit_type, " / ", metric, " / leave-one-", factor_name,
    ": min DID=", sprintf("%.3f", min_did),
    ", worst omitted=", worst_omitted
  )) |>
  pull(line)

writeLines(
  c(
    paste0("Completed 58_own_family_reciprocal_reservoir_test.R at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("N_PERM=", N_PERM),
    paste0("Rows=", nrow(results))
  ),
  file.path(OUT_LOGS, "58_own_family_reciprocal_reservoir_test.log")
)

cat("Own-family reciprocal reservoir test complete.\n")
