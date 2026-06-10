
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
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

set.seed(20260605L + 34L)
N_PERM <- as.integer(Sys.getenv("WIDER_I2_N_PERM", unset = "500"))
if (!is.finite(N_PERM) || N_PERM < 100L) {
  stop("WIDER_I2_N_PERM must be at least 100.", call. = FALSE)
}

NEAR_KM <- 30

is_i2_unit <- function(d) {
  d$phylo_group %in% c("PH908_primary", "PH908_map_nonprimary") | d$broad_family == "nonPH908_I2"
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

empirical_p_greater <- function(null_vals, obs) {
  null_vals <- null_vals[is.finite(null_vals)]
  (1 + sum(null_vals >= obs)) / (length(null_vals) + 1)
}

ci95 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  as.numeric(stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
}

add_mutual_i2_fields <- function(d) {
  refs_all <- d |>
    filter(
      phylo_group %in% c("PH908_primary", "PH908_map_nonprimary") | broad_family == "nonPH908_I2",
      branch_class == "Relic"
    ) |>
    mutate(slava_region = paste(Slava, Region, sep = " | "))

  refs_nonph <- refs_all |>
    filter(broad_family == "nonPH908_I2")

  if (nrow(refs_all) < 10 || nrow(refs_nonph) < 8) {
    stop("Not enough I2 relic references.", call. = FALSE)
  }

  rows <- lapply(seq_len(nrow(d)), function(i) {
    row <- d[i, ]
    score_ref <- function(refs) {
      ref <- refs |> filter(phylo_group != row$phylo_group)
      if (nrow(ref) < 5) {
        return(tibble(
          social = NA,
          near = NA,
          strict = NA,
          cloud_combo = NA,
          min_km = NA_real_
        ))
      }
      cell <- paste(row$Slava, row$Region, sep = " | ")
      dist <- haversine_km(row$lat, row$long, ref$lat, ref$long)
      social <- cell %in% ref$slava_region
      near <- min(dist, na.rm = TRUE) <= NEAR_KM
      strict <- any(ref$slava_region == cell & dist <= NEAR_KM, na.rm = TRUE)
      tibble(
        social = social,
        near = near,
        strict = strict,
        cloud_combo = social & near,
        min_km = min(dist, na.rm = TRUE)
      )
    }
    all_score <- score_ref(refs_all)
    nonph_score <- score_ref(refs_nonph)
    tibble(
      other_i2_social = all_score$social,
      other_i2_near_30km = all_score$near,
      other_i2_strict_local_cell = all_score$strict,
      other_i2_cloud_social_spatial = all_score$cloud_combo,
      min_km_to_other_i2_relic = all_score$min_km,
      nonph_i2_social = nonph_score$social,
      nonph_i2_near_30km = nonph_score$near,
      nonph_i2_strict_local_cell = nonph_score$strict,
      nonph_i2_cloud_social_spatial = nonph_score$cloud_combo,
      min_km_to_nonph_i2_relic = nonph_score$min_km
    )
  })

  bind_cols(d, bind_rows(rows))
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

prepare_pair <- function(d, target_mode, control_set, metric_col) {
  target_filter <- if (target_mode == "I2_all_mutual") {
    is_i2_unit(d)
  } else if (target_mode == "PH908_primary") {
    d$phylo_group == "PH908_primary"
  } else if (target_mode == "nonPH908_I2") {
    d$broad_family == "nonPH908_I2"
  } else {
    stop("Unknown target mode.", call. = FALSE)
  }

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
      .metric_yes = .data[[metric_col]]
    ) |>
    filter(
      (.is_target_candidate | .is_control_candidate),
      branch_class %in% c("Relic", "Founder"),
      !is.na(.metric_yes)
    ) |>
    mutate(is_target = .is_target_candidate) |>
    select(-.is_target_candidate, -.is_control_candidate)
}

observed_row <- function(d, variant, unit_type, target_mode, control_set, metric_col, n_perm) {
  pair <- prepare_pair(d, target_mode, control_set, metric_col) |>
    filter(!is.na(Region), Region != "")
  if (nrow(pair) < 20) return(NULL)
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
      variant = variant,
      unit_type = unit_type,
      target_mode = target_mode,
      control_set = control_set,
      metric = metric_col,
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

leave_one_region_rows <- function(d, variant, unit_type, target_mode, control_set, metric_col) {
  pair <- prepare_pair(d, target_mode, control_set, metric_col) |>
    filter(!is.na(Region), Region != "")
  bind_rows(lapply(sort(unique(pair$Region)), function(reg) {
    x <- pair |> filter(Region != reg)
    st <- did_detail(x$is_target, x$branch_class == "Relic", as.logical(x$.metric_yes))
    if (is.null(st)) return(NULL)
    st |> mutate(
      variant = variant,
      unit_type = unit_type,
      target_mode = target_mode,
      control_set = control_set,
      metric = metric_col,
      omitted_region = reg
    )
  }))
}

units_path <- file.path(OUT_TABLES, "css_did_units.csv")
if (!file.exists(units_path)) {
  stop("Run 03_build_class_specific_units.R before this step.", call. = FALSE)
}
base_units <- readr::read_csv(units_path, show_col_types = FALSE)

variants <- c("high_resolution_upstream_excluded")
unit_types <- c("cluster", "surname_region_terminal", "location_terminal")
target_modes <- c("I2_all_mutual", "PH908_primary", "nonPH908_I2")
control_sets <- c("non_I2_excl_R1a", "all_non_I2")
metrics <- c(
  "other_i2_cloud_social_spatial",
  "other_i2_strict_local_cell",
  "nonph_i2_cloud_social_spatial"
)

scored <- bind_rows(lapply(variants, function(v) {
  bind_rows(lapply(unit_types, function(ut) {
    cat("Scoring mutual I2 fields ", v, " / ", ut, "\n", sep = "")
    base_units |> filter(variant == v, unit_type == ut) |> add_mutual_i2_fields()
  }))
}))
readr::write_csv(scored, file.path(OUT_TABLES, "mutual_i2_coretention_units.csv"))

result_grid <- bind_rows(lapply(variants, function(v) {
  bind_rows(lapply(unit_types, function(ut) {
    d <- scored |> filter(variant == v, unit_type == ut)
    bind_rows(lapply(target_modes, function(tm) {
      bind_rows(lapply(control_sets, function(cs) {
        bind_rows(lapply(metrics, function(m) {
          cat("Testing ", v, " / ", ut, " / ", tm, " / ", cs, " / ", m, "\n", sep = "")
          observed_row(d, v, ut, tm, cs, m, N_PERM)
        }))
      }))
    }))
  }))
})) |>
  select(variant, unit_type, target_mode, control_set, metric, n_perm, everything()) |>
  arrange(target_region_branch_perm_p, target_region_size_branch_perm_p, class_family_region_perm_p)

region_delete <- bind_rows(lapply(variants, function(v) {
  bind_rows(lapply(unit_types, function(ut) {
    d <- scored |> filter(variant == v, unit_type == ut)
    bind_rows(lapply(target_modes, function(tm) {
      bind_rows(lapply(control_sets, function(cs) {
        bind_rows(lapply(metrics, function(m) leave_one_region_rows(d, v, ut, tm, cs, m)))
      }))
    }))
  }))
}))

robust <- result_grid |>
  left_join(
    region_delete |>
      group_by(variant, unit_type, target_mode, control_set, metric) |>
      summarise(
        min_leave_one_region_did = min(did_lift, na.rm = TRUE),
        weakest_region = omitted_region[which.min(did_lift)],
        .groups = "drop"
      ),
    by = c("variant", "unit_type", "target_mode", "control_set", "metric")
  ) |>
  arrange(target_region_branch_perm_p, target_region_size_branch_perm_p, class_family_region_perm_p)

group_detail <- scored |>
  filter(
    phylo_group %in% c("PH908_primary", "PH908_map_nonprimary") | broad_family == "nonPH908_I2",
    branch_class %in% c("Relic", "Founder")
  ) |>
  group_by(unit_type, phylo_group, branch_class) |>
  summarise(
    n = n(),
    other_i2_cloud_yes = sum(other_i2_cloud_social_spatial, na.rm = TRUE),
    other_i2_strict_yes = sum(other_i2_strict_local_cell, na.rm = TRUE),
    nonph_cloud_yes = sum(nonph_i2_cloud_social_spatial, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(result_grid, file.path(OUT_TABLES, "mutual_i2_coretention_results.csv"))
readr::write_csv(region_delete, file.path(OUT_TABLES, "mutual_i2_coretention_leave_region.csv"))
readr::write_csv(robust, file.path(OUT_TABLES, "mutual_i2_coretention_robust.csv"))
readr::write_csv(group_detail, file.path(OUT_TABLES, "mutual_i2_coretention_group_detail.csv"))

best_lines <- robust |>
  slice_head(n = 25) |>
  mutate(line = paste0(
    "- ", unit_type, " / ", target_mode, " / ", control_set, " / ", metric,
    ": target relic=", target_relic_yes, "/", target_relic_n,
    " (", sprintf("%.3f", target_relic_rate), ")",
    ", target founder=", target_founder_yes, "/", target_founder_n,
    " (", sprintf("%.3f", target_founder_rate), ")",
    ", control relic=", control_relic_yes, "/", control_relic_n,
    " (", sprintf("%.3f", control_relic_rate), ")",
    ", control founder=", control_founder_yes, "/", control_founder_n,
    " (", sprintf("%.3f", control_founder_rate), ")",
    ", DID=", sprintf("%.3f", did_lift),
    ", target-region-branch p=", sprintf("%.4f", target_region_branch_perm_p),
    ", target-region-size-branch p=", sprintf("%.4f", target_region_size_branch_perm_p),
    ", class-family-region p=", sprintf("%.4f", class_family_region_perm_p),
    ", min leave-one-region DID=", sprintf("%.3f", min_leave_one_region_did)
  )) |>
  pull(line)

writeLines(
  c(
    paste0("Completed 04_build_mutual_i2_units.R at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("N_PERM=", N_PERM),
    paste0("Rows=", nrow(result_grid))
  ),
  file.path(OUT_LOGS, "04_build_mutual_i2_units.log")
)
