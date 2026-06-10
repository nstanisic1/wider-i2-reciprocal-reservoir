
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

UNIT_TYPES <- c("cluster", "surname_region_terminal", "location_terminal")
CONTROL_SETS <- c("non_I2_excl_R1a", "all_non_I2")
I2_GROUPS <- c(
  "PH908_primary",
  "PH908_map_nonprimary",
  "Y3120_deep_mid_nonPH908",
  "Y3120_shallow_nonPH908",
  "I2_unmapped_nonPH908"
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

norm_snp <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("^R-", "", x)
  x <- gsub("^I-", "", x)
  x <- gsub("^E-", "", x)
  x <- gsub("^J-", "", x)
  x <- gsub("^G-", "", x)
  x <- gsub("^N-", "", x)
  x
}

coord_key <- function(lat, lon) sprintf("%.6f,%.6f", lat, lon)

make_region_size_bin <- function(region) {
  counts <- table(region)
  n <- as.numeric(counts[region])
  qs <- unique(stats::quantile(n, probs = c(0, 0.33, 0.66, 1), na.rm = TRUE))
  if (length(qs) < 2) return(rep("bin_1", length(region)))
  out <- paste0("bin_", cut(n, breaks = qs, include.lowest = TRUE, labels = FALSE))
  out[is.na(out)] <- "bin_1"
  out
}

binary_did <- function(d, target, feature) {
  relic <- d$branch_class == "Relic"
  yes <- as.logical(feature)
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

continuous_did <- function(d, target, value) {
  relic <- d$branch_class == "Relic"
  ok <- is.finite(value)
  if (sum(target & relic & ok) < 5 ||
      sum(target & !relic & ok) < 5 ||
      sum(!target & relic & ok) < 5 ||
      sum(!target & !relic & ok) < 5) return(NULL)
  tr <- mean(value[target & relic], na.rm = TRUE)
  tf <- mean(value[target & !relic], na.rm = TRUE)
  cr <- mean(value[!target & relic], na.rm = TRUE)
  cf <- mean(value[!target & !relic], na.rm = TRUE)
  tibble(
    target_relic_mean = tr,
    target_founder_mean = tf,
    control_relic_mean = cr,
    control_founder_mean = cf,
    target_relic_n = sum(target & relic & ok),
    target_founder_n = sum(target & !relic & ok),
    control_relic_n = sum(!target & relic & ok),
    control_founder_n = sum(!target & !relic & ok),
    target_relic_minus_founder = tr - tf,
    control_relic_minus_founder = cr - cf,
    did_lift = (tr - tf) - (cr - cf)
  )
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

add_relic_ecology <- function(d) {
  refs <- d |>
    filter(
      branch_class == "Relic",
      is.finite(lat),
      is.finite(long)
    ) |>
    mutate(
      is_i2_ref = phylo_group %in% I2_GROUPS,
      is_non_i2_ref = broad_family == "non_I2",
      is_r1a_ref = phylo_group == "R1a",
      is_j2_ref = phylo_group == "J2",
      is_g_ref = phylo_group == "G"
    )

  if (nrow(refs) < 50) stop("Not enough relic references.", call. = FALSE)

  rows <- lapply(seq_len(nrow(d)), function(i) {
    row <- d[i, ]
    if (!is.finite(row$lat) || !is.finite(row$long)) {
      return(tibble(
        relic_count_30 = NA_integer_,
        relic_count_50 = NA_integer_,
        i2_relic_count_20 = NA_integer_,
        i2_relic_count_30 = NA_integer_,
        i2_relic_count_50 = NA_integer_,
        i2_terminal_count_30 = NA_integer_,
        i2_group_count_30 = NA_integer_,
        non_i2_relic_count_30 = NA_integer_,
        r1a_relic_count_30 = NA_integer_,
        j2_relic_count_30 = NA_integer_,
        g_relic_count_30 = NA_integer_,
        i2_relic_share_30 = NA_real_,
        i2_relic_share_50 = NA_real_,
        min_km_to_any_relic = NA_real_,
        min_km_to_i2_relic_diff_terminal = NA_real_,
        min_km_to_non_i2_relic = NA_real_,
        min_km_to_r1a_relic = NA_real_,
        nearest_relic_is_i2 = NA,
        nearest_relic_family = NA_character_
      ))
    }
    ref <- refs |> filter(unit_id != row$unit_id)
    dist <- haversine_km(row$lat, row$long, ref$lat, ref$long)
    near30 <- dist <= 30
    near50 <- dist <= 50
    i2_diff_terminal <- ref$is_i2_ref & ref$terminal_snp_norm != row$terminal_snp_norm
    any_idx <- which(is.finite(dist))
    nearest_idx <- if (length(any_idx)) any_idx[which.min(dist[any_idx])] else NA_integer_
    nearest_family <- if (is.na(nearest_idx)) NA_character_ else ref$phylo_group[nearest_idx]
    tibble(
      relic_count_30 = sum(near30, na.rm = TRUE),
      relic_count_50 = sum(near50, na.rm = TRUE),
      i2_relic_count_20 = sum(dist <= 20 & i2_diff_terminal, na.rm = TRUE),
      i2_relic_count_30 = sum(near30 & i2_diff_terminal, na.rm = TRUE),
      i2_relic_count_50 = sum(near50 & i2_diff_terminal, na.rm = TRUE),
      i2_terminal_count_30 = dplyr::n_distinct(ref$terminal_snp_norm[near30 & i2_diff_terminal]),
      i2_group_count_30 = dplyr::n_distinct(ref$phylo_group[near30 & i2_diff_terminal]),
      non_i2_relic_count_30 = sum(near30 & ref$is_non_i2_ref, na.rm = TRUE),
      r1a_relic_count_30 = sum(near30 & ref$is_r1a_ref, na.rm = TRUE),
      j2_relic_count_30 = sum(near30 & ref$is_j2_ref, na.rm = TRUE),
      g_relic_count_30 = sum(near30 & ref$is_g_ref, na.rm = TRUE),
      i2_relic_share_30 = ifelse(sum(near30, na.rm = TRUE) > 0, sum(near30 & i2_diff_terminal, na.rm = TRUE) / sum(near30, na.rm = TRUE), NA_real_),
      i2_relic_share_50 = ifelse(sum(near50, na.rm = TRUE) > 0, sum(near50 & i2_diff_terminal, na.rm = TRUE) / sum(near50, na.rm = TRUE), NA_real_),
      min_km_to_any_relic = ifelse(length(any_idx), min(dist[any_idx], na.rm = TRUE), NA_real_),
      min_km_to_i2_relic_diff_terminal = ifelse(any(i2_diff_terminal, na.rm = TRUE), min(dist[i2_diff_terminal], na.rm = TRUE), NA_real_),
      min_km_to_non_i2_relic = ifelse(any(ref$is_non_i2_ref, na.rm = TRUE), min(dist[ref$is_non_i2_ref], na.rm = TRUE), NA_real_),
      min_km_to_r1a_relic = ifelse(any(ref$is_r1a_ref, na.rm = TRUE), min(dist[ref$is_r1a_ref], na.rm = TRUE), NA_real_),
      nearest_relic_is_i2 = ifelse(is.na(nearest_idx), NA, ref$is_i2_ref[nearest_idx]),
      nearest_relic_family = nearest_family
    )
  })

  bind_cols(d, bind_rows(rows))
}

units <- readr::read_csv(file.path(OUT_TABLES, "mutual_i2_coretention_units.csv"), show_col_types = FALSE) |>
  filter(
    variant == "high_resolution_upstream_excluded",
    unit_type %in% UNIT_TYPES,
    branch_class %in% c("Relic", "Founder")
  )

elevation_path <- file.path(OUT_TABLES, "elevation_cache_open_meteo.csv")
elevation <- if (file.exists(elevation_path)) {
  readr::read_csv(elevation_path, show_col_types = FALSE) |>
    mutate(coord_key = coord_key(lat, long)) |>
    select(coord_key, elevation_m)
} else {
  tibble(coord_key = character(), elevation_m = numeric())
}

branch_ages <- readr::read_csv(file.path(PROJECT_DIR, "branch_ages.csv"), show_col_types = FALSE) |>
  mutate(terminal_snp_norm_join = norm_snp(terminal_snp)) |>
  group_by(terminal_snp_norm_join) |>
  summarise(tmrca_ybp = median(tmrca_ybp, na.rm = TRUE), .groups = "drop")

slava_counts <- units |>
  group_by(unit_type, Slava) |>
  summarise(slava_n = n(), .groups = "drop")
slava_region_counts <- units |>
  mutate(slava_region = paste(Slava, Region, sep = " | ")) |>
  group_by(unit_type, slava_region) |>
  summarise(slava_region_n = n(), .groups = "drop")
surname_counts <- units |>
  mutate(surname_ascii_tmp = clean_primary_surname(Surname)) |>
  group_by(unit_type, surname_ascii_tmp) |>
  summarise(surname_n = n(), .groups = "drop")
surname_region_counts <- units |>
  mutate(
    surname_ascii_tmp = clean_primary_surname(Surname),
    surname_region_tmp = paste(surname_ascii_tmp, Region, sep = " | ")
  ) |>
  group_by(unit_type, surname_region_tmp) |>
  summarise(surname_region_n = n(), .groups = "drop")

features <- bind_rows(lapply(UNIT_TYPES, function(ut) {
  units |> filter(unit_type == ut) |> add_relic_ecology()
})) |>
  mutate(
    coord_key = coord_key(lat, long),
    terminal_snp_norm_join = norm_snp(terminal_snp_norm),
    surname_ascii = clean_primary_surname(Surname),
    slava_ascii = to_ascii_latin(Slava),
    slava_region = paste(Slava, Region, sep = " | "),
    surname_region_tmp = paste(surname_ascii, Region, sep = " | "),
    surname_ic = grepl("ic$", surname_ascii),
    surname_ovic_evic = grepl("(ovic|evic)$", surname_ascii),
    surname_non_ovic_ic = surname_ic & !surname_ovic_evic,
    surname_vic = grepl("vic$", surname_ascii),
    surname_ski = grepl("ski$", surname_ascii),
    slava_informative = !grepl("nepozn|nema|unknown|none|bez", slava_ascii) & nchar(slava_ascii) > 1
  ) |>
  left_join(elevation, by = "coord_key") |>
  left_join(branch_ages, by = "terminal_snp_norm_join") |>
  left_join(slava_counts, by = c("unit_type", "Slava")) |>
  left_join(slava_region_counts, by = c("unit_type", "slava_region")) |>
  left_join(surname_counts, by = c("unit_type", "surname_ascii" = "surname_ascii_tmp")) |>
  left_join(surname_region_counts, by = c("unit_type", "surname_region_tmp")) |>
  group_by(unit_type, Region) |>
  mutate(
    elevation_region_median = median(elevation_m, na.rm = TRUE),
    elevation_resid_region = elevation_m - elevation_region_median
  ) |>
  ungroup() |>
  mutate(
    other_i2_dist_le_20 = min_km_to_other_i2_relic <= 20,
    other_i2_dist_le_50 = min_km_to_other_i2_relic <= 50,
    i2_ecology_near30 = i2_relic_count_30 >= 1,
    i2_ecology_dense30 = i2_relic_count_30 >= 2,
    i2_ecology_diverse30 = i2_group_count_30 >= 2,
    i2_ecology_terminal_diverse30 = i2_terminal_count_30 >= 3,
    i2_ecology_share30_gt50 = i2_relic_share_30 > 0.50,
    i2_ecology_share50_gt50 = i2_relic_share_50 > 0.50,
    i2_ecology_no_r1a30 = i2_relic_count_30 >= 1 & r1a_relic_count_30 == 0,
    i2_ecology_i2_over_non_i2_30 = i2_relic_count_30 > non_i2_relic_count_30,
    nearest_i2_and_i2_share30_gt50 = nearest_relic_is_i2 & i2_ecology_share30_gt50,
    geography_plus_suffix_ic = other_i2_near_30km & surname_ic,
    geography_plus_ovic_evic = other_i2_near_30km & surname_ovic_evic,
    cloud_plus_suffix_ic = other_i2_cloud_social_spatial & surname_ic,
    cloud_plus_non_ovic_ic = other_i2_cloud_social_spatial & surname_non_ovic_ic,
    geo_cloud_suffix_triad = other_i2_near_30km & other_i2_cloud_social_spatial & surname_ic,
    rare_slava_region_2 = slava_region_n <= 2,
    rare_slava_region_5 = slava_region_n <= 5,
    rare_surname_region_1 = surname_region_n <= 1,
    rare_surname_region_2 = surname_region_n <= 2,
    elevation_ge_700 = elevation_m >= 700,
    elevation_ge_900 = elevation_m >= 900,
    elevation_resid_ge_100 = elevation_resid_region >= 100,
    mountain_i2_near30 = other_i2_near_30km & elevation_ge_700,
    mountain_cloud_i2 = other_i2_cloud_social_spatial & elevation_ge_700,
    old_age_known = is.finite(tmrca_ybp),
    age_ge_1500 = tmrca_ybp >= 1500,
    age_ge_2000 = tmrca_ybp >= 2000,
    i2_vs_non_i2_distance_gap = min_km_to_non_i2_relic - min_km_to_i2_relic_diff_terminal,
    i2_vs_r1a_distance_gap = min_km_to_r1a_relic - min_km_to_i2_relic_diff_terminal
  )

binary_features <- c(
  "other_i2_social",
  "other_i2_near_30km",
  "other_i2_dist_le_20",
  "other_i2_dist_le_50",
  "other_i2_cloud_social_spatial",
  "nonph_i2_social",
  "nonph_i2_near_30km",
  "nonph_i2_cloud_social_spatial",
  "surname_ic",
  "surname_ovic_evic",
  "surname_non_ovic_ic",
  "surname_vic",
  "surname_ski",
  "slava_informative",
  "rare_slava_region_2",
  "rare_slava_region_5",
  "rare_surname_region_1",
  "rare_surname_region_2",
  "i2_ecology_near30",
  "i2_ecology_dense30",
  "i2_ecology_diverse30",
  "i2_ecology_terminal_diverse30",
  "i2_ecology_share30_gt50",
  "i2_ecology_share50_gt50",
  "i2_ecology_no_r1a30",
  "i2_ecology_i2_over_non_i2_30",
  "nearest_relic_is_i2",
  "nearest_i2_and_i2_share30_gt50",
  "geography_plus_suffix_ic",
  "geography_plus_ovic_evic",
  "cloud_plus_suffix_ic",
  "cloud_plus_non_ovic_ic",
  "geo_cloud_suffix_triad",
  "elevation_ge_700",
  "elevation_ge_900",
  "elevation_resid_ge_100",
  "mountain_i2_near30",
  "mountain_cloud_i2",
  "old_age_known",
  "age_ge_1500",
  "age_ge_2000"
)

continuous_features <- c(
  "min_km_to_other_i2_relic",
  "min_km_to_nonph_i2_relic",
  "min_km_to_i2_relic_diff_terminal",
  "min_km_to_non_i2_relic",
  "min_km_to_r1a_relic",
  "i2_vs_non_i2_distance_gap",
  "i2_vs_r1a_distance_gap",
  "i2_relic_count_30",
  "i2_terminal_count_30",
  "i2_group_count_30",
  "i2_relic_share_30",
  "i2_relic_share_50",
  "elevation_m",
  "elevation_resid_region",
  "tmrca_ybp",
  "slava_n",
  "slava_region_n",
  "surname_n",
  "surname_region_n"
)

targets <- list(
  I2_all = I2_GROUPS,
  PH908_primary = "PH908_primary",
  nonPH908_I2 = c("PH908_map_nonprimary", "Y3120_deep_mid_nonPH908", "Y3120_shallow_nonPH908", "I2_unmapped_nonPH908"),
  R1a = "R1a",
  R1b = "R1b",
  J2 = "J2",
  G = "G",
  E = "E",
  I1 = "I1"
)

binary_results <- bind_rows(lapply(UNIT_TYPES, function(ut) {
  d0 <- features |> filter(unit_type == ut)
  bind_rows(lapply(names(targets), function(target_name) {
    target0 <- d0$phylo_group %in% targets[[target_name]]
    bind_rows(lapply(CONTROL_SETS, function(control_set) {
      keep <- target0 | control_mask(d0, control_set, target0)
      d <- d0[keep, ]
      target <- d$phylo_group %in% targets[[target_name]]
      bind_rows(lapply(binary_features, function(feature_name) {
        st <- binary_did(d, target, d[[feature_name]])
        if (is.null(st)) return(NULL)
        st |> mutate(
          unit_type = ut,
          target_name = target_name,
          control_set = control_set,
          feature = feature_name,
          feature_type = "binary"
        )
      }))
    }))
  }))
}))

continuous_results <- bind_rows(lapply(UNIT_TYPES, function(ut) {
  d0 <- features |> filter(unit_type == ut)
  bind_rows(lapply(names(targets), function(target_name) {
    target0 <- d0$phylo_group %in% targets[[target_name]]
    bind_rows(lapply(CONTROL_SETS, function(control_set) {
      keep <- target0 | control_mask(d0, control_set, target0)
      d <- d0[keep, ]
      target <- d$phylo_group %in% targets[[target_name]]
      bind_rows(lapply(continuous_features, function(feature_name) {
        st <- continuous_did(d, target, d[[feature_name]])
        if (is.null(st)) return(NULL)
        st |> mutate(
          unit_type = ut,
          target_name = target_name,
          control_set = control_set,
          feature = feature_name,
          feature_type = "continuous"
        )
      }))
    }))
  }))
}))

group_summary <- features |>
  group_by(unit_type, phylo_group, broad_family, branch_class) |>
  summarise(
    n = n(),
    other_i2_cloud_rate = mean(other_i2_cloud_social_spatial, na.rm = TRUE),
    other_i2_near30_rate = mean(other_i2_near_30km, na.rm = TRUE),
    i2_ecology_share30_mean = mean(i2_relic_share_30, na.rm = TRUE),
    nearest_relic_i2_rate = mean(nearest_relic_is_i2, na.rm = TRUE),
    surname_ic_rate = mean(surname_ic, na.rm = TRUE),
    elevation_median = median(elevation_m, na.rm = TRUE),
    tmrca_median = median(tmrca_ybp, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(features, file.path(OUT_TABLES, "coordinate_surname_slava_ecology_units.csv"))
readr::write_csv(binary_results, file.path(OUT_TABLES, "coordinate_surname_slava_ecology_binary_did.csv"))
readr::write_csv(continuous_results, file.path(OUT_TABLES, "coordinate_surname_slava_ecology_continuous_did.csv"))
readr::write_csv(group_summary, file.path(OUT_TABLES, "coordinate_surname_slava_ecology_group_summary.csv"))

top_i2_binary <- binary_results |>
  filter(target_name %in% c("I2_all", "PH908_primary", "nonPH908_I2")) |>
  arrange(desc(did_lift)) |>
  slice_head(n = 40)

top_i2_cont <- continuous_results |>
  filter(target_name %in% c("I2_all", "PH908_primary", "nonPH908_I2")) |>
  arrange(desc(abs(did_lift))) |>
  slice_head(n = 30)

contrast_lines <- binary_results |>
  filter(
    unit_type == "cluster",
    control_set == "non_I2_excl_R1a",
    target_name %in% c("I2_all", "R1a", "R1b", "J2", "G"),
    feature %in% c("surname_ic", "other_i2_near_30km", "other_i2_cloud_social_spatial", "i2_ecology_share30_gt50", "nearest_relic_is_i2", "geography_plus_suffix_ic")
  ) |>
  arrange(feature, desc(did_lift)) |>
  mutate(line = paste0(
    "- ", feature, " / ", target_name,
    ": DID=", sprintf("%.3f", did_lift),
    ", target relic/founder=", target_relic_yes, "/", target_relic_n,
    " vs ", target_founder_yes, "/", target_founder_n
  )) |>
  pull(line)

binary_lines <- top_i2_binary |>
  mutate(line = paste0(
    "- ", unit_type, " / ", target_name, " / ", control_set, " / ", feature,
    ": DID=", sprintf("%.3f", did_lift),
    ", target relic/founder=", target_relic_yes, "/", target_relic_n,
    " vs ", target_founder_yes, "/", target_founder_n,
    ", control relic/founder=", control_relic_yes, "/", control_relic_n,
    " vs ", control_founder_yes, "/", control_founder_n
  )) |>
  pull(line)

cont_lines <- top_i2_cont |>
  mutate(line = paste0(
    "- ", unit_type, " / ", target_name, " / ", control_set, " / ", feature,
    ": DID=", sprintf("%.3f", did_lift),
    ", target relic mean=", sprintf("%.2f", target_relic_mean),
    ", target founder mean=", sprintf("%.2f", target_founder_mean),
    ", control relic mean=", sprintf("%.2f", control_relic_mean),
    ", control founder mean=", sprintf("%.2f", control_founder_mean)
  )) |>
  pull(line)

writeLines(
  c(
    paste0("Completed 10_coordinate_surname_slava_ecology.R at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Feature rows=", nrow(features)),
    paste0("Binary result rows=", nrow(binary_results)),
    paste0("Continuous result rows=", nrow(continuous_results))
  ),
  file.path(OUT_LOGS, "10_coordinate_surname_slava_ecology.log")
)

cat("Coordinate/surname/Slava/elevation ecology scan complete.\n")
