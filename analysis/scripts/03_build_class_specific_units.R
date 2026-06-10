
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

set.seed(20260605L + 27L)
N_PERM <- as.integer(Sys.getenv("WIDER_I2_N_PERM", unset = "500"))
if (!is.finite(N_PERM) || N_PERM < 100L) {
  stop("WIDER_I2_N_PERM must be at least 100.", call. = FALSE)
}

RELIC_MAX_N <- 3L
FOUNDER_MIN_N <- 15L
NEAR_KM <- 30

normalize_snp <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^I-|^R-", "", x)
  gsub("\\*$", "", x)
}

classify_branch <- function(n) {
  dplyr::case_when(
    n <= RELIC_MAX_N ~ "Relic",
    n >= FOUNDER_MIN_N ~ "Founder",
    TRUE ~ "Middle"
  )
}

dominant_region_fun <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
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

safe_median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
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

make_unit_table <- function(d, branch_tbl, unit_type) {
  joined <- d |>
    left_join(branch_tbl |> select(phylo_group, terminal_snp_norm, branch_size, branch_class), by = c("phylo_group", "terminal_snp_norm")) |>
    filter(branch_class %in% c("Relic", "Founder"))

  keys <- switch(
    unit_type,
    cluster = c("cluster_id"),
    surname_region_terminal = c("Surname", "Region", "terminal_snp_norm", "phylo_group", "branch_class"),
    location_terminal = c("Location", "terminal_snp_norm", "phylo_group", "branch_class"),
    slava_region_terminal = c("Slava", "Region", "terminal_snp_norm", "phylo_group", "branch_class"),
    stop("Unknown unit type.", call. = FALSE)
  )

  joined |>
    mutate(unit_type = unit_type, unit_id = do.call(paste, c(across(all_of(keys)), sep = " | "))) |>
    group_by(unit_type, unit_id, phylo_group, broad_family, branch_class, terminal_snp_norm) |>
    summarise(
      Surname = first(Surname),
      Slava = first(Slava),
      Region = first(Region),
      Location = first(Location),
      branch_size = first(branch_size),
      lat = safe_median(lat),
      long = safe_median(long),
      n_source_rows = n(),
      .groups = "drop"
    ) |>
    filter(is.finite(lat), is.finite(long))
}

add_ph908_syndrome <- function(units) {
  ph_ref <- units |>
    filter(phylo_group == "PH908_primary", branch_class == "Relic") |>
    mutate(slava_region = paste(Slava, Region, sep = " | "))

  if (nrow(ph_ref) < 8) stop("Not enough PH908 relic reference units.", call. = FALSE)

  rows <- lapply(seq_len(nrow(units)), function(i) {
    row <- units[i, ]
    dists <- haversine_km(row$lat, row$long, ph_ref$lat, ph_ref$long)
    min_km <- min(dists, na.rm = TRUE)
    tibble(
      ph908_relic_slava_region = paste(row$Slava, row$Region, sep = " | ") %in% ph_ref$slava_region,
      min_km_to_ph908_relic_cloud = min_km
    )
  })

  bind_cols(units, bind_rows(rows)) |>
    mutate(
      near_ph908_relic_cloud_30km = min_km_to_ph908_relic_cloud <= NEAR_KM,
      social_spatial_syndrome = ph908_relic_slava_region & near_ph908_relic_cloud_30km
    )
}

make_variant_units <- function(pre2, variant) {
  d <- pre2
  if (variant == "high_resolution_only") {
    d <- d |> filter(resolution_class == "high")
  } else if (variant == "high_resolution_upstream_excluded") {
    upstream <- c("PH908", "Y3120", "S17250", "M423", "L621", "CTS10228", "M172", "M170")
    d <- d |> filter(resolution_class == "high", !terminal_snp_norm %in% upstream)
  } else {
    stop("Unknown variant.", call. = FALSE)
  }

  branch_tbl <- d |>
    group_by(phylo_group, broad_family, terminal_snp_norm) |>
    summarise(branch_size = n(), dominant_region = dominant_region_fun(Region), .groups = "drop") |>
    mutate(branch_class = classify_branch(branch_size))

  bind_rows(lapply(c("cluster", "surname_region_terminal", "location_terminal", "slava_region_terminal"), function(ut) {
    make_unit_table(d, branch_tbl, ut) |> add_ph908_syndrome()
  })) |>
    mutate(variant = variant)
}

prepare_pair <- function(d, control_set) {
  target_filter <- d$broad_family == "nonPH908_I2"
  control_filter <- if (control_set == "all_non_I2") {
    d$broad_family == "non_I2"
  } else if (control_set == "non_I2_excl_R1a") {
    d$broad_family == "non_I2" & d$phylo_group != "R1a"
  } else if (control_set == "R1a_only") {
    d$phylo_group == "R1a"
  } else if (control_set == "R1a_I1") {
    d$phylo_group %in% c("R1a", "I1")
  } else {
    stop("Unknown control set.", call. = FALSE)
  }

  d |>
    filter(target_filter | control_filter) |>
    mutate(is_target = target_filter[target_filter | control_filter])
}

did_stat <- function(x) {
  shares <- x |>
    group_by(is_target, branch_class) |>
    summarise(
      n = n(),
      yes = sum(social_spatial_syndrome, na.rm = TRUE),
      rate = mean(social_spatial_syndrome, na.rm = TRUE),
      .groups = "drop"
    )
  needed <- expand.grid(is_target = c(FALSE, TRUE), branch_class = c("Relic", "Founder"), stringsAsFactors = FALSE)
  shares <- needed |> left_join(shares, by = c("is_target", "branch_class"))
  if (any(is.na(shares$n)) || any(shares$n < 5)) return(NULL)

  tr <- shares$rate[shares$is_target & shares$branch_class == "Relic"]
  tf <- shares$rate[shares$is_target & shares$branch_class == "Founder"]
  cr <- shares$rate[!shares$is_target & shares$branch_class == "Relic"]
  cf <- shares$rate[!shares$is_target & shares$branch_class == "Founder"]

  tibble(
    target_relic_rate = tr,
    target_founder_rate = tf,
    control_relic_rate = cr,
    control_founder_rate = cf,
    target_relic_n = shares$n[shares$is_target & shares$branch_class == "Relic"],
    target_founder_n = shares$n[shares$is_target & shares$branch_class == "Founder"],
    control_relic_n = shares$n[!shares$is_target & shares$branch_class == "Relic"],
    control_founder_n = shares$n[!shares$is_target & shares$branch_class == "Founder"],
    target_relic_yes = shares$yes[shares$is_target & shares$branch_class == "Relic"],
    target_founder_yes = shares$yes[shares$is_target & shares$branch_class == "Founder"],
    control_relic_yes = shares$yes[!shares$is_target & shares$branch_class == "Relic"],
    control_founder_yes = shares$yes[!shares$is_target & shares$branch_class == "Founder"],
    target_relic_minus_founder = tr - tf,
    control_relic_minus_founder = cr - cf,
    did_lift = (tr - tf) - (cr - cf)
  )
}

make_region_size_bin <- function(x) {
  region_n <- x |>
    group_by(Region) |>
    summarise(region_units = n(), .groups = "drop")
  out <- x |> left_join(region_n, by = "Region")
  qs <- unique(stats::quantile(out$region_units, probs = c(0, 0.33, 0.66, 1), na.rm = TRUE))
  if (length(qs) < 2) {
    out$region_size_bin <- "bin_1"
  } else {
    out$region_size_bin <- paste0("bin_", cut(out$region_units, breaks = qs, include.lowest = TRUE, labels = FALSE))
  }
  out$region_size_bin[is.na(out$region_size_bin)] <- "bin_1"
  out
}

permute_by_strata <- function(x, strata, n_perm) {
  null_vals <- numeric(n_perm)
  split_idx <- split(seq_along(strata), as.character(strata))
  for (i in seq_len(n_perm)) {
    y <- x
    for (idx in split_idx) {
      if (length(idx) > 1) y$is_target[idx] <- sample(y$is_target[idx], length(idx), replace = FALSE)
    }
    st <- did_stat(y)
    null_vals[i] <- if (is.null(st)) NA_real_ else st$did_lift
  }
  null_vals
}

observed_row <- function(d, variant, unit_type, control_set, n_perm) {
  pair <- prepare_pair(d, control_set) |>
    filter(branch_class %in% c("Relic", "Founder"), !is.na(Region), Region != "") |>
    make_region_size_bin()
  st <- did_stat(pair)
  if (is.null(st)) return(NULL)

  null_branch <- permute_by_strata(pair, pair$branch_class, n_perm)
  null_region_branch <- permute_by_strata(pair, paste(pair$Region, pair$branch_class, sep = " | "), n_perm)
  null_region_size_branch <- permute_by_strata(pair, paste(pair$Region, pair$region_size_bin, pair$branch_class, sep = " | "), n_perm)
  q_branch <- ci95(null_branch)
  q_region <- ci95(null_region_branch)
  q_region_size <- ci95(null_region_size_branch)

  st |>
    mutate(
      variant = variant,
      unit_type = unit_type,
      control_set = control_set,
      metric = "class_specific_social_spatial_did",
      n_perm = n_perm,
      branchclass_perm_p = empirical_p_greater(null_branch, did_lift),
      region_branchclass_perm_p = empirical_p_greater(null_region_branch, did_lift),
      region_size_branchclass_perm_p = empirical_p_greater(null_region_size_branch, did_lift),
      null_branch_mean = mean(null_branch, na.rm = TRUE),
      null_region_branch_mean = mean(null_region_branch, na.rm = TRUE),
      null_region_size_branch_mean = mean(null_region_size_branch, na.rm = TRUE),
      null_branch_ci_low = q_branch[1],
      null_branch_ci_high = q_branch[2],
      null_region_branch_ci_low = q_region[1],
      null_region_branch_ci_high = q_region[2],
      null_region_size_branch_ci_low = q_region_size[1],
      null_region_size_branch_ci_high = q_region_size[2]
    )
}

leave_one_region_rows <- function(d, variant, unit_type, control_set) {
  pair <- prepare_pair(d, control_set) |>
    filter(branch_class %in% c("Relic", "Founder"), !is.na(Region), Region != "")
  bind_rows(lapply(sort(unique(pair$Region)), function(reg) {
    x <- pair |> filter(Region != reg)
    st <- did_stat(x)
    if (is.null(st)) return(NULL)
    st |> mutate(variant = variant, unit_type = unit_type, control_set = control_set, omitted_region = reg)
  }))
}

pre <- readRDS(file.path(PROJECT_DIR, "results", "preprocessed_dataset.rds"))
branch_map <- readr::read_tsv(file.path(PROJECT_DIR, "ph908_branch_map.tsv"), show_col_types = FALSE) |>
  mutate(node_norm = normalize_snp(node))
age_tiers <- readr::read_tsv(file.path(PROJECT_DIR, "ph908_node_age_tiers.tsv"), show_col_types = FALSE) |>
  mutate(node_norm = normalize_snp(node), age_ybp = as.numeric(age_ybp))

ph908_map_nodes <- unique(branch_map$node_norm)
age_nodes <- unique(age_tiers$node_norm)

pre2 <- pre |>
  mutate(
    terminal_snp_norm = normalize_snp(terminal_snp),
    geo_primary_ok = exclude_geo_primary == FALSE,
    in_ph908_map = terminal_snp_norm %in% ph908_map_nodes,
    in_age_tiers = terminal_snp_norm %in% age_nodes
  ) |>
  left_join(age_tiers |> select(terminal_snp_norm = node_norm, age_ybp), by = "terminal_snp_norm") |>
  mutate(
    phylo_group = case_when(
      major_hg == "I2" & is_ph908_primary ~ "PH908_primary",
      major_hg == "I2" & !is_ph908_primary & in_ph908_map ~ "PH908_map_nonprimary",
      major_hg == "I2" & !is_ph908_primary & !in_ph908_map & !is.na(age_ybp) & age_ybp >= 1450 ~ "Y3120_deep_mid_nonPH908",
      major_hg == "I2" & !is_ph908_primary & !in_ph908_map & !is.na(age_ybp) & age_ybp < 1450 ~ "Y3120_shallow_nonPH908",
      major_hg == "I2" & !is_ph908_primary & !in_ph908_map & is.na(age_ybp) ~ "I2_unmapped_nonPH908",
      major_hg %in% c("R1a", "E", "J2", "I1", "R1b", "G", "N") ~ major_hg,
      TRUE ~ NA_character_
    ),
    broad_family = case_when(
      phylo_group %in% c("PH908_primary", "PH908_map_nonprimary") ~ "PH908_family",
      phylo_group %in% c("Y3120_deep_mid_nonPH908", "Y3120_shallow_nonPH908", "I2_unmapped_nonPH908") ~ "nonPH908_I2",
      !is.na(phylo_group) ~ "non_I2",
      TRUE ~ NA_character_
    )
  ) |>
  filter(
    matched == TRUE,
    geo_primary_ok == TRUE,
    !is.na(phylo_group),
    !is.na(terminal_snp_norm),
    terminal_snp_norm != "",
    !is.na(Region),
    Region != "",
    !is.na(Slava),
    Slava != "",
    is.finite(lat),
    is.finite(long)
  )

variants <- c("high_resolution_only", "high_resolution_upstream_excluded")
all_units <- bind_rows(lapply(variants, function(v) make_variant_units(pre2, v)))
readr::write_csv(all_units, file.path(OUT_TABLES, "css_did_units.csv"))

control_sets <- c("non_I2_excl_R1a", "all_non_I2", "R1a_only", "R1a_I1")

result_grid <- bind_rows(lapply(variants, function(v) {
  bind_rows(lapply(unique(all_units$unit_type), function(ut) {
    d <- all_units |> filter(variant == v, unit_type == ut)
    bind_rows(lapply(control_sets, function(cs) {
      cat("Testing ", v, " / ", ut, " / ", cs, "\n", sep = "")
      observed_row(d, v, ut, cs, N_PERM)
    }))
  }))
})) |>
  select(variant, unit_type, control_set, metric, n_perm, everything()) |>
  arrange(region_branchclass_perm_p, region_size_branchclass_perm_p, branchclass_perm_p)

region_delete <- bind_rows(lapply(variants, function(v) {
  bind_rows(lapply(unique(all_units$unit_type), function(ut) {
    d <- all_units |> filter(variant == v, unit_type == ut)
    bind_rows(lapply(control_sets, function(cs) leave_one_region_rows(d, v, ut, cs)))
  }))
}))

robust <- result_grid |>
  left_join(
    region_delete |>
      group_by(variant, unit_type, control_set) |>
      summarise(
        min_leave_one_region_did = min(did_lift, na.rm = TRUE),
        weakest_region = omitted_region[which.min(did_lift)],
        .groups = "drop"
      ),
    by = c("variant", "unit_type", "control_set")
  ) |>
  arrange(region_branchclass_perm_p, region_size_branchclass_perm_p, branchclass_perm_p)

readr::write_csv(result_grid, file.path(OUT_TABLES, "css_did_results.csv"))
readr::write_csv(region_delete, file.path(OUT_TABLES, "css_did_leave_region.csv"))
readr::write_csv(robust, file.path(OUT_TABLES, "css_did_robust.csv"))

lines <- robust |>
  slice_head(n = 24) |>
  mutate(line = paste0(
    "- ", variant, " / ", unit_type, " / ", control_set,
    ": target relic=", target_relic_yes, "/", target_relic_n,
    " (", sprintf("%.3f", target_relic_rate), ")",
    ", target founder=", target_founder_yes, "/", target_founder_n,
    " (", sprintf("%.3f", target_founder_rate), ")",
    ", control relic=", control_relic_yes, "/", control_relic_n,
    " (", sprintf("%.3f", control_relic_rate), ")",
    ", control founder=", control_founder_yes, "/", control_founder_n,
    " (", sprintf("%.3f", control_founder_rate), ")",
    ", DID=", sprintf("%.3f", did_lift),
    ", branch p=", sprintf("%.4f", branchclass_perm_p),
    ", region+branch p=", sprintf("%.4f", region_branchclass_perm_p),
    ", region-size+branch p=", sprintf("%.4f", region_size_branchclass_perm_p),
    ", min leave-one-region DID=", sprintf("%.3f", min_leave_one_region_did)
  )) |>
  pull(line)

writeLines(
  c(
    paste0("Completed 03_build_class_specific_units.R at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("N_PERM=", N_PERM),
    paste0("Rows=", nrow(result_grid))
  ),
  file.path(OUT_LOGS, "03_build_class_specific_units.log")
)
