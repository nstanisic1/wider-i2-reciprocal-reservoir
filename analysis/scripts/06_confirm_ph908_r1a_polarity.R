
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
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

set.seed(20260605L + 5L)
N_PERM <- as.integer(Sys.getenv("WIDER_I2_N_PERM", unset = "1000"))
if (!is.finite(N_PERM) || N_PERM < 100L) {
  stop("WIDER_I2_N_PERM must be at least 100 for confirmation.", call. = FALSE)
}

UNIT_FILTER <- Sys.getenv("WIDER_I2_UNITS", unset = "")

normalize_snp <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^I-|^R-", "", x)
  x <- gsub("\\*$", "", x)
  x
}

make_bins <- function(x, n = 4L) {
  x <- as.numeric(x)
  if (sum(is.finite(x)) < 2) return(rep("bin_1", length(x)))
  qs <- unique(stats::quantile(x[is.finite(x)], probs = seq(0, 1, length.out = n + 1), na.rm = TRUE))
  if (length(qs) < 2) return(rep("bin_1", length(x)))
  out <- cut(x, breaks = qs, include.lowest = TRUE, labels = FALSE)
  out[is.na(out)] <- 1L
  paste0("bin_", out)
}

make_group_indices <- function(groups) {
  split(seq_along(groups), as.character(groups))
}

permute_by_indices <- function(x, idx_list) {
  out <- x
  for (idx in idx_list) {
    if (length(idx) > 1) out[idx] <- sample(out[idx], length(idx), replace = FALSE)
  }
  out
}

class_polarity_stat <- function(hg, cls, core) {
  share <- function(h, c) {
    idx <- hg == h & cls == c
    n <- sum(idx, na.rm = TRUE)
    if (n == 0) return(NA_real_)
    sum(core[idx], na.rm = TRUE) / n
  }
  ph908_relic <- share("PH908", "Relic")
  ph908_founder <- share("PH908", "Founder")
  r1a_relic <- share("R1a", "Relic")
  r1a_founder <- share("R1a", "Founder")
  (ph908_relic - ph908_founder) - (r1a_relic - r1a_founder)
}

class_polarity_detail <- function(hg, cls, core) {
  share_n <- function(h, c) {
    idx <- hg == h & cls == c
    n <- sum(idx, na.rm = TRUE)
    core_n <- sum(core[idx], na.rm = TRUE)
    tibble(n = n, core_n = core_n, share = ifelse(n > 0, core_n / n, NA_real_))
  }
  pr <- share_n("PH908", "Relic")
  pf <- share_n("PH908", "Founder")
  rr <- share_n("R1a", "Relic")
  rf <- share_n("R1a", "Founder")
  tibble(
    class_polarity_contrast = (pr$share - pf$share) - (rr$share - rf$share),
    ph908_relic_core_share = pr$share,
    ph908_founder_core_share = pf$share,
    r1a_relic_core_share = rr$share,
    r1a_founder_core_share = rf$share,
    ph908_relic_n = pr$n,
    ph908_founder_n = pf$n,
    r1a_relic_n = rr$n,
    r1a_founder_n = rf$n,
    ph908_relic_core_n = pr$core_n,
    ph908_founder_core_n = pf$core_n,
    r1a_relic_core_n = rr$core_n,
    r1a_founder_core_n = rf$core_n
  )
}

empirical_p_greater <- function(null_vals, obs) {
  null_vals <- null_vals[is.finite(null_vals)]
  (1 + sum(null_vals >= obs)) / (length(null_vals) + 1)
}

ci95 <- function(x) {
  as.numeric(stats::quantile(x[is.finite(x)], probs = c(0.025, 0.975), na.rm = TRUE))
}

summarise_model <- function(vals, obs, model) {
  q <- ci95(vals)
  tibble(
    null_model = model,
    observed = obs,
    null_mean = mean(vals, na.rm = TRUE),
    null_sd = stats::sd(vals, na.rm = TRUE),
    null_ci_low = q[1],
    null_ci_high = q[2],
    p_one_sided = empirical_p_greater(vals, obs)
  )
}

make_unit_table <- function(focus_joined, branch_tbl, unit_type) {
  if (unit_type == "terminal_branch") {
    return(branch_tbl |>
      filter(hg_group %in% c("PH908", "R1a")) |>
      transmute(
        unit_type = unit_type,
        unit_id = terminal_snp_norm,
        hg_group,
        branch_class = class,
        in_core = in_top4_core,
        source_region = dominant_region,
        terminal_snp_norm
      ) |>
      distinct())
  }
  keys <- switch(
    unit_type,
    cluster = c("cluster_id"),
    surname_slava_region_terminal = c("Surname", "Slava", "Region", "terminal_snp_norm", "hg_group", "branch_class", "in_core"),
    surname_region_terminal = c("Surname", "Region", "terminal_snp_norm", "hg_group", "branch_class", "in_core"),
    slava_region_terminal = c("Slava", "Region", "terminal_snp_norm", "hg_group", "branch_class", "in_core"),
    location_terminal = c("Location", "terminal_snp_norm", "hg_group", "branch_class", "in_core"),
    stop("Unknown unit_type: ", unit_type, call. = FALSE)
  )
  focus_joined |>
    mutate(
      unit_type = unit_type,
      source_region = Region,
      unit_id = do.call(paste, c(across(all_of(keys)), sep = " | "))
    ) |>
    select(unit_type, unit_id, hg_group, branch_class, in_core, source_region, terminal_snp_norm) |>
    distinct(unit_type, unit_id, hg_group, branch_class, in_core, terminal_snp_norm, .keep_all = TRUE)
}

run_polarity_nulls <- function(unit_tbl, unit_type, n_perm) {
  hg <- as.character(unit_tbl$hg_group)
  cls <- as.character(unit_tbl$branch_class)
  core <- as.logical(unit_tbl$in_core)
  region_bin <- as.character(unit_tbl$region_sample_bin)
  
  obs <- class_polarity_stat(hg, cls, core)
  
  idx_class_within_lineage <- make_group_indices(hg)
  idx_lineage_within_class <- make_group_indices(cls)
  idx_lineage_within_class_region <- make_group_indices(interaction(cls, region_bin, drop = TRUE))
  idx_class_within_lineage_region <- make_group_indices(interaction(hg, region_bin, drop = TRUE))
  
  null_class_within_lineage <- numeric(n_perm)
  null_lineage_within_class <- numeric(n_perm)
  null_lineage_within_class_region <- numeric(n_perm)
  null_class_within_lineage_region <- numeric(n_perm)
  
  for (i in seq_len(n_perm)) {
    cls_1 <- permute_by_indices(cls, idx_class_within_lineage)
    hg_2 <- permute_by_indices(hg, idx_lineage_within_class)
    hg_3 <- permute_by_indices(hg, idx_lineage_within_class_region)
    cls_4 <- permute_by_indices(cls, idx_class_within_lineage_region)
    
    null_class_within_lineage[i] <- class_polarity_stat(hg, cls_1, core)
    null_lineage_within_class[i] <- class_polarity_stat(hg_2, cls, core)
    null_lineage_within_class_region[i] <- class_polarity_stat(hg_3, cls, core)
    null_class_within_lineage_region[i] <- class_polarity_stat(hg, cls_4, core)
  }
  
  summary <- bind_rows(
    summarise_model(null_class_within_lineage, obs, "class_within_lineage"),
    summarise_model(null_lineage_within_class, obs, "lineage_within_class"),
    summarise_model(null_lineage_within_class_region, obs, "lineage_within_class_region_intensity"),
    summarise_model(null_class_within_lineage_region, obs, "class_within_lineage_region_intensity")
  ) |>
    mutate(unit_type = unit_type, n_perm = n_perm, n_units = nrow(unit_tbl)) |>
    select(unit_type, n_units, null_model, observed, everything())
  
  null_draws <- tibble(
    perm_id = seq_len(n_perm),
    unit_type = unit_type,
    class_within_lineage = null_class_within_lineage,
    lineage_within_class = null_lineage_within_class,
    lineage_within_class_region_intensity = null_lineage_within_class_region,
    class_within_lineage_region_intensity = null_class_within_lineage_region
  )
  
  list(summary = summary, null_draws = null_draws)
}

pre <- readRDS(file.path(PROJECT_DIR, "results", "preprocessed_dataset.rds"))
branch_tbl <- readr::read_csv(file.path(PROJECT_DIR, "tables", "table_p5_branch_table.csv"), show_col_types = FALSE)

pre <- pre |>
  mutate(
    terminal_snp_norm = normalize_snp(terminal_snp),
    hg_group = case_when(
      is_ph908_primary ~ "PH908",
      major_hg == "R1a" ~ "R1a",
      TRUE ~ major_hg
    ),
    geo_primary_ok = exclude_geo_primary == FALSE
  )

core_regions <- branch_tbl |>
  filter(hg_group == "PH908", class == "Relic") |>
  count(dominant_region, name = "n_relic_branches") |>
  arrange(desc(n_relic_branches), dominant_region) |>
  slice_head(n = 4) |>
  pull(dominant_region)

branch_tbl <- branch_tbl |>
  filter(hg_group %in% c("PH908", "R1a")) |>
  mutate(
    terminal_snp_norm = normalize_snp(terminal_snp_norm),
    in_top4_core = dominant_region %in% core_regions
  )

focus_joined <- pre |>
  filter(
    matched == TRUE,
    geo_primary_ok == TRUE,
    hg_group %in% c("PH908", "R1a")
  ) |>
  left_join(
    branch_tbl |>
      select(hg_group, terminal_snp_norm, branch_class = class, dominant_region),
    by = c("hg_group", "terminal_snp_norm")
  ) |>
  filter(!is.na(branch_class)) |>
  mutate(in_core = Region %in% core_regions)

region_intensity <- focus_joined |>
  count(Region, name = "region_focus_n") |>
  mutate(region_sample_bin = make_bins(region_focus_n, n = 4L))

unit_types <- c(
  "terminal_branch",
  "cluster",
  "surname_slava_region_terminal",
  "surname_region_terminal",
  "slava_region_terminal",
  "location_terminal"
)
if (nzchar(UNIT_FILTER)) {
  requested <- trimws(strsplit(UNIT_FILTER, ",")[[1]])
  unit_types <- intersect(unit_types, requested)
  if (!length(unit_types)) stop("WIDER_I2_UNITS did not match any known unit type.", call. = FALSE)
}

unit_tables <- lapply(unit_types, function(ut) {
  make_unit_table(focus_joined, branch_tbl, ut) |>
    left_join(region_intensity, by = c("source_region" = "Region")) |>
    mutate(
      region_focus_n = coalesce(region_focus_n, 1L),
      region_sample_bin = coalesce(region_sample_bin, "bin_1")
    )
})
names(unit_tables) <- unit_types

observed <- bind_rows(lapply(names(unit_tables), function(ut) {
  d <- unit_tables[[ut]]
  class_polarity_detail(as.character(d$hg_group), as.character(d$branch_class), as.logical(d$in_core)) |>
    mutate(unit_type = ut, n_units = nrow(d)) |>
    select(unit_type, n_units, everything())
}))
readr::write_csv(observed, file.path(OUT_TABLES, "signal_class_polarity_observed.csv"))

all_results <- lapply(names(unit_tables), function(ut) {
  cat("Polarity nulls for ", ut, " (", nrow(unit_tables[[ut]]), " units)\n", sep = "")
  run_polarity_nulls(unit_tables[[ut]], ut, N_PERM)
})
names(all_results) <- unit_types

summary_tbl <- bind_rows(lapply(all_results, `[[`, "summary")) |>
  arrange(unit_type, null_model)
draws_tbl <- bind_rows(lapply(all_results, `[[`, "null_draws"))

readr::write_csv(summary_tbl, file.path(OUT_TABLES, "signal_class_polarity_null_summary.csv"))
readr::write_csv(draws_tbl, file.path(OUT_TABLES, "signal_class_polarity_null_draws.csv"))

best_summary <- summary_tbl |>
  group_by(unit_type) |>
  summarise(
    observed = first(observed),
    max_p = max(p_one_sided, na.rm = TRUE),
    min_p = min(p_one_sided, na.rm = TRUE),
    most_conservative_null = null_model[which.max(p_one_sided)],
    n_perm = first(n_perm),
    n_units = first(n_units),
    .groups = "drop"
  ) |>
  arrange(max_p, desc(observed))
readr::write_csv(best_summary, file.path(OUT_TABLES, "signal_class_polarity_summary.csv"))

log_lines <- c(
  paste0("project_dir=", PROJECT_DIR),
  paste0("n_perm=", N_PERM),
  paste0("unit_types=", paste(unit_types, collapse = ", ")),
  paste0("summary_rows=", nrow(summary_tbl)),
  paste0("draw_rows=", nrow(draws_tbl))
)
writeLines(log_lines, file.path(OUT_LOGS, "06_confirm_ph908_r1a_polarity.log"), useBytes = TRUE)
cat(paste(log_lines, collapse = "\n"), "\n")
