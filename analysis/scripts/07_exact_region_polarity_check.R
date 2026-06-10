
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

set.seed(20260605L + 8L)
N_PERM <- as.integer(Sys.getenv("WIDER_I2_N_PERM", unset = "500"))
if (!is.finite(N_PERM) || N_PERM < 100L) {
  stop("WIDER_I2_N_PERM must be at least 100.", call. = FALSE)
}

normalize_snp <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^I-|^R-", "", x)
  x <- gsub("\\*$", "", x)
  x
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
  (share("PH908", "Relic") - share("PH908", "Founder")) -
    (share("R1a", "Relic") - share("R1a", "Founder"))
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
  if (!length(null_vals) || !is.finite(obs)) return(NA_real_)
  (1 + sum(null_vals >= obs)) / (length(null_vals) + 1)
}

ci95 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  as.numeric(stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
}

strata_diagnostics <- function(unit_tbl) {
  lineage_exact <- unit_tbl |>
    group_by(branch_class, source_region) |>
    summarise(
      n_units = n(),
      n_ph908 = sum(hg_group == "PH908"),
      n_r1a = sum(hg_group == "R1a"),
      swappable = n_ph908 > 0 & n_r1a > 0 & n_units > 1,
      .groups = "drop"
    )
  class_exact <- unit_tbl |>
    group_by(hg_group, source_region) |>
    summarise(
      n_units = n(),
      n_relic = sum(branch_class == "Relic"),
      n_founder = sum(branch_class == "Founder"),
      n_middle = sum(branch_class == "Middle"),
      swappable = sum(c(n_relic > 0, n_founder > 0, n_middle > 0)) > 1 & n_units > 1,
      .groups = "drop"
    )
  tibble(
    lineage_exact_strata = nrow(lineage_exact),
    lineage_exact_swappable_strata = sum(lineage_exact$swappable),
    lineage_exact_units_in_swappable = sum(lineage_exact$n_units[lineage_exact$swappable]),
    lineage_exact_fraction_units_swappable = lineage_exact_units_in_swappable / nrow(unit_tbl),
    class_exact_strata = nrow(class_exact),
    class_exact_swappable_strata = sum(class_exact$swappable),
    class_exact_units_in_swappable = sum(class_exact$n_units[class_exact$swappable]),
    class_exact_fraction_units_swappable = class_exact_units_in_swappable / nrow(unit_tbl)
  )
}

run_exact_nulls <- function(unit_tbl, n_perm) {
  hg <- as.character(unit_tbl$hg_group)
  cls <- as.character(unit_tbl$branch_class)
  core <- as.logical(unit_tbl$in_core)
  region <- as.character(unit_tbl$source_region)
  obs <- class_polarity_stat(hg, cls, core)
  
  idx_lineage_exact <- make_group_indices(interaction(cls, region, drop = TRUE))
  idx_class_exact <- make_group_indices(interaction(hg, region, drop = TRUE))
  
  null_lineage_exact <- numeric(n_perm)
  null_class_exact <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    hg_1 <- permute_by_indices(hg, idx_lineage_exact)
    cls_2 <- permute_by_indices(cls, idx_class_exact)
    null_lineage_exact[i] <- class_polarity_stat(hg_1, cls, core)
    null_class_exact[i] <- class_polarity_stat(hg, cls_2, core)
  }
  summarise_one <- function(vals, model) {
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
  bind_rows(
    summarise_one(null_lineage_exact, "lineage_within_class_exact_region"),
    summarise_one(null_class_exact, "class_within_lineage_exact_region")
  )
}

make_unit_table <- function(focus_joined, branch_tbl, unit_type) {
  if (unit_type == "terminal_branch") {
    return(branch_tbl |>
      filter(hg_group %in% c("PH908", "R1a")) |>
      transmute(
        unit_type = unit_type,
        unit_id = paste(hg_group, terminal_snp_norm, sep = " | "),
        hg_group,
        branch_class = class,
        in_core = dominant_region %in% core_regions,
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

pre <- readRDS(file.path(PROJECT_DIR, "results", "preprocessed_dataset.rds"))
branch_tbl <- readr::read_csv(file.path(PROJECT_DIR, "tables", "table_p5_branch_table.csv"), show_col_types = FALSE)

pre2 <- pre |>
  mutate(
    terminal_snp_norm = normalize_snp(terminal_snp),
    hg_group = case_when(
      is_ph908_primary ~ "PH908",
      major_hg == "R1a" ~ "R1a",
      TRUE ~ NA_character_
    ),
    geo_primary_ok = exclude_geo_primary == FALSE
  ) |>
  filter(
    matched == TRUE,
    geo_primary_ok == TRUE,
    hg_group %in% c("PH908", "R1a"),
    !is.na(terminal_snp_norm),
    terminal_snp_norm != "",
    !is.na(Region),
    Region != ""
  )

branch_tbl <- branch_tbl |>
  filter(hg_group %in% c("PH908", "R1a")) |>
  mutate(terminal_snp_norm = normalize_snp(terminal_snp_norm))

core_regions <- branch_tbl |>
  filter(hg_group == "PH908", class == "Relic") |>
  count(dominant_region, name = "n_relic_branches") |>
  arrange(desc(n_relic_branches), dominant_region) |>
  slice_head(n = 4) |>
  pull(dominant_region)

focus_joined <- pre2 |>
  left_join(
    branch_tbl |>
      select(hg_group, terminal_snp_norm, branch_class = class, dominant_region),
    by = c("hg_group", "terminal_snp_norm")
  ) |>
  filter(!is.na(branch_class)) |>
  mutate(in_core = Region %in% core_regions)

unit_types <- c("terminal_branch", "cluster", "surname_slava_region_terminal", "surname_region_terminal", "slava_region_terminal", "location_terminal")

observed_all <- list()
null_all <- list()
diag_all <- list()

for (unit_type in unit_types) {
  cat("Exact-region nulls: ", unit_type, "\n", sep = "")
  unit_tbl <- make_unit_table(focus_joined, branch_tbl, unit_type)
  observed_all[[unit_type]] <- class_polarity_detail(as.character(unit_tbl$hg_group), as.character(unit_tbl$branch_class), as.logical(unit_tbl$in_core)) |>
    mutate(unit_type = unit_type, n_units = nrow(unit_tbl), n_perm = N_PERM)
  null_all[[unit_type]] <- run_exact_nulls(unit_tbl, N_PERM) |>
    mutate(unit_type = unit_type, n_units = nrow(unit_tbl), n_perm = N_PERM)
  diag_all[[unit_type]] <- strata_diagnostics(unit_tbl) |>
    mutate(unit_type = unit_type, n_units = nrow(unit_tbl))
}

observed_tbl <- bind_rows(observed_all) |>
  select(unit_type, n_units, n_perm, everything())
null_tbl <- bind_rows(null_all) |>
  select(unit_type, n_units, n_perm, null_model, observed, everything())
diag_tbl <- bind_rows(diag_all) |>
  select(unit_type, n_units, everything())
conservative_tbl <- null_tbl |>
  group_by(unit_type) |>
  summarise(
    observed = first(observed),
    max_p = max(p_one_sided, na.rm = TRUE),
    min_p = min(p_one_sided, na.rm = TRUE),
    most_conservative_null = null_model[which.max(p_one_sided)],
    n_units = first(n_units),
    n_perm = first(n_perm),
    .groups = "drop"
  ) |>
  left_join(diag_tbl, by = c("unit_type", "n_units")) |>
  arrange(max_p, desc(observed))

readr::write_csv(observed_tbl, file.path(OUT_TABLES, "exact_region_polarity_obs.csv"))
readr::write_csv(null_tbl, file.path(OUT_TABLES, "exact_region_polarity_null.csv"))
readr::write_csv(diag_tbl, file.path(OUT_TABLES, "exact_region_polarity_diag.csv"))
readr::write_csv(conservative_tbl, file.path(OUT_TABLES, "exact_region_polarity_cons.csv"))

log_lines <- c(
  paste0("project_dir=", PROJECT_DIR),
  paste0("n_perm=", N_PERM),
  paste0("unit_types=", paste(unit_types, collapse = ", ")),
  paste0("observed_rows=", nrow(observed_tbl)),
  paste0("null_rows=", nrow(null_tbl)),
  paste0("diagnostic_rows=", nrow(diag_tbl))
)
writeLines(log_lines, file.path(OUT_LOGS, "08_exact_region_polarity.log"), useBytes = TRUE)
cat(paste(log_lines, collapse = "\n"), "\n")
