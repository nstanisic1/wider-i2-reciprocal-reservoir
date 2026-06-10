suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
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
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
TABLES_DIR <- file.path(PROJECT_DIR, "tables")
FACTS_DIR <- file.path(PROJECT_DIR, "facts")
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FACTS_DIR, recursive = TRUE, showWarnings = FALSE)

normalize_snp <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^I-|^R-", "", x)
  gsub("\\*$", "", x)
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_first_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  x[[1]]
}

dominant_region_fun <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

classify_branch <- function(n) {
  dplyr::case_when(
    n <= 3L ~ "Relic",
    n >= 15L ~ "Founder",
    TRUE ~ "Middle"
  )
}

pre_path <- file.path(RESULTS_DIR, "preprocessed_dataset.rds")
age_path <- file.path(PROJECT_DIR, "branch_ages.csv")
if (!file.exists(pre_path)) stop("Missing preprocessed dataset: ", pre_path, call. = FALSE)
if (!file.exists(age_path)) stop("Missing branch age table: ", age_path, call. = FALSE)

cl <- readRDS(pre_path)
if (!"is_balkan_country" %in% names(cl)) {
  if (!"country_code" %in% names(cl)) stop("Missing is_balkan_country and country_code.", call. = FALSE)
  balkan_iso2 <- c("RS", "BA", "ME", "HR", "SI", "MK", "AL", "BG", "RO", "GR", "XK")
  cl <- cl |>
    mutate(is_balkan_country = toupper(.data$country_code) %in% balkan_iso2)
}
required_cols <- c(
  "matched", "exclude_geo_primary", "is_balkan_country", "terminal_snp",
  "Region", "lat", "long", "major_hg", "tmrca", "is_ph908_primary"
)
missing_cols <- setdiff(required_cols, names(cl))
if (length(missing_cols)) stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)

age_tbl <- readr::read_csv(age_path, show_col_types = FALSE) |>
  mutate(
    terminal_snp_norm = normalize_snp(.data$terminal_snp),
    tmrca_ybp = as.numeric(.data$tmrca_ybp)
  ) |>
  filter(is.finite(.data$tmrca_ybp), .data$tmrca_ybp > 0) |>
  group_by(.data$terminal_snp_norm) |>
  summarise(tmrca_ybp = safe_first_non_na(.data$tmrca_ybp), .groups = "drop")

base <- cl |>
  filter(
    .data$matched %in% TRUE,
    .data$exclude_geo_primary %in% FALSE,
    .data$is_balkan_country %in% TRUE,
    !is.na(.data$terminal_snp),
    .data$terminal_snp != "",
    !is.na(.data$Region),
    .data$Region != "",
    is.finite(.data$lat),
    is.finite(.data$long)
  ) |>
  mutate(
    terminal_snp_norm = normalize_snp(.data$terminal_snp),
    hg_group = case_when(
      .data$is_ph908_primary %in% TRUE ~ "PH908",
      .data$major_hg == "R1a" ~ "R1a",
      TRUE ~ "Other"
    )
  ) |>
  filter(.data$hg_group %in% c("PH908", "R1a"))

if (!nrow(base)) stop("No usable PH908/R1a rows after filtering.", call. = FALSE)

branch_tbl <- base |>
  group_by(.data$hg_group, .data$terminal_snp_norm) |>
  summarise(
    branch_size = n(),
    lat_c = mean(.data$lat, na.rm = TRUE),
    lon_c = mean(.data$long, na.rm = TRUE),
    dominant_region = dominant_region_fun(.data$Region),
    dataset_tmrca = safe_mean(.data$tmrca),
    .groups = "drop"
  ) |>
  left_join(age_tbl, by = "terminal_snp_norm") |>
  mutate(
    tmrca_final = coalesce(.data$tmrca_ybp, .data$dataset_tmrca),
    age_source = case_when(
      is.finite(.data$tmrca_ybp) ~ "branch_ages_csv",
      is.finite(.data$dataset_tmrca) ~ "dataset_tmrca",
      TRUE ~ "missing"
    ),
    class = classify_branch(.data$branch_size)
  ) |>
  filter(is.finite(.data$tmrca_final))

if (!nrow(branch_tbl)) stop("No branches with usable TMRCA after merging age information.", call. = FALSE)

readr::write_csv(branch_tbl, file.path(TABLES_DIR, "table_p5_branch_table.csv"))

diagnostics <- branch_tbl |>
  group_by(.data$hg_group) |>
  summarise(
    n_branches = n(),
    total_mass = sum(.data$branch_size, na.rm = TRUE),
    n_relic = sum(.data$class == "Relic", na.rm = TRUE),
    n_middle = sum(.data$class == "Middle", na.rm = TRUE),
    n_founder = sum(.data$class == "Founder", na.rm = TRUE),
    age_from_age_file = sum(.data$age_source == "branch_ages_csv", na.rm = TRUE),
    age_from_dataset = sum(.data$age_source == "dataset_tmrca", na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(diagnostics, file.path(TABLES_DIR, "table_branch_reference_diagnostics.csv"))

writeLines(
  c(
    paste0("completed_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("branches: ", nrow(branch_tbl)),
    paste0("diagnostic_rows: ", nrow(diagnostics))
  ),
  file.path(FACTS_DIR, "facts_branch_reference_table.txt"),
  useBytes = TRUE
)
