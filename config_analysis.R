
SEED_MAIN <- 20260326L


resolve_project_dir <- function() {
  env_root <- Sys.getenv("WIDER_I2_PROJECT_DIR", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = TRUE))
  }
  
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_path <- NULL
  
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
  }
  
  candidates <- unique(na.omit(c(
    if (!is.null(script_path)) file.path(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)), "config_analysis.R") else NULL,
    file.path(getwd(), "config_analysis.R")
  )))
  
  hit <- candidates[file.exists(candidates)][1]
  
  if (!is.na(hit) && length(hit) > 0) {
    return(normalizePath(dirname(hit), winslash = "/", mustWork = TRUE))
  }
  
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

PROJECT_DIR <- resolve_project_dir()

if (!dir.exists(PROJECT_DIR)) {
  stop(
    paste0(
      "Project directory does not exist:\n",
      PROJECT_DIR,
      "\n\nSet WIDER_I2_PROJECT_DIR or place config_analysis.R in the project directory."
    ),
    call. = FALSE
  )
}

DIR_ROOT_LONG <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE)
DIR_ROOT <- if (.Platform$OS.type == "windows") utils::shortPathName(DIR_ROOT_LONG) else DIR_ROOT_LONG

DIR_INPUT <- DIR_ROOT
DIR_RESULTS <- file.path(DIR_ROOT, "results")
DIR_TABLES <- file.path(DIR_ROOT, "tables")
DIR_FIGUREDATA <- file.path(DIR_ROOT, "figuredata")
DIR_FACTS <- file.path(DIR_ROOT, "facts")
DIR_FIGURES <- file.path(DIR_ROOT, "figures")
DIR_SUPP <- file.path(DIR_ROOT, "supplement")

dir.create(DIR_RESULTS, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_TABLES, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIGUREDATA, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FACTS, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIGURES, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_SUPP, showWarnings = FALSE, recursive = TRUE)

FILE_DNK <- file.path(DIR_INPUT, "dnk_source_snapshot.csv")
FILE_BRANCH_AGES <- file.path(DIR_INPUT, "branch_ages.csv")
FILE_GEOCODE_CACHE <- file.path(DIR_INPUT, "geocoded_locations_cache.rds")

FILE_RAW_RDS <- file.path(DIR_RESULTS, "dnk_samples_raw.rds")
FILE_CLEAN_RDS <- file.path(DIR_RESULTS, "dnk_samples_clean.rds")
FILE_CLUSTERS_RDS <- file.path(DIR_RESULTS, "clusters.rds")
FILE_TMRCA_MAP_CSV <- file.path(DIR_TABLES, "tmrca_map.csv")
FILE_NODE_AGE_TIERS_CSV <- file.path(DIR_TABLES, "node_age_tiers.csv")
FILE_ANNOTATED_RDS <- file.path(DIR_RESULTS, "clusters_ph908_annotated.rds")
FILE_GEOCODED_RDS <- file.path(DIR_RESULTS, "clusters_ph908_geocoded.rds")
FILE_PREPROCESSED_RDS <- file.path(DIR_RESULTS, "preprocessed_dataset.rds")

FILE_QC_SUMMARY_CSV <- file.path(DIR_TABLES, "phase0_qc_summary.csv")
FILE_QC_EXAMPLES_CSV <- file.path(DIR_TABLES, "phase0_qc_examples.csv")
FILE_GEOCODE_QC_CSV <- file.path(DIR_TABLES, "geocode_qc_locations.csv")
FILE_GEOCODE_VALIDATION_CSV <- file.path(DIR_TABLES, "geocode_cache_validation_summary.csv")
FILE_GEOCODE_PROBLEM_CSV <- file.path(DIR_TABLES, "geocode_cache_problem_rows.csv")
FILE_GEOCODE_MISSING_CSV <- file.path(DIR_TABLES, "unexpected_locations_missing_from_cache.csv")

FILE_FIGUREDATA_PREPROCESSING_FLOW <- file.path(DIR_FIGUREDATA, "figuredata_preprocessing_flow.csv")
FILE_FIGUREDATA_PREPROCESSING_GEOCODES <- file.path(DIR_FIGUREDATA, "figuredata_preprocessing_geocodes.csv")
FILE_FIGUREDATA_PREPROCESSING_METRICS <- file.path(DIR_FIGUREDATA, "figuredata_preprocessing_metrics.csv")

FILE_GEOCODE_MAP_PNG <- file.path(DIR_SUPP, "Figure_S1B_geocode_qc_overview_map.png")
FILE_TABLE_S1_CSV <- file.path(DIR_TABLES, "Table_S1_preprocessing_metrics.csv")
FILE_FACTS_PREPROCESSING <- file.path(DIR_FACTS, "facts_preprocessing.txt")
FILE_TABLE_PREPROCESSING_SUPP <- file.path(DIR_TABLES, "table_preprocessing_supp.csv")
FILE_RESULTS_PREPROCESSING <- file.path(DIR_RESULTS, "results_preprocessing.rds")

BBOX_LAT_MIN <- 38
BBOX_LAT_MAX <- 49
BBOX_LON_MIN <- 11
BBOX_LON_MAX <- 31

BALKAN_ISO2 <- c("RS", "BA", "HR", "ME", "MK", "AL", "BG", "RO", "GR", "SI", "HU", "XK")

REQUIRED_DOWNSTREAM_COLS <- c(
  "major_hg",
  "terminal_snp",
  "is_ph908_primary",
  "Region",
  "lat",
  "long",
  "matched",
  "exclude_geo_primary"
)

COL_PH908 <- "#B22222"
COL_R1A <- "#1F77B4"

COL_RELIC <- "#8B1E3F"
COL_MIDDLE <- "#C98910"
COL_FOUNDER <- "#2E8B57"

COL_QC_ACCEPTED <- "grey45"
COL_QC_FLAGGED <- "red3"
COL_QC_MISSING <- "grey80"

stop_msg <- function(...) stop(paste0(...), call. = FALSE)

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0) b else a
}

nz <- function(x) !is.na(x) & x != ""

write_facts_header <- function(path, script_name, input_files = character()) {
  lines <- c(
    paste0("script: ", script_name),
    paste0("timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("seed: ", SEED_MAIN),
    paste0("root_dir: ", DIR_ROOT),
    if (length(input_files)) c("inputs:", paste0("  - ", input_files)) else "inputs: none",
    ""
  )
  writeLines(lines, path, useBytes = TRUE)
}

append_facts_lines <- function(path, lines) {
  cat(paste0(lines, collapse = "\n"), "\n", file = path, append = TRUE)
}
