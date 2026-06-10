args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[[1]] else "full"
if (!mode %in% c("full", "quick", "figures")) {
  stop("Mode must be one of: full, quick, figures.", call. = FALSE)
}

this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(this_file)) {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else file.path(getwd(), "run_all.R")
}

ROOT_DIR <- normalizePath(dirname(this_file), winslash = "/", mustWork = TRUE)
EXEC_DIR <- if (.Platform$OS.type == "windows") utils::shortPathName(ROOT_DIR) else ROOT_DIR
setwd(EXEC_DIR)

DIR_ANALYSIS_SCRIPTS <- file.path(EXEC_DIR, "analysis", "scripts")
DIR_ANALYSIS_TABLES <- file.path(EXEC_DIR, "analysis", "outputs", "tables")
DIR_ANALYSIS_LOGS <- file.path(EXEC_DIR, "analysis", "outputs", "logs")
DIR_RESULT_TABLES <- file.path(EXEC_DIR, "result_tables")
DIR_LOGS <- file.path(EXEC_DIR, "logs")
DIR_DATA <- file.path(EXEC_DIR, "data")
DIR_FIGURES <- file.path(EXEC_DIR, "figures")
DIR_FIGUREDATA <- file.path(EXEC_DIR, "figuredata")
DIR_RESULTS <- file.path(EXEC_DIR, "results")
DIR_TABLES <- file.path(EXEC_DIR, "tables")
DIR_FACTS <- file.path(EXEC_DIR, "facts")

for (d in c(DIR_ANALYSIS_TABLES, DIR_ANALYSIS_LOGS, DIR_RESULT_TABLES, DIR_LOGS, DIR_DATA, DIR_FIGURES, DIR_FIGUREDATA, DIR_RESULTS, DIR_TABLES, DIR_FACTS)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
LOG_FILE <- file.path(DIR_LOGS, paste0("run_all_", run_stamp, ".log"))
SUMMARY_FILE <- file.path(DIR_LOGS, paste0("run_all_summary_", run_stamp, ".csv"))

log_line <- function(...) {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = LOG_FILE, append = TRUE)
}

clear_directory <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(path, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  if (length(entries)) unlink(entries, recursive = TRUE, force = TRUE)
}

if (mode %in% c("full", "quick")) {
  for (d in c(DIR_ANALYSIS_TABLES, DIR_ANALYSIS_LOGS, DIR_RESULT_TABLES, DIR_FIGURES, DIR_FIGUREDATA, DIR_RESULTS, DIR_TABLES, DIR_FACTS)) {
    clear_directory(d)
  }
}

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

required_r_packages <- c("dplyr", "readr", "tidyr", "tibble", "purrr", "stringr", "stringi", "ggplot2", "scales", "patchwork")
missing_r_packages <- required_r_packages[!vapply(required_r_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_r_packages)) {
  stop("Missing required R package(s): ", paste(missing_r_packages, collapse = ", "), call. = FALSE)
}

run_cmd <- function(label, exe, args, env = character()) {
  log_line("START ", label)
  log_line("COMMAND ", exe, " ", paste(shQuote(args), collapse = " "))
  start_time <- Sys.time()
  cmd_args <- if (.Platform$OS.type == "windows") shQuote(args, type = "cmd") else args
  old_env <- list()
  if (length(env)) {
    env_split <- strsplit(env, "=", fixed = TRUE)
    env_names <- vapply(env_split, `[`, character(1), 1)
    env_vals <- vapply(env_split, function(x) paste(x[-1], collapse = "="), character(1))
    old_env <- as.list(Sys.getenv(env_names, unset = NA_character_))
    names(old_env) <- env_names
    do.call(Sys.setenv, as.list(stats::setNames(env_vals, env_names)))
    on.exit({
      for (nm in names(old_env)) {
        if (is.na(old_env[[nm]])) {
          Sys.unsetenv(nm)
        } else {
          do.call(Sys.setenv, as.list(stats::setNames(old_env[[nm]], nm)))
        }
      }
    }, add = TRUE)
  }
  out <- system2(exe, args = cmd_args, stdout = TRUE, stderr = TRUE)
  if (length(out)) cat(paste0(out, collapse = "\n"), "\n", file = LOG_FILE, append = TRUE)
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
  log_line("END ", label, " status=", status, " elapsed_sec=", elapsed)
  list(label = label, status = as.integer(status), elapsed_sec = elapsed)
}

fail_if_bad <- function(result) {
  if (!identical(result$status, 0L)) stop("Step failed: ", result$label, ". See log: ", LOG_FILE, call. = FALSE)
}

run_r <- function(label, script_name, env = character()) {
  path <- file.path(DIR_ANALYSIS_SCRIPTS, script_name)
  if (!file.exists(path)) stop("Missing script: ", path, call. = FALSE)
  res <- run_cmd(label, rscript, path, env)
  fail_if_bad(res)
  res
}

run_root_r <- function(label, path, env = character()) {
  if (!file.exists(path)) stop("Missing script: ", path, call. = FALSE)
  res <- run_cmd(label, rscript, path, env)
  fail_if_bad(res)
  res
}

copy_required <- function(from, to) {
  if (!file.exists(from)) stop("Missing expected output: ", from, call. = FALSE)
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  file.copy(from, to, overwrite = TRUE)
}

write_manifest <- function() {
  files <- list.files(".", recursive = TRUE, full.names = FALSE, all.files = TRUE, no.. = TRUE)
  rel <- gsub("\\\\", "/", files)
  drop <- grepl("^\\.git/", rel) |
    grepl("^\\.Rproj\\.user/", rel) |
    rel %in% c(".Rhistory", ".RData") |
    grepl("^analysis/outputs/", rel) |
    grepl("^results/", rel) |
    grepl("^tables/", rel) |
    grepl("^facts/", rel) |
    grepl("^figuredata/", rel) |
    grepl("^logs/", rel)
  keep <- !drop
  files <- files[keep]
  info <- file.info(files)
  manifest <- data.frame(
    path = gsub("\\\\", "/", files),
    size_bytes = info$size,
    modified = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(manifest, "pipeline_manifest.csv", row.names = FALSE, fileEncoding = "UTF-8")
}

prepare_ph908_source <- function() {
  units_path <- file.path(DIR_ANALYSIS_TABLES, "own_family_reciprocal_units.csv")
  if (!file.exists(units_path)) stop("Missing own-family unit table: ", units_path, call. = FALSE)
  units <- readr::read_csv(units_path, show_col_types = FALSE)
  required <- c("unit_type", "phylo_group", "broad_family", "branch_class", "terminal_snp_norm", "Surname", "Slava", "Region", "Location", "lat", "long")
  missing <- setdiff(required, names(units))
  if (length(missing)) stop("Missing columns for PH908 source table: ", paste(missing, collapse = ", "), call. = FALSE)
  out <- units |>
    dplyr::transmute(
      unit_type = .data$unit_type,
      phylo_group = .data$phylo_group,
      broad_family = .data$broad_family,
      branch_class = .data$branch_class,
      terminal_snp_norm = .data$terminal_snp_norm,
      surname_id = .data$Surname,
      slava_id = .data$Slava,
      region_id = .data$Region,
      location_id = .data$Location,
      lat = .data$lat,
      long = .data$long
    )
  readr::write_csv(out, file.path(DIR_DATA, "ph908_sensitivity_unit_source.csv"))
}

sync_result_tables <- function() {
  names_to_copy <- c(
    "coordinate_surname_slava_ecology_binary_did.csv",
    "coordinate_surname_slava_ecology_continuous_did.csv",
    "coordinate_surname_slava_ecology_group_summary.csv",
    "cross_surname_location_i2_results.csv",
    "cross_surname_location_i2_robust_summary.csv",
    "dual_channel_reciprocal_i2_joint_results.csv",
    "exact_region_polarity_cons.csv",
    "exact_region_polarity_null.csv",
    "exact_region_polarity_obs.csv",
    "geography_only_reciprocal_i2_results.csv",
    "reciprocal_i2_evidence_summary.csv",
    "own_family_frequency_conditioned_cloud_5000.csv",
    "own_family_frequency_conditioned_social_5000.csv",
    "own_family_reciprocal_family_summary.csv",
    "own_family_reciprocal_results.csv",
    "signal_class_polarity_summary.csv"
  )
  for (name in names_to_copy) {
    copy_required(file.path(DIR_ANALYSIS_TABLES, name), file.path(DIR_RESULT_TABLES, name))
  }
}

confirm_perm <- if (mode == "quick") "1000" else "5000"
initial_perm <- if (mode == "quick") "100" else "500"
initial_env <- c(paste0("WIDER_I2_N_PERM=", initial_perm))
confirm_env <- c(paste0("WIDER_I2_N_PERM=", confirm_perm))
summary_rows <- list()
add_step <- function(x) {
  summary_rows[[length(summary_rows) + 1L]] <<- x
}

required_inputs <- c(
  "dnk_source_snapshot.csv",
  "dnk_source_table.csv",
  "branch_ages.csv",
  "geocoded_locations_cache.rds",
  "ph908_branch_map.tsv",
  "ph908_node_age_tiers.tsv",
  "config_analysis.R"
)
missing_inputs <- required_inputs[!file.exists(file.path(EXEC_DIR, required_inputs))]
if (length(missing_inputs)) stop("Missing required input file(s): ", paste(missing_inputs, collapse = ", "), call. = FALSE)

log_line("START wider-I2 public replication pipeline")
log_line("mode=", mode)
log_line("root=", ROOT_DIR)
log_line("exec_dir=", normalizePath(getwd(), winslash = "/", mustWork = TRUE))
log_line("initial_permutations=", initial_perm)
log_line("confirmation_permutations=", confirm_perm)
log_line("rscript=", rscript)

if (mode %in% c("full", "quick")) {
  add_step(run_r("01 preprocess registry data", "01_preprocess_registry_data.R", initial_env))
  add_step(run_r("02 build branch reference table", "02_build_branch_reference_table.R", initial_env))
  add_step(run_r("03 build class-specific units", "03_build_class_specific_units.R", initial_env))
  add_step(run_r("04 build mutual I2 units", "04_build_mutual_i2_units.R", initial_env))
  add_step(run_r("05 confirm mutual I2 co-retention", "05_confirm_mutual_i2_coretention.R", confirm_env))
  add_step(run_r("06 confirm PH908/R1a polarity", "06_confirm_ph908_r1a_polarity.R", confirm_env))
  add_step(run_r("07 exact-region polarity check", "07_exact_region_polarity_check.R", initial_env))
  add_step(run_r("08 geography-only reciprocal I2 test", "08_geography_only_reciprocal_i2_test.R", initial_env))
  add_step(run_r("09 dual-channel reciprocal I2 test", "09_dual_channel_reciprocal_i2_test.R", initial_env))
  add_step(run_r("10 coordinate-surname-Slava ecology", "10_coordinate_surname_slava_ecology.R", initial_env))
  add_step(run_r("11 cross-surname/location I2 co-retention", "11_cross_surname_location_i2_coretention.R", initial_env))
  add_step(run_r("12 own-family reciprocal reservoir test", "12_own_family_reciprocal_reservoir_test.R", confirm_env))
  add_step(run_r(
    "13a own-family social hard nulls",
    "13_own_family_frequency_conditioned_nulls.R",
    c(confirm_env, "WIDER_I2_METRICS=own_family_social_xsurname_xlocation", "WIDER_I2_UNIT_TYPES=surname_region_terminal", "WIDER_I2_CONTROL_SETS=non_I2_excl_R1a")
  ))
  copy_required(
    file.path(DIR_ANALYSIS_TABLES, "own_family_frequency_conditioned_results.csv"),
    file.path(DIR_ANALYSIS_TABLES, "own_family_frequency_conditioned_social_5000.csv")
  )
  add_step(run_r(
    "13b own-family coordinate hard nulls",
    "13_own_family_frequency_conditioned_nulls.R",
    c(confirm_env, "WIDER_I2_METRICS=own_family_cloud_xsurname_xlocation", "WIDER_I2_UNIT_TYPES=surname_region_terminal", "WIDER_I2_CONTROL_SETS=non_I2_excl_R1a")
  ))
  copy_required(
    file.path(DIR_ANALYSIS_TABLES, "own_family_frequency_conditioned_results.csv"),
    file.path(DIR_ANALYSIS_TABLES, "own_family_frequency_conditioned_cloud_5000.csv")
  )
  add_step(run_r("14 summarise evidence summary", "14_summarise_evidence_summary.R", confirm_env))
  prepare_ph908_source()
  sync_result_tables()
}

add_step(run_root_r("15 render figures 1-6", file.path(EXEC_DIR, "scripts", "01_make_wider_i2_figures.R"), initial_env))
add_step(run_root_r("16 PH908-removed sensitivity", file.path(EXEC_DIR, "scripts", "02_make_ph908_removed_sensitivity.R"), c(paste0("PH908_SENSITIVITY_N_PERM=", if (mode == "quick") "500" else "5000"))))
add_step(run_root_r("17 registry-bias sensitivity", file.path(EXEC_DIR, "scripts", "09_make_registry_bias_sensitivity.R"), initial_env))

summary_df <- do.call(rbind, lapply(summary_rows, as.data.frame))
utils::write.csv(summary_df, SUMMARY_FILE, row.names = FALSE)
file.copy(SUMMARY_FILE, file.path(DIR_LOGS, "run_all_summary_latest.csv"), overwrite = TRUE)
write_manifest()

log_line("summary=", SUMMARY_FILE)
log_line("manifest=pipeline_manifest.csv")
log_line("DONE wider-I2 public replication pipeline")
