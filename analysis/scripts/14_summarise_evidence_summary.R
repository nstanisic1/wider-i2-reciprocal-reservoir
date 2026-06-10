
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

read_one <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE)
}

social_hard <- read_one(file.path(OUT_TABLES, "own_family_frequency_conditioned_social_5000.csv"))
cloud_hard <- read_one(file.path(OUT_TABLES, "own_family_frequency_conditioned_cloud_5000.csv"))
own_results <- read_one(file.path(OUT_TABLES, "own_family_reciprocal_results.csv"))
own_robust <- read_one(file.path(OUT_TABLES, "own_family_reciprocal_robust_summary.csv"))
mutual <- read_one(file.path(OUT_TABLES, "mutual_i2_coretention_confirmation_5000.csv"))
mutual_robust <- read_one(file.path(OUT_TABLES, "mutual_i2_coretention_confirmation_robust.csv"))
cross_loc <- read_one(file.path(OUT_TABLES, "cross_surname_location_i2_results.csv"))
cross_loc_robust <- read_one(file.path(OUT_TABLES, "cross_surname_location_i2_robust_summary.csv"))

robust_min <- function(tbl, unit_type, control_set, metric, factor_name) {
  out <- tbl |>
    filter(
      .data$unit_type == .env$unit_type,
      .data$control_set == .env$control_set,
      .data$metric == .env$metric,
      .data$factor_name == .env$factor_name
    ) |>
    pull(min_did)
  if (length(out) == 0) NA_real_ else out[1]
}

mutual_row <- mutual |>
  filter(unit_type == "cluster", control_set == "non_I2_excl_R1a", metric == "other_i2_cloud_social_spatial") |>
  slice(1)
mutual_robust_row <- mutual_robust |>
  filter(unit_type == "cluster", control_set == "non_I2_excl_R1a") |>
  slice(1)

cross_row <- cross_loc |>
  filter(unit_type == "cluster", control_set == "non_I2_excl_R1a", metric == "other_i2_cloud_xsurname_xlocation") |>
  slice(1)

own_social_uncond <- own_results |>
  filter(unit_type == "surname_region_terminal", control_set == "non_I2_excl_R1a", metric == "own_family_social_xsurname_xlocation") |>
  slice(1)
own_cloud_uncond <- own_results |>
  filter(unit_type == "surname_region_terminal", control_set == "non_I2_excl_R1a", metric == "own_family_cloud_xsurname_xlocation") |>
  slice(1)

summary <- bind_rows(
  tibble(
    rank = 1L,
    evidence = "Own-family reciprocal reservoir, social cross-surname/location",
    unit_type = social_hard$unit_type,
    control_set = social_hard$control_set,
    metric = social_hard$metric,
    did_lift = social_hard$did_lift,
    base_target_p = own_social_uncond$target_region_branch_p,
    base_target_size_p = own_social_uncond$target_region_size_branch_p,
    base_class_p = own_social_uncond$class_target_region_p,
    hard_target_slavafreq_p = social_hard$target_slavafreq_p,
    hard_target_cellfreq_p = social_hard$target_cellfreq_p,
    hard_class_slavafreq_p = social_hard$class_slavafreq_p,
    hard_class_cellfreq_p = social_hard$class_cellfreq_p,
    min_leave_one_region_did = robust_min(own_robust, social_hard$unit_type, social_hard$control_set, social_hard$metric, "Region"),
    min_leave_one_slava_did = robust_min(own_robust, social_hard$unit_type, social_hard$control_set, social_hard$metric, "Slava"),
    min_leave_one_surname_did = robust_min(own_robust, social_hard$unit_type, social_hard$control_set, social_hard$metric, "Surname"),
    min_leave_one_i2_subgroup_did = robust_min(own_robust, social_hard$unit_type, social_hard$control_set, social_hard$metric, "reservoir_subgroup"),
    interpretation = "Every lineage scored against its own relic field; same surname and same location excluded."
  ),
  tibble(
    rank = 2L,
    evidence = "Own-family reciprocal reservoir, coordinate-bearing cloud",
    unit_type = cloud_hard$unit_type,
    control_set = cloud_hard$control_set,
    metric = cloud_hard$metric,
    did_lift = cloud_hard$did_lift,
    base_target_p = own_cloud_uncond$target_region_branch_p,
    base_target_size_p = own_cloud_uncond$target_region_size_branch_p,
    base_class_p = own_cloud_uncond$class_target_region_p,
    hard_target_slavafreq_p = cloud_hard$target_slavafreq_p,
    hard_target_cellfreq_p = cloud_hard$target_cellfreq_p,
    hard_class_slavafreq_p = cloud_hard$class_slavafreq_p,
    hard_class_cellfreq_p = cloud_hard$class_cellfreq_p,
    min_leave_one_region_did = robust_min(own_robust, cloud_hard$unit_type, cloud_hard$control_set, cloud_hard$metric, "Region"),
    min_leave_one_slava_did = robust_min(own_robust, cloud_hard$unit_type, cloud_hard$control_set, cloud_hard$metric, "Slava"),
    min_leave_one_surname_did = robust_min(own_robust, cloud_hard$unit_type, cloud_hard$control_set, cloud_hard$metric, "Surname"),
    min_leave_one_i2_subgroup_did = robust_min(own_robust, cloud_hard$unit_type, cloud_hard$control_set, cloud_hard$metric, "reservoir_subgroup"),
    interpretation = "Same own-family design plus coordinate proximity to the own-family relic cloud."
  ),
  tibble(
    rank = 3L,
    evidence = "Cross-surname/location mutual I2 co-retention",
    unit_type = cross_row$unit_type,
    control_set = cross_row$control_set,
    metric = cross_row$metric,
    did_lift = cross_row$did_lift,
    base_target_p = cross_row$target_region_branch_p,
    base_target_size_p = cross_row$target_region_size_branch_p,
    base_class_p = cross_row$class_target_region_p,
    hard_target_slavafreq_p = NA_real_,
    hard_target_cellfreq_p = NA_real_,
    hard_class_slavafreq_p = NA_real_,
    hard_class_cellfreq_p = NA_real_,
    min_leave_one_region_did = robust_min(cross_loc_robust, cross_row$unit_type, cross_row$control_set, cross_row$metric, "Region"),
    min_leave_one_slava_did = robust_min(cross_loc_robust, cross_row$unit_type, cross_row$control_set, cross_row$metric, "Slava"),
    min_leave_one_surname_did = robust_min(cross_loc_robust, cross_row$unit_type, cross_row$control_set, cross_row$metric, "Surname"),
    min_leave_one_i2_subgroup_did = robust_min(cross_loc_robust, cross_row$unit_type, cross_row$control_set, cross_row$metric, "phylo_group"),
    interpretation = "Mutual I2 co-retention with same surname and same location excluded."
  ),
  tibble(
    rank = 4L,
    evidence = "Earlier mutual I2 co-retention baseline",
    unit_type = mutual_row$unit_type,
    control_set = mutual_row$control_set,
    metric = mutual_row$metric,
    did_lift = mutual_row$did_lift,
    base_target_p = mutual_row$target_region_branch_perm_p,
    base_target_size_p = mutual_row$target_region_size_branch_perm_p,
    base_class_p = mutual_row$class_family_region_perm_p,
    hard_target_slavafreq_p = NA_real_,
    hard_target_cellfreq_p = NA_real_,
    hard_class_slavafreq_p = NA_real_,
    hard_class_cellfreq_p = NA_real_,
    min_leave_one_region_did = mutual_robust_row$min_leave_one_region_did,
    min_leave_one_slava_did = mutual_robust_row$min_leave_one_target_factor_did_Slava,
    min_leave_one_surname_did = NA_real_,
    min_leave_one_i2_subgroup_did = mutual_robust_row$min_leave_one_target_factor_did_phylo_group,
    interpretation = "Earlier broad mutual I2 social-spatial co-retention baseline."
  )
)

readr::write_csv(summary, file.path(OUT_TABLES, "reciprocal_i2_evidence_summary.csv"))

writeLines(
  c(
    paste0("Completed 14_summarise_evidence_summary.R at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Summary rows=", nrow(summary))
  ),
  file.path(OUT_LOGS, "14_summarise_evidence_summary.log")
)

cat("Evidence summary complete.\n")
