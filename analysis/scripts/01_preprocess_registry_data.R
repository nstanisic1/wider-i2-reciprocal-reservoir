
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

locate_config_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- NULL
  idx <- grep(paste0("^", file_arg), args)
  if (length(idx) > 0) {
    script_path <- normalizePath(sub(file_arg, "", args[idx[1]]), winslash = "/", mustWork = FALSE)
  }
  
  candidates <- unique(c(
    file.path(getwd(), "config_analysis.R"),
    if (!is.null(script_path)) file.path(dirname(script_path), "config_analysis.R") else NULL,
    if (!is.null(script_path)) file.path(dirname(script_path), "..", "config_analysis.R") else NULL
  ))
  
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit) || length(hit) == 0) {
    stop(
      "Could not locate config_analysis.R. Place this script inside the project tree or run it from the project root.",
      call. = FALSE
    )
  }
  
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

CONFIG_FILE <- locate_config_file()
source(CONFIG_FILE)
setwd(DIR_ROOT)
set.seed(SEED_MAIN)

stop_msg <- function(...) stop(paste0(...), call. = FALSE)

rel_path <- function(path, root = DIR_ROOT) {
  escape_rx <- function(x) gsub("([][{}()+*^$|\\?.])", "\\\\\\1", x)
  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root_norm <- normalizePath(root, winslash = "/", mustWork = FALSE)
  root_long <- if (exists("DIR_ROOT_LONG", inherits = TRUE)) normalizePath(DIR_ROOT_LONG, winslash = "/", mustWork = FALSE) else root_norm
  roots <- unique(c(root_norm, root_long))
  for (candidate_root in roots) {
    prefix <- paste0(candidate_root, "/")
    out <- sub(paste0("^", escape_rx(prefix)), "", path_norm)
    if (!identical(out, path_norm)) return(out)
    if (identical(path_norm, candidate_root)) return(".")
  }
  basename(path_norm)
}

append_facts_lines <- function(path, lines) {
  cat(paste(lines, collapse = "\n"), "\n", file = path, append = TRUE, sep = "")
}

write_facts_header <- function(path, script_name, input_files = character()) {
  git_commit <- tryCatch(
    system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  git_commit <- if (length(git_commit)) git_commit[1] else "NA"
  
  lines <- c(
    paste0("script: ", script_name),
    paste0("timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("seed: ", SEED_MAIN),
    paste0("config_file: ", rel_path(CONFIG_FILE)),
    "project_root: .",
    paste0("git_commit: ", git_commit),
    if (length(input_files)) c("inputs:", paste0("  - ", vapply(input_files, rel_path, character(1)))) else "inputs: none",
    ""
  )
  writeLines(lines, path, useBytes = TRUE)
}

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop_msg(msg)
}

safe_n_distinct <- function(x) {
  dplyr::n_distinct(x, na.rm = TRUE)
}

norm_ws <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\\s+", " ")
  stringr::str_trim(x)
}

mojibake_rate <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(0)
  mean(stringr::str_detect(x, "Ã|Ã‘|Ãƒ"))
}

fix_mojibake <- function(x) {
  y1 <- iconv(x, from = "latin1", to = "UTF-8")
  y2 <- iconv(x, from = "CP1252", to = "UTF-8")
  ifelse(!is.na(y1), y1, y2)
}

extract_major_hg <- function(snp) {
  stringr::str_extract(dplyr::coalesce(as.character(snp), ""), "^[A-Za-z0-9]+")
}

extract_terminal_snp <- function(snp) {
  x <- dplyr::coalesce(as.character(snp), "")
  parts <- stringr::str_split(x, ">", simplify = TRUE)
  out <- parts[, ncol(parts)]
  out <- stringr::str_trim(out)
  out <- gsub("^I-", "", out)
  out <- gsub("^R-", "", out)
  out <- gsub("\\*$", "", out)
  out
}

in_bbox <- function(lat, long) {
  ifelse(
    !is.na(lat) & !is.na(long),
    lat >= BBOX_LAT_MIN & lat <= BBOX_LAT_MAX &
      long >= BBOX_LON_MIN & long <= BBOX_LON_MAX,
    FALSE
  )
}

read_ph908_branch_map <- function(path) {
  out <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE) |>
    dplyr::mutate(
      node = stringr::str_trim(as.character(node)),
      parent = dplyr::na_if(stringr::str_trim(as.character(parent)), ""),
      path_string = stringr::str_trim(as.character(path_string)),
      imputation_rule = stringr::str_trim(as.character(imputation_rule)),
      tmrca = suppressWarnings(as.numeric(tmrca))
    ) |>
    dplyr::filter(!is.na(node), node != "", !is.na(tmrca))
  
  assert_true(!anyDuplicated(out$node), "Duplicate nodes found in PH908 branch map file")
  out
}

read_node_age_tiers <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE) |>
    dplyr::mutate(
      node = stringr::str_trim(as.character(node)),
      age_ybp = suppressWarnings(as.numeric(age_ybp)),
      tier = stringr::str_trim(as.character(tier)),
      notes = stringr::str_trim(as.character(notes))
    ) |>
    dplyr::filter(!is.na(node), node != "", !is.na(age_ybp), !is.na(tier), tier != "")
}

FILE_DNK_RECON_CSV <- file.path(DIR_INPUT, "dnk_source_table.csv")
FILE_PH908_MAP_TSV <- file.path(DIR_INPUT, "ph908_branch_map.tsv")
FILE_PH908_NODE_AGE_TIERS_TSV <- file.path(DIR_INPUT, "ph908_node_age_tiers.tsv")

required_inputs <- c(
  FILE_DNK_RECON_CSV,
  FILE_PH908_MAP_TSV,
  FILE_PH908_NODE_AGE_TIERS_TSV,
  FILE_GEOCODE_CACHE,
  FILE_BRANCH_AGES
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs) > 0) {
  stop(
    paste0(
      "Missing required input file(s):\n  - ",
      paste(missing_inputs, collapse = "\n  - ")
    ),
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("01 PREPROCESSING\n")
cat("============================================================\n")
cat("Root: ", DIR_ROOT, "\n", sep = "")
cat("Seed: ", SEED_MAIN, "\n", sep = "")
cat("Inputs:\n")
cat("  - ", FILE_DNK_RECON_CSV, "\n", sep = "")
cat("  - ", FILE_PH908_MAP_TSV, "\n", sep = "")
cat("  - ", FILE_PH908_NODE_AGE_TIERS_TSV, "\n", sep = "")
cat("  - ", FILE_GEOCODE_CACHE, "\n", sep = "")
cat("  - ", FILE_BRANCH_AGES, "\n\n", sep = "")

write_facts_header(
  FILE_FACTS_PREPROCESSING,
  script_name = "01_preprocess_registry_data.R",
  input_files = c(
    FILE_DNK_RECON_CSV,
    FILE_PH908_MAP_TSV,
    FILE_PH908_NODE_AGE_TIERS_TSV,
    FILE_GEOCODE_CACHE,
    FILE_BRANCH_AGES
  )
)

cat("STEP 0: LOAD SOURCE TABLE\n")

df <- readr::read_csv(FILE_DNK_RECON_CSV, show_col_types = FALSE) |>
  dplyr::mutate(
    Id = suppressWarnings(as.integer(Id)),
    Surname = dplyr::na_if(stringr::str_trim(Surname), ""),
    Slava = dplyr::na_if(stringr::str_trim(Slava), ""),
    Region = dplyr::na_if(stringr::str_trim(Region), ""),
    Location = dplyr::na_if(stringr::str_trim(Location), ""),
    SNP = dplyr::na_if(stringr::str_trim(SNP), "")
  )

assert_true(nrow(df) >= 1, "Source table is empty")

raw_summary <- list(
  n_rows = nrow(df),
  id_na = sum(is.na(df$Id)),
  missing_location = sum(is.na(df$Location)),
  missing_snp = sum(is.na(df$SNP)),
  unique_surnames = safe_n_distinct(df$Surname),
  unique_locations = safe_n_distinct(df$Location),
  snp_has_arrow_pct = round(100 * mean(stringr::str_detect(dplyr::coalesce(df$SNP, ""), ">")), 2)
)

append_facts_lines(FILE_FACTS_PREPROCESSING, c(
  "[raw_source_summary]",
  paste0("n_rows: ", raw_summary$n_rows),
  paste0("id_na: ", raw_summary$id_na),
  paste0("missing_location: ", raw_summary$missing_location),
  paste0("missing_snp: ", raw_summary$missing_snp),
  paste0("unique_surnames: ", raw_summary$unique_surnames),
  paste0("unique_locations: ", raw_summary$unique_locations),
  paste0("snp_has_arrow_pct: ", raw_summary$snp_has_arrow_pct),
  ""
))

cat("  Raw rows: ", raw_summary$n_rows, "\n\n", sep = "")

cat("STEP 1: ENCODING QC\n")

df_clean <- df

fields <- c("Surname", "Slava", "Region", "Location", "SNP")
rates <- sapply(fields, function(f) mojibake_rate(df_clean[[f]]))
needs_fix <- any(rates[c("Surname", "Slava", "Region", "Location")] > 0.05)

if (needs_fix) {
  for (f in c("Surname", "Slava", "Region", "Location")) {
    df_clean[[f]] <- fix_mojibake(df_clean[[f]])
  }
}

encoding_summary <- list(
  mojibake_rate_surname = unname(rates["Surname"]),
  mojibake_rate_slava = unname(rates["Slava"]),
  mojibake_rate_region = unname(rates["Region"]),
  mojibake_rate_location = unname(rates["Location"]),
  mojibake_rate_snp = unname(rates["SNP"]),
  fixes_applied = needs_fix
)

append_facts_lines(FILE_FACTS_PREPROCESSING, c(
  "[encoding_qc_summary]",
  paste0("mojibake_rate_surname: ", round(encoding_summary$mojibake_rate_surname, 6)),
  paste0("mojibake_rate_slava: ", round(encoding_summary$mojibake_rate_slava, 6)),
  paste0("mojibake_rate_region: ", round(encoding_summary$mojibake_rate_region, 6)),
  paste0("mojibake_rate_location: ", round(encoding_summary$mojibake_rate_location, 6)),
  paste0("mojibake_rate_snp: ", round(encoding_summary$mojibake_rate_snp, 6)),
  paste0("fixes_applied: ", encoding_summary$fixes_applied),
  ""
))

cat("STEP 2: BUILD CLUSTERS\n")

df2 <- df_clean |>
  dplyr::mutate(
    Surname_n = norm_ws(Surname),
    Slava_n = norm_ws(Slava),
    Location_n = norm_ws(Location),
    SNP_n = norm_ws(SNP)
  ) |>
  dplyr::mutate(
    cluster_id = stringr::str_c(
      dplyr::coalesce(Surname_n, ""),
      dplyr::coalesce(Slava_n, ""),
      dplyr::coalesce(Location_n, ""),
      sep = " | "
    )
  )

clusters <- df2 |>
  dplyr::group_by(cluster_id) |>
  dplyr::summarise(
    n_members = dplyr::n(),
    Surname = dplyr::first(Surname_n),
    Slava = dplyr::first(Slava_n),
    Region = dplyr::first(Region),
    Location = dplyr::first(Location_n),
    SNP_any = names(sort(table(SNP_n), decreasing = TRUE))[1],
    .groups = "drop"
  )

assert_true(!anyDuplicated(clusters$cluster_id), "Duplicated cluster_id found after clustering")
assert_true(nrow(clusters) > 0, "Cluster table is empty")

cluster_summary <- list(
  n_clusters = nrow(clusters),
  singleton_clusters = sum(clusters$n_members == 1, na.rm = TRUE),
  max_cluster_size = max(clusters$n_members, na.rm = TRUE),
  mean_cluster_size = round(mean(clusters$n_members, na.rm = TRUE), 4),
  median_cluster_size = stats::median(clusters$n_members, na.rm = TRUE),
  missing_snp_any = sum(is.na(clusters$SNP_any))
)

append_facts_lines(FILE_FACTS_PREPROCESSING, c(
  "[cluster_summary]",
  paste0("n_clusters: ", cluster_summary$n_clusters),
  paste0("singleton_clusters: ", cluster_summary$singleton_clusters),
  paste0("max_cluster_size: ", cluster_summary$max_cluster_size),
  paste0("mean_cluster_size: ", cluster_summary$mean_cluster_size),
  paste0("median_cluster_size: ", cluster_summary$median_cluster_size),
  paste0("missing_snp_any: ", cluster_summary$missing_snp_any),
  ""
))

cat("  Clusters: ", nrow(clusters), "\n\n", sep = "")

cat("STEP 3: ANNOTATE PH908 STRUCTURE\n")

map_ph908_only <- read_ph908_branch_map(FILE_PH908_MAP_TSV)
age_df <- read_node_age_tiers(FILE_PH908_NODE_AGE_TIERS_TSV)

required_anchor_nodes <- c("PH908", "FT14506", "FT16449", "Z16983")
assert_true(all(required_anchor_nodes %in% map_ph908_only$node), "Missing required PH908 anchor nodes")

major_hg <- extract_major_hg(clusters$SNP_any)
terminal_snp <- extract_terminal_snp(clusters$SNP_any)

is_i2 <- major_hg == "I2"
node_set <- unique(map_ph908_only$node)
is_ph908_wide <- is_i2 & (terminal_snp %in% node_set | terminal_snp == "PH908")

parent_lookup <- stats::setNames(map_ph908_only$parent, map_ph908_only$node)

infer_branch <- function(node) {
  if (is.na(node) || node == "") return(NA_character_)
  n <- node
  for (i in 1:30) {
    if (n %in% c("FT14506", "FT16449", "Z16983", "PH908")) return(n)
    p <- parent_lookup[[n]]
    if (is.na(p) || p == "") break
    n <- p
  }
  if (node %in% node_set || stringr::str_detect(node, "^PH908")) return("PH908")
  NA_character_
}

ph908_branch <- rep(NA_character_, length(terminal_snp))
ph908_branch[is_ph908_wide] <- vapply(terminal_snp[is_ph908_wide], infer_branch, character(1))

tmrca <- map_ph908_only$tmrca[match(terminal_snp, map_ph908_only$node)]
tmrca_source <- ifelse(!is.na(tmrca), "terminal_map", NA_character_)

branch_tmrca <- map_ph908_only$tmrca[match(ph908_branch, map_ph908_only$node)]
idx <- which(is.na(tmrca) & !is.na(branch_tmrca) & is_ph908_wide)
tmrca[idx] <- branch_tmrca[idx]
tmrca_source[idx] <- "branch_fallback"

snp_any_clean <- stringr::str_trim(dplyr::coalesce(clusters$SNP_any, ""))

resolution_class <- dplyr::case_when(
  stringr::str_detect(snp_any_clean, ">") & !stringr::str_detect(snp_any_clean, "\\*\\s*$") ~ "high",
  TRUE ~ "low"
)

cl2 <- clusters |>
  dplyr::mutate(
    major_hg = major_hg,
    terminal_snp = terminal_snp,
    is_ph908_wide = is_ph908_wide,
    ph908_branch = ph908_branch,
    tmrca = tmrca,
    tmrca_source = tmrca_source,
    resolution_class = resolution_class
  ) |>
  dplyr::mutate(
    is_ph908_primary = is_ph908_wide &
      ph908_branch %in% c("FT14506", "FT16449") &
      !terminal_snp %in% c("PH908", "Y3120", "Z16983") &
      resolution_class == "high",
    is_ph908_accurate = is_ph908_wide &
      ph908_branch %in% c("FT14506", "FT16449", "Z16983") &
      !terminal_snp %in% c("PH908", "Y3120"),
    gating_issue = dplyr::case_when(
      is_ph908_primary & resolution_class != "high" ~ "primary_but_not_highres",
      is_ph908_wide & resolution_class == "high" & !is_ph908_primary &
        ph908_branch %in% c("FT14506", "FT16449") &
        !terminal_snp %in% c("PH908", "Y3120", "Z16983") ~ "wide_but_excluded_downstream",
      TRUE ~ "ok"
    )
  )

assert_true(all(!cl2$is_ph908_primary | cl2$is_ph908_wide), "Found is_ph908_primary rows not marked PH908-wide")

annotation_summary <- list(
  total_clusters = nrow(cl2),
  total_i2 = sum(cl2$major_hg == "I2", na.rm = TRUE),
  total_ph908_wide = sum(cl2$is_ph908_wide, na.rm = TRUE),
  total_ph908_primary = sum(cl2$is_ph908_primary, na.rm = TRUE),
  total_ph908_accurate = sum(cl2$is_ph908_accurate, na.rm = TRUE),
  tmrca_source_terminal_map = sum(cl2$tmrca_source == "terminal_map", na.rm = TRUE),
  tmrca_source_branch_fallback = sum(cl2$tmrca_source == "branch_fallback", na.rm = TRUE),
  resolution_high = sum(cl2$resolution_class == "high", na.rm = TRUE),
  resolution_low = sum(cl2$resolution_class == "low", na.rm = TRUE),
  ph908_branch_ft14506 = sum(cl2$ph908_branch == "FT14506", na.rm = TRUE),
  ph908_branch_ft16449 = sum(cl2$ph908_branch == "FT16449", na.rm = TRUE),
  ph908_branch_z16983 = sum(cl2$ph908_branch == "Z16983", na.rm = TRUE),
  ph908_branch_ph908 = sum(cl2$ph908_branch == "PH908", na.rm = TRUE),
  gating_issue_ok = sum(cl2$gating_issue == "ok", na.rm = TRUE)
)

append_facts_lines(FILE_FACTS_PREPROCESSING, c(
  "[annotation_summary]",
  paste0("total_clusters: ", annotation_summary$total_clusters),
  paste0("total_i2: ", annotation_summary$total_i2),
  paste0("total_ph908_wide: ", annotation_summary$total_ph908_wide),
  paste0("total_ph908_primary: ", annotation_summary$total_ph908_primary),
  paste0("total_ph908_accurate: ", annotation_summary$total_ph908_accurate),
  paste0("tmrca_source_terminal_map: ", annotation_summary$tmrca_source_terminal_map),
  paste0("tmrca_source_branch_fallback: ", annotation_summary$tmrca_source_branch_fallback),
  paste0("resolution_high: ", annotation_summary$resolution_high),
  paste0("resolution_low: ", annotation_summary$resolution_low),
  paste0("ph908_branch_ft14506: ", annotation_summary$ph908_branch_ft14506),
  paste0("ph908_branch_ft16449: ", annotation_summary$ph908_branch_ft16449),
  paste0("ph908_branch_z16983: ", annotation_summary$ph908_branch_z16983),
  paste0("ph908_branch_ph908: ", annotation_summary$ph908_branch_ph908),
  paste0("gating_issue_ok: ", annotation_summary$gating_issue_ok),
  ""
))

cat("STEP 4: LOAD AND VALIDATE REVIEWED GEOCODING CACHE\n")

assert_true(file.exists(FILE_GEOCODE_CACHE), paste0("Missing geocoding cache: ", FILE_GEOCODE_CACHE))

locs <- cl2 |>
  dplyr::distinct(Location, Region) |>
  dplyr::mutate(
    Location = stringr::str_trim(Location),
    Region = stringr::str_trim(Region)
  ) |>
  dplyr::filter(nz(Location)) |>
  dplyr::mutate(
    settlement = stringr::str_trim(stringr::str_split_fixed(Location, "/", 2)[, 1]),
    municipality = stringr::str_trim(stringr::str_split_fixed(Location, "/", 2)[, 2]),
    has_slash = stringr::str_detect(Location, "/"),
    municipality = ifelse(municipality == "", NA_character_, municipality)
  )

expected_cols <- c(
  "Location", "Region", "settlement", "municipality", "has_slash",
  "query_used", "query_text", "lat", "long", "display_name", "importance", "place_rank",
  "osm_type", "osm_id", "class", "type", "country_code", "country", "state", "county",
  "matched", "note"
)

cache <- readRDS(FILE_GEOCODE_CACHE)

for (nm in setdiff(expected_cols, names(cache))) {
  cache[[nm]] <- NA
}
cache <- cache |>
  dplyr::select(dplyr::all_of(expected_cols))

todo <- locs |>
  dplyr::anti_join(
    cache |>
      dplyr::distinct(Location, Region),
    by = c("Location", "Region")
  )

cat("  Unique locations missing from cache: ", nrow(todo), "\n", sep = "")

if (nrow(todo) > 0) {
  stop(
    paste0(
      "Geocoding cache is incomplete: ",
      nrow(todo),
      " location(s) missing."
    ),
    call. = FALSE
  )
}

cache_problem_rows <- cache |>
  dplyr::mutate(
    coords_present = !is.na(lat) & !is.na(long),
    in_bbox = in_bbox(lat, long),
    is_balkan_country = toupper(country_code) %in% BALKAN_ISO2,
    flagged_problem = dplyr::case_when(
      coords_present & in_bbox & !is_balkan_country ~ TRUE,
      coords_present & !in_bbox & is_balkan_country ~ TRUE,
      TRUE ~ FALSE
    )
  ) |>
  dplyr::filter(flagged_problem)

validation_summary <- tibble::tibble(
  metric = c(
    "unique_location_region_pairs_checked",
    "cache_rows",
    "missing_required_pairs_in_cache",
    "cache_rows_with_coordinates",
    "cache_rows_missing_coordinates",
    "cache_rows_inside_study_bbox",
    "cache_rows_in_balkan_countries",
    "cache_problem_rows"
  ),
  value = c(
    nrow(locs),
    nrow(cache),
    nrow(todo),
    sum(!is.na(cache$lat) & !is.na(cache$long)),
    sum(is.na(cache$lat) | is.na(cache$long)),
    sum(in_bbox(cache$lat, cache$long), na.rm = TRUE),
    sum(toupper(cache$country_code) %in% BALKAN_ISO2, na.rm = TRUE),
    nrow(cache_problem_rows)
  )
)

cl_geo <- cl2 |>
  dplyr::left_join(
    cache |>
      dplyr::mutate(
        in_bbox = in_bbox(lat, long),
        is_balkan_country = toupper(country_code) %in% BALKAN_ISO2
      ),
    by = c("Location", "Region")
  ) |>
  dplyr::mutate(
    exclude_geo_primary = is.na(lat) | is.na(long) | !in_bbox | !is_balkan_country
  )

geocode_summary <- list(
  unique_location_region_pairs_checked = nrow(locs),
  cache_rows = nrow(cache),
  missing_required_pairs_in_cache = nrow(todo),
  cache_rows_with_coordinates = sum(!is.na(cache$lat) & !is.na(cache$long)),
  cache_rows_missing_coordinates = sum(is.na(cache$lat) | is.na(cache$long)),
  cache_rows_inside_study_bbox = sum(in_bbox(cache$lat, cache$long), na.rm = TRUE),
  cache_rows_in_balkan_countries = sum(toupper(cache$country_code) %in% BALKAN_ISO2, na.rm = TRUE),
  cache_problem_rows = nrow(cache_problem_rows)
)

append_facts_lines(FILE_FACTS_PREPROCESSING, c(
  "[geocode_summary]",
  paste0("unique_location_region_pairs_checked: ", geocode_summary$unique_location_region_pairs_checked),
  paste0("cache_rows: ", geocode_summary$cache_rows),
  paste0("missing_required_pairs_in_cache: ", geocode_summary$missing_required_pairs_in_cache),
  paste0("cache_rows_with_coordinates: ", geocode_summary$cache_rows_with_coordinates),
  paste0("cache_rows_missing_coordinates: ", geocode_summary$cache_rows_missing_coordinates),
  paste0("cache_rows_inside_study_bbox: ", geocode_summary$cache_rows_inside_study_bbox),
  paste0("cache_rows_in_balkan_countries: ", geocode_summary$cache_rows_in_balkan_countries),
  paste0("cache_problem_rows: ", geocode_summary$cache_problem_rows),
  ""
))

cat("STEP 5: WRITE ANALYSIS-READY DATASET\n")

validation_df <- cl_geo |>
  dplyr::mutate(
    lat = as.numeric(lat),
    long = as.numeric(long),
    is_balkan_country = toupper(country_code) %in% BALKAN_ISO2
  )

missing_required <- setdiff(REQUIRED_DOWNSTREAM_COLS, names(validation_df))
assert_true(
  length(missing_required) == 0,
  paste0("Required columns missing before final serialization: ", paste(missing_required, collapse = ", "))
)

final_cols_30 <- c(
  "cluster_id",
  "n_members",
  "Surname",
  "Slava",
  "Region",
  "Location",
  "SNP_any",
  "major_hg",
  "terminal_snp",
  "is_ph908_wide",
  "ph908_branch",
  "is_ph908_primary",
  "is_ph908_accurate",
  "gating_issue",
  "tmrca",
  "tmrca_source",
  "resolution_class",
  "settlement",
  "municipality",
  "has_slash",
  "lat",
  "long",
  "display_name",
  "country_code",
  "country",
  "state",
  "county",
  "matched",
  "note",
  "exclude_geo_primary"
)

missing_final_cols <- setdiff(final_cols_30, names(validation_df))
assert_true(
  length(missing_final_cols) == 0,
  paste0("Missing final serialized columns: ", paste(missing_final_cols, collapse = ", "))
)

final_df <- validation_df |>
  dplyr::select(dplyr::all_of(final_cols_30))

assert_true(
  ncol(final_df) == 30,
  paste0("Final serialized dataset must have 30 columns; found ", ncol(final_df))
)

saveRDS(final_df, FILE_PREPROCESSED_RDS)

final_validation_summary <- list(
  final_rows = nrow(final_df),
  final_columns = ncol(final_df),
  excluded_geo_primary = sum(final_df$exclude_geo_primary, na.rm = TRUE),
  non_missing_coordinates = sum(!is.na(final_df$lat) & !is.na(final_df$long)),
  rows_in_balkan_country = sum(validation_df$is_balkan_country, na.rm = TRUE),
  ph908_primary_rows = sum(final_df$is_ph908_primary, na.rm = TRUE)
)

append_facts_lines(FILE_FACTS_PREPROCESSING, c(
  "[final_validation_summary]",
  paste0("final_rows: ", final_validation_summary$final_rows),
  paste0("final_columns: ", final_validation_summary$final_columns),
  paste0("excluded_geo_primary: ", final_validation_summary$excluded_geo_primary),
  paste0("non_missing_coordinates: ", final_validation_summary$non_missing_coordinates),
  paste0("rows_in_balkan_country: ", final_validation_summary$rows_in_balkan_country),
  paste0("ph908_primary_rows: ", final_validation_summary$ph908_primary_rows),
  ""
))

cat("  Final dataset rows: ", final_validation_summary$final_rows, "\n", sep = "")
cat("  Final dataset columns: ", final_validation_summary$final_columns, "\n\n", sep = "")

cat("STEP 6: WRITE FROZEN OUTPUTS\n")

flow_df <- tibble::tibble(
  step_order = c(1, 2, 3, 4, 5, 6),
  step_label = c(
    "Source table",
    "Encoding-QC source table",
    "Surnameâ€“slavaâ€“location clusters",
    "PH908-annotated clusters",
    "Geocoded clusters",
    "Final analysis-ready dataset"
  ),
  n_rows = c(
    nrow(df),
    nrow(df_clean),
    nrow(clusters),
    nrow(cl2),
    nrow(cl_geo),
    final_validation_summary$final_rows
  ),
  n_unique_surnames = c(
    safe_n_distinct(df$Surname),
    safe_n_distinct(df_clean$Surname),
    safe_n_distinct(clusters$Surname),
    safe_n_distinct(cl2$Surname),
    safe_n_distinct(cl_geo$Surname),
    safe_n_distinct(final_df$Surname)
  ),
  n_unique_locations = c(
    safe_n_distinct(df$Location),
    safe_n_distinct(df_clean$Location),
    safe_n_distinct(clusters$Location),
    safe_n_distinct(cl2$Location),
    safe_n_distinct(cl_geo$Location),
    safe_n_distinct(final_df$Location)
  ),
  note = c(
    "Loaded source table",
    "Minimal encoding repair only",
    "Modal SNP per cluster",
    "Broad and strict PH908 gates",
    "Reviewed geocoding cache joined",
    "Final serialized object"
  )
)

geocode_fig_df <- cache |>
  dplyr::mutate(
    coords_present = !is.na(lat) & !is.na(long),
    in_bbox = in_bbox(lat, long),
    is_balkan_country = toupper(country_code) %in% BALKAN_ISO2,
    flagged = dplyr::case_when(
      coords_present & in_bbox & !is_balkan_country ~ TRUE,
      coords_present & !in_bbox & is_balkan_country ~ TRUE,
      TRUE ~ FALSE
    ),
    qc_group = dplyr::case_when(
      flagged ~ "flagged_problem",
      !coords_present ~ "missing_coordinates",
      coords_present & in_bbox & is_balkan_country ~ "accepted",
      TRUE ~ "other_reviewed"
    )
  ) |>
  dplyr::select(
    Location, Region, settlement, municipality, lat, long,
    country_code, country, matched, note, coords_present,
    in_bbox, is_balkan_country, flagged, qc_group
  )

table_preprocessing_supp <- tibble::tibble(
  metric = c(
    "Source rows",
    "Missing Location fields",
    "Missing SNP fields",
    "Unique surnames",
    "Unique location strings",
    "Unique Location Ã— Region pairs checked",
    "Reviewed cache rows",
    "Reviewed cache rows with coordinates",
    "Reviewed cache rows without coordinates",
    "Reviewed cache rows inside study bounding box",
    "Reviewed cache rows in Balkan countries",
    "Reviewed cache problem rows",
    "Final analysis-ready rows",
    "Final analysis-ready columns"
  ),
  value = c(
    raw_summary$n_rows,
    raw_summary$missing_location,
    raw_summary$missing_snp,
    raw_summary$unique_surnames,
    raw_summary$unique_locations,
    geocode_summary$unique_location_region_pairs_checked,
    geocode_summary$cache_rows,
    geocode_summary$cache_rows_with_coordinates,
    geocode_summary$cache_rows_missing_coordinates,
    geocode_summary$cache_rows_inside_study_bbox,
    geocode_summary$cache_rows_in_balkan_countries,
    geocode_summary$cache_problem_rows,
    final_validation_summary$final_rows,
    final_validation_summary$final_columns
  ),
  note = c(
    "Loaded source table",
    "Expected to be zero",
    "Expected to be zero",
    "Distinct source-table surnames",
    "Distinct source-table location strings",
    "Pairs checked against reviewed cache",
    "Rows present in reviewed cache",
    "Rows with non-missing latitude and longitude",
    "Rows lacking coordinates in reviewed cache",
    "Rows within fixed study bounding box",
    "Rows within Balkan-country whitelist",
    "Rows flagged during cache validation",
    "Saved analysis-ready dataset rows",
    "Saved analysis-ready dataset columns"
  )
)

readr::write_csv(flow_df, FILE_FIGUREDATA_PREPROCESSING_FLOW)
readr::write_csv(geocode_fig_df, FILE_FIGUREDATA_PREPROCESSING_GEOCODES)
readr::write_csv(table_preprocessing_supp, FILE_TABLE_PREPROCESSING_SUPP)

results_preprocessing <- list(
    script = "01_preprocess_registry_data.R",
  seed = SEED_MAIN,
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  input_files = list(
    dnk_reconstructed = FILE_DNK_RECON_CSV,
    ph908_branch_map = FILE_PH908_MAP_TSV,
    ph908_node_age_tiers = FILE_PH908_NODE_AGE_TIERS_TSV,
    geocode_cache = FILE_GEOCODE_CACHE,
    branch_ages = FILE_BRANCH_AGES
  ),
  raw_source = df,
  clean_source = df_clean,
  clusters = clusters,
  ph908_map = map_ph908_only,
  node_age_tiers = age_df,
  annotated_clusters = cl2,
  geocoded_clusters = cl_geo,
  preprocessed_dataset = final_df,
  summaries = list(
    raw_summary = raw_summary,
    encoding_summary = encoding_summary,
    cluster_summary = cluster_summary,
    annotation_summary = annotation_summary,
    geocode_summary = geocode_summary,
    final_validation_summary = final_validation_summary
  ),
  outputs = list(
    preprocessed_rds = FILE_PREPROCESSED_RDS,
    results_rds = FILE_RESULTS_PREPROCESSING,
    table_preprocessing_supp = FILE_TABLE_PREPROCESSING_SUPP,
    figuredata_flow = FILE_FIGUREDATA_PREPROCESSING_FLOW,
    figuredata_geocodes = FILE_FIGUREDATA_PREPROCESSING_GEOCODES,
    facts = FILE_FACTS_PREPROCESSING
  )
)

saveRDS(results_preprocessing, FILE_RESULTS_PREPROCESSING)

append_facts_lines(FILE_FACTS_PREPROCESSING, c(
  "[wording_ready_lines]",
  paste0(
    "methods_source_sentence: The source table contained ",
    raw_summary$n_rows,
    " rows, with no missing Location fields and no missing SNP fields."
  ),
  paste0(
    "methods_geocode_sentence: A total of ",
    geocode_summary$unique_location_region_pairs_checked,
    " unique Location Ã— Region pairs were checked against the reviewed geocoding cache, and no required pairs were absent."
  ),
  paste0(
    "methods_final_sentence: The final analysis-ready dataset contained ",
    final_validation_summary$final_rows,
    " rows and ",
    final_validation_summary$final_columns,
    " columns, and all required downstream columns were present at the final validation step."
  ),
  paste0(
    "figure_s1_caption_line: Supplementary preprocessing flow and geocoding outputs should be rendered from figuredata_preprocessing_flow.csv and figuredata_preprocessing_geocodes.csv."
  ),
  paste0(
    "table_s1_caption_line: table_preprocessing_supp.csv summarizes source-table, reviewed-cache validation, and final dataset metrics for the deterministic preprocessing workflow."
  ),
  ""
))

cat("============================================================\n")
cat("PREPROCESSING COMPLETE\n")
cat("============================================================\n")
cat("Final rows: ", final_validation_summary$final_rows, "\n", sep = "")
cat("Final columns: ", final_validation_summary$final_columns, "\n", sep = "")
cat("Results object: ", FILE_RESULTS_PREPROCESSING, "\n", sep = "")
cat("Facts file: ", FILE_FACTS_PREPROCESSING, "\n", sep = "")
cat("Supplement table: ", FILE_TABLE_PREPROCESSING_SUPP, "\n", sep = "")
cat("Figuredata flow: ", FILE_FIGUREDATA_PREPROCESSING_FLOW, "\n", sep = "")
cat("Figuredata geocodes: ", FILE_FIGUREDATA_PREPROCESSING_GEOCODES, "\n", sep = "")
cat("Serialized dataset: ", FILE_PREPROCESSED_RDS, "\n", sep = "")
cat("============================================================\n")
