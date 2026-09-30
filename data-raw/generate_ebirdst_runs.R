library(auk)
library(dplyr)
library(ebirdst)
library(glue)
library(googlesheets4)
library(jsonlite)
library(lubridate)
library(readr)
library(stringr)

# release/prediction year
prediction_year <- ebirdst_version()[["status_version_year"]]
# S3 bucket for data products
s3_bucket <- Sys.getenv("EBIRDST_S3_BUCKET")

# species with results uploaded to the S3 bucket
species_codes <- glue("aws s3 ls {s3_bucket}/{prediction_year}/") |>
  paste("awk '{print $2}'", sep = " | ") |>
  system(intern = TRUE) |>
  str_remove_all("/")

# reviews
gs_key <- Sys.getenv("EBIRDST_STATUS_GS_KEY")
status_review <- glue(
  "https://docs.google.com/spreadsheets/d/{gs_key}/",
  "export?format=csv"
) |>
  read_csv(show_col_types = FALSE) |>
  rename_with(tolower) |>
  filter(status == "REVIEWED", full_year_quality > 0) |>
  rename(resident_quality = full_year_quality) |>
  filter(!is.na(species_code))

# species passing review but missing from S3
filter(status_review, !species_code %in% species_codes)
# species on S3 but not passing review
setdiff(species_codes, status_review$species_code)

# correctly na season dates
seasons <- c(
  "breeding",
  "nonbreeding",
  "prebreeding_migration",
  "postbreeding_migration",
  "resident"
)
is_resident <- status_review$summarize_as_resident
for (s in seasons) {
  s_fail <- status_review[[paste0(s, "_quality")]] == 0 |
    is.na(status_review[[paste0(s, "_quality")]])
  status_review[[paste0(s, "_start")]][s_fail] <- NA_character_
  status_review[[paste0(s, "_end")]][s_fail] <- NA_character_
  if (s == "resident") {
    status_review[[paste0(s, "_start")]][!is_resident] <- NA_character_
    status_review[[paste0(s, "_end")]][!is_resident] <- NA_character_
    status_review[[paste0(s, "_quality")]][!is_resident] <- NA_character_
  } else {
    # remove any seasonal information for residents
    status_review[[paste0(s, "_start")]][is_resident] <- NA_character_
    status_review[[paste0(s, "_end")]][is_resident] <- NA_character_
    status_review[[paste0(s, "_quality")]][is_resident] <- NA_character_
    status_review[[paste0(s, "_range_modeled")]][is_resident] <- NA_character_
  }
}

# default residents to full year
fy_resident <- is_resident &
  status_review$resident_quality > 0 &
  is.na(status_review$resident_start) &
  is.na(status_review$resident_end)
status_review$resident_start[fy_resident] <- "01-04"
status_review$resident_end[fy_resident] <- "12-28"

# clean up
convert_to_date <- function(x) {
  ymd(ifelse(is.na(x), NA_character_, paste0(prediction_year, "-", x)))
}
status_review <- status_review |>
  select(!c(common_name, taxon_order)) |>
  inner_join(ebird_taxonomy, by = "species_code") |>
  arrange(taxonomic_order) |>
  mutate(
    across(ends_with("start"), convert_to_date),
    across(ends_with("end"), convert_to_date),
    status_version_year = ebirdst_version()[["status_version_year"]]
  ) |>
  select(
    species_code,
    taxon_concept_id,
    scientific_name,
    common_name,
    taxonomic_order,
    is_resident = summarize_as_resident,
    breeding_quality,
    breeding_start,
    breeding_end,
    nonbreeding_quality,
    nonbreeding_start,
    nonbreeding_end,
    postbreeding_migration_quality,
    postbreeding_migration_start,
    postbreeding_migration_end,
    prebreeding_migration_quality,
    prebreeding_migration_start,
    prebreeding_migration_end,
    resident_quality,
    resident_start,
    resident_end,
    status_version_year
  )

# trends reviews
trends_review <- read_csv(
  "data-raw/ebird-trends_runs_2022.csv",
  show_col_types = FALSE
) |>
  mutate(
    species_code = recode_values(
      species_code,
      "norgos2" ~ "norgos",
      default = species_code
    )
  ) |>
  mutate(
    has_trends = TRUE,
    trends_version_year = ebirdst_version()[["trends_version_year"]]
  ) |>
  select(
    has_trends,
    species_code,
    trends_season = season,
    trends_region = modeled_region,
    trends_start_year = start_year,
    trends_end_year = end_year,
    trends_start_date = start_date,
    trends_end_date = end_date,
    rsquared,
    beta0,
    trends_version_year
  )

# species codes in trends but not status
# these are likely taxonomy changes since trends is quite outdated
setdiff(trends_review$species_code, status_review$species_code)

# combine
ebirdst_runs <- left_join(
  status_review,
  trends_review,
  by = "species_code"
) |>
  mutate(has_trends = coalesce(has_trends, FALSE)) |>
  arrange(taxonomic_order)

# add a row for yebsap example
ebirdst_runs <- ebirdst_runs |>
  filter(species_code == "yebsap") |>
  mutate(species_code = "yebsap-example") |>
  bind_rows(ebirdst_runs)

usethis::use_data(ebirdst_runs, overwrite = TRUE)
