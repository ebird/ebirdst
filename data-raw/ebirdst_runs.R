library(tidyverse)
library(auk)
library(lubridate)
library(jsonlite)
library(glue)

pred_year <- file.path("data-raw", "config_status.json") |>
  read_json(simplifyVector = TRUE) |>
  pluck("SRD_PRED_YEAR")

# s3 species
s3_bucket <- Sys.getenv("EBIRDST_S3_BUCKET")
species_codes <- glue("aws s3 ls {s3_bucket}/{pred_year}/") |>
  paste("awk '{print $2}'", sep = " | ") |>
  system(intern = TRUE) |>
  str_remove_all("/")

# reviews
gs_key <- Sys.getenv("EBIRDST_STATUS_GS_KEY")
runs <- glue("https://docs.google.com/spreadsheets/d/{gs_key}/",
             "export?format=csv") |>
  read_csv(show_col_types = FALSE) |>
  rename_with(tolower) |>
  filter(status == "REVIEWED", full_year_quality > 0) |>
  rename(resident_quality = full_year_quality) |>
  filter(!is.na(species_code))

# species passing review but missing from s3
filter(runs, !species_code %in% species_codes)
# species on s3 but not passing review
setdiff(species_codes, runs$species_code)

# correctly na season dates
seasons <- c("breeding", "nonbreeding",
             "prebreeding_migration", "postbreeding_migration",
             "resident")
is_resident <- runs$summarize_as_resident
for (s in seasons) {
  s_fail <- runs[[paste0(s, "_quality")]] == 0 |
    is.na(runs[[paste0(s, "_quality")]])
  runs[[paste0(s, "_start")]][s_fail] <- NA_character_
  runs[[paste0(s, "_end")]][s_fail] <- NA_character_
  if (s == "resident") {
    runs[[paste0(s, "_start")]][!is_resident] <- NA_character_
    runs[[paste0(s, "_end")]][!is_resident] <- NA_character_
    runs[[paste0(s, "_quality")]][!is_resident] <- NA_character_
  } else {
    # remove any seasonal information for residents
    runs[[paste0(s, "_start")]][is_resident] <- NA_character_
    runs[[paste0(s, "_end")]][is_resident] <- NA_character_
    runs[[paste0(s, "_quality")]][is_resident] <- NA_character_
    runs[[paste0(s, "_range_modeled")]][is_resident] <- NA_character_
  }
}

# default residents to full year
fy_resident <- is_resident & runs$resident_quality > 0 &
  is.na(runs$resident_start) & is.na(runs$resident_end)
runs$resident_start[fy_resident] <- "01-04"
runs$resident_end[fy_resident] <- "12-28"

# clean up
convert_to_date <- function(x) {
  ymd(ifelse(is.na(x), NA_character_, paste0(pred_year, "-", x)))
}
ebirdst_runs <- runs |>
  select(-common_name) |>
  inner_join(ebird_taxonomy, by = "species_code") |>
  arrange(taxon_order) |>
  mutate(across(ends_with("start"), convert_to_date),
         across(ends_with("end"), convert_to_date)) |>
  select(species_code, scientific_name, common_name,
         is_resident = summarize_as_resident,
         breeding_quality, breeding_start, breeding_end,
         nonbreeding_quality, nonbreeding_start, nonbreeding_end,
         postbreeding_migration_quality,
         postbreeding_migration_start, postbreeding_migration_end,
         prebreeding_migration_quality,
         prebreeding_migration_start, prebreeding_migration_end,
         resident_quality, resident_start, resident_end)

# trends runs
trends <- read_csv("data-raw/ebird-trends_runs_2022.csv",
                   show_col_types = FALSE) |>
  mutate(species_code = case_match(species_code, "norgos2" ~ "norgos",
                                   .default = species_code)) |>
  transmute(has_trends = TRUE,
            species_code,
            trends_season = season,
            trends_region = modeled_region,
            trends_start_year = start_year,
            trends_end_year = end_year,
            trends_start_date = start_date,
            trends_end_date = end_date,
            rsquared, beta0)

# combine
ebirdst_runs <- left_join(ebirdst_runs, trends, by = "species_code") |>
  mutate(has_trends = coalesce(has_trends, FALSE)) |>
  arrange(species_code)

# add a row for yebsap example
ebirdst_runs <- ebirdst_runs |>
  filter(species_code == "yebsap") |>
  mutate(species_code = "yebsap-example") |>
  bind_rows(ebirdst_runs)

usethis::use_data(ebirdst_runs, overwrite = TRUE)
