library(tidyverse)
library(auk)
library(lubridate)
library(jsonlite)
library(glue)

pred_year <- file.path("data-raw", "config.json") %>%
  read_json(simplifyVector = TRUE) %>%
  pluck("SRD_PRED_YEAR")

# s3 species
s3_bucket <- Sys.getenv("EBIRDST_S3_BUCKET")
species_codes <- glue("aws s3 ls {s3_bucket}/{pred_year}/") %>%
  paste("awk '{print $2}'", sep = " | ") %>%
  system(intern = TRUE) %>%
  str_remove_all("/")

# reviews
gs_key <- Sys.getenv("EBIRDST_GS_KEY")
runs <- glue("https://docs.google.com/spreadsheets/d/{gs_key}/",
                "export?format=csv") %>%
  read_csv(show_col_types = FALSE) %>%
  rename_with(tolower) %>%
  filter(full_year_quality > 0) %>%
  rename(resident_quality = full_year_quality) %>%
  filter(!is.na(species_code)) %>%
  # these were included mistakely
  filter(!species_code %in% c("gobwar2", "auspip1"))

# species passing review but missing from s3
filter(runs, !species_code %in% species_codes)
# species on s3 but not passing review
setdiff(species_codes, runs$species_code)

# correctly na season dates
seasons <- c("breeding", "nonbreeding",
             "prebreeding_migration", "postbreeding_migration",
             "resident")
for (s in seasons) {
  s_fail <- runs[[paste0(s, "_quality")]] == 0 |
    is.na(runs[[paste0(s, "_quality")]])
  runs[[paste0(s, "_start")]][s_fail] <- NA_character_
  runs[[paste0(s, "_end")]][s_fail] <- NA_character_
  if (s == "resident") {
    runs[[paste0(s, "_start")]][!runs$resident] <- NA_character_
    runs[[paste0(s, "_end")]][!runs$resident] <- NA_character_
  } else {
    runs[[paste0(s, "_start")]][runs$resident] <- NA_character_
    runs[[paste0(s, "_end")]][runs$resident] <- NA_character_
  }
}

# default residents to full year
fy_resident <- runs$resident & runs$resident_quality > 0 &
  is.na(runs$resident_start) & is.na(runs$resident_end)
runs$resident_start[fy_resident] <- "01-04"
runs$resident_end[fy_resident] <- "12-28"

# clean up
convert_to_date <- function(x) {
  ymd(ifelse(is.na(x), NA_character_, paste0(pred_year, "-", x)))
}
ebirdst_runs <- runs %>%
  select(-common_name) %>%
  inner_join(ebird_taxonomy, by = "species_code") %>%
  arrange(taxon_order) %>%
  mutate(across(ends_with("start"), convert_to_date),
         across(ends_with("end"), convert_to_date)) %>%
  select(species_code, scientific_name, common_name, resident,
         breeding_quality, breeding_range_modeled,
         breeding_start, breeding_end,
         nonbreeding_quality, nonbreeding_range_modeled,
         nonbreeding_start, nonbreeding_end,
         postbreeding_migration_quality,
         postbreeding_migration_range_modeled,
         postbreeding_migration_start, postbreeding_migration_end,
         prebreeding_migration_quality,
         prebreeding_migration_range_modeled,
         prebreeding_migration_start, prebreeding_migration_end,
         resident_quality, resident_start, resident_end)

usethis::use_data(ebirdst_runs, overwrite = TRUE)
