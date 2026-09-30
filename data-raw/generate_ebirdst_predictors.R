library(dplyr)
library(ebirdst)
library(glue)
library(jsonlite)
library(purrr)
library(readr)
library(stringi)
library(stringr)
library(tidyr)

# release/prediction year
prediction_year <- ebirdst_version()[["status_version_year"]]
# S3 bucket for data products
s3_bucket <- Sys.getenv("EBIRDST_S3_BUCKET")

# grab an example config.json file
example_species <- "yebsap"
src_uri <- glue("{s3_bucket}/{prediction_year}/{example_species}/config.json")
dst_uri <- file.path(tempdir(), "config.json")
glue("aws s3 cp {src_uri} {dst_uri}") |>
  system()
pred_list <- read_json(dst_uri, simplifyVector = TRUE)[["PREDICTOR_LIST"]]
# add in trends predictors
pred_list <- c(
  "longitude",
  "latitude",
  pred_list,
  "mcd12q1_lccs2_c9_ed",
  "mcd12q1_lccs2_c9_pland"
)

# categories
gs_key <- Sys.getenv("EBIRDST_FEATURES_GS_KEY")
p <- glue(
  "https://docs.google.com/spreadsheets/d/{gs_key}/",
  "export?format=csv&sheet=predictors"
) |>
  read_csv(show_col_types = FALSE) |>
  mutate(row = row_number())

# don't need to split
p_nosplit <- filter(p, !str_detect(predictor, "\\{"))

# split to generate all predictors
p_split <- filter(p, str_detect(predictor, "\\{")) |>
  mutate(
    suffix = str_extract(predictor, "\\{.*\\}") |>
      str_remove_all("[\\{\\}]") |>
      map(~ data.frame(suffix = str_split_1(., "/"))),
    prefix = str_remove(predictor, "\\{.*\\}")
  ) |>
  unnest(suffix) |>
  mutate(
    predictor = paste0(prefix, suffix),
    label = paste(
      label,
      recode(
        suffix,
        median = "(median)",
        mean = "(mean)",
        sd = "(SD)",
        pland = "(% cover)",
        ed = "(edge density)"
      )
    ),
    label = str_replace(label, "µg/L", "g/1000L")
  ) |>
  select(-prefix, -suffix)

# only keep predictors we use in status or trends models
ebirdst_predictors <- bind_rows(p_nosplit, p_split) |>
  arrange(row) |>
  select(-row) |>
  filter(predictor %in% pred_list) |>
  as_tibble()

usethis::use_data(ebirdst_predictors, overwrite = TRUE)

# predictor datasets
ebirdst_predictor_descriptions <- glue(
  "https://docs.google.com/spreadsheets/d/{gs_key}/",
  "export?format=csv&sheet=predictor_datasets"
) |>
  read_csv(show_col_types = FALSE) |>
  filter(str_detect(predictor, "\\{") | predictor %in% pred_list) |>
  as_tibble()

usethis::use_data(ebirdst_predictor_descriptions, overwrite = TRUE)
