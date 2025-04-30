library(jsonlite)
library(stringi)
library(tidyverse)

# feature set - status + trends
pred_list <- file.path("data-raw", "config_status.json") |>
  read_json(simplifyVector = TRUE) |>
  pluck("PREDICTOR_LIST")
# add in trends predictors
pred_list <- c("longitude", "latitude", pred_list,
               "mcd12q1_lccs2_c9_ed", "mcd12q1_lccs2_c9_pland")

# categories
p <- read_csv("data-raw/ebirdst_features_2023 - predictors.csv") |>
  mutate(row = row_number())

# don't need to split
p_nosplit <- filter(p, !str_detect(predictor, "\\{"))

# split to generate all predictors
p_split <- filter(p, str_detect(predictor, "\\{")) |>
  mutate(suffix = str_extract(predictor, "\\{.*\\}") |>
           str_remove_all("[\\{\\}]") |>
           map(~ data.frame(suffix = str_split_1(., "/"))),
         prefix = str_remove(predictor, "\\{.*\\}")) |>
  unnest(suffix) |>
  mutate(predictor = paste0(prefix, suffix),
         label = paste(label,
                       recode(suffix,
                              median = "(median)",
                              mean = "(mean)",
                              sd = "(SD)",
                              pland = "(% cover)",
                              ed = "(edge density)")
         ),
         label = str_replace(label, "µg/L", "g/1000L")) |>
  select(-prefix, -suffix)

# only keep predictors we use in status or trends models
ebirdst_predictors <- bind_rows(p_nosplit, p_split) |>
  arrange(row) |>
  select(-row) |>
  filter(predictor %in% pred_list) |>
  as_tibble()

usethis::use_data(ebirdst_predictors, overwrite = TRUE)

# predictor datasets
ebirdst_predictor_descriptions <- read_csv("data-raw/ebirdst_features_2023 - predictor_datasets.csv") |>
  select(!index) |>
  filter(str_detect(predictor, "\\{") | predictor %in% pred_list) |>
  as_tibble()

usethis::use_data(ebirdst_predictor_descriptions, overwrite = TRUE)
