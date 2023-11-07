library(tidyverse)
library(jsonlite)

# feature set - status + trends
pred_list <- file.path("data-raw", "config_status.json") %>%
  read_json(simplifyVector = TRUE) %>%
  pluck("PREDICTOR_LIST") %>%
  # add in trends predictors
  c("longitude", "latitude", .,
    "mcd12q1_lccs2_c9_ed", "mcd12q1_lccs2_c9_pland")

# categories
p <- read_csv("data-raw/ebirdst_features_2022 - predictors.csv") %>%
  mutate(row = row_number())

# don't need to split
p_nosplit <- filter(p, !str_detect(predictor, "\\{"))

# split to generate all predictors
p_split <- filter(p, str_detect(predictor, "\\{")) %>%
  mutate(suffix = str_extract(predictor, "\\{.*\\}") %>%
           str_remove_all("[\\{\\}]") %>%
           map(~ data.frame(suffix = str_split_1(., "/"))),
         prefix = str_remove(predictor, "\\{.*\\}")) %>%
  unnest(suffix) %>%
  mutate(predictor = paste0(prefix, suffix)) %>%
  select(-prefix, -suffix)

# only keep predictors we use in status or trends models
ebirdst_predictors <- bind_rows(p_nosplit, p_split) %>%
  arrange(row) %>%
  select(-row) %>%
  filter(predictor %in% pred_list)

usethis::use_data(ebirdst_predictors, overwrite = TRUE)

# predictor datasets
ebirdst_predictor_datasets <- read_csv("data-raw/ebirdst_features_2022 - predictor_datasets.csv") %>%
  filter(str_detect(predictor, "\\{") | predictor %in% pred_list)

usethis::use_data(ebirdst_predictor_datasets, overwrite = TRUE)
