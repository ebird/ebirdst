ebirdst_trends_runs <- readr::read_csv(
  "data-raw/ebird-trends_2021_model-summary_production.csv",
  na = ""
)
ebirdst_trends_runs <- dplyr::arrange(ebirdst_trends_runs, species_code, season)
usethis::use_data(ebirdst_trends_runs, overwrite = TRUE)
