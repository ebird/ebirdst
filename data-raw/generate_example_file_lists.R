library(ebirdst)
library(glue)
library(purrr)
library(readr)
library(stringr)

status_year <- ebirdst_version()[["status_version_year"]]
trends_year <- ebirdst_version()[["trends_version_year"]]

repo <- "https://raw.githubusercontent.com/ebird/ebirdst_example-data/main/"
file.path(repo, glue("example-data/{trends_year}/file-list.txt")) |>
  read_lines() |>
  keep(str_detect, pattern = "trends/") |>
  write_lines("inst/extdata/example-data_file-list_trends.txt")
file.path(repo, glue("example-data/{status_year}/file-list.txt")) |>
  read_lines() |>
  discard(str_detect, pattern = "trends/") |>
  write_lines("inst/extdata/example-data_file-list_status.txt")
