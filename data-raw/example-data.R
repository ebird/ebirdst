library(readr)
library(purrr)
library(stringr)

repo <- "https://raw.githubusercontent.com/ebird/ebirdst_example-data/main/"
read_lines(file.path(repo, "example-data/2022/file-list.txt")) |>
  keep(str_detect, pattern = "trends/") |>
  write_lines("inst/extdata/example-data_file-list_trends.txt")
read_lines(file.path(repo, "example-data/2023/file-list.txt")) |>
  discard(str_detect, pattern = "trends/") |>
  write_lines("inst/extdata/example-data_file-list_status.txt")
