context("Fetch")

skip_on_cran()

# these keys are independent of any real download, so no network or access
# key is required to test them
keys <- c(
  "2023/woothr/config.json",
  "2023/woothr/weekly/woothr_abundance_median_3km_2023.tif",
  "2023/woothr/weekly/woothr_proportion-population_median_3km_2023.tif",
  "2023/woothr/weekly/woothr_occurrence_median_3km_2023.tif",
  "2023/woothr/weekly/woothr_count_median_3km_2023.tif",
  "2023/woothr/ranges/woothr_range_smooth_27km_2023.gpkg",
  "2023/woothr/regional_stats.csv",
  "2023/woothr/pis/woothr_pi_occurrence_gsw_c2_pland_27km_2023.tif",
  "2023/woothr/ppms/woothr_ppm_binary-f1_mean_27km_2023.tif"
)

test_that("status_key() and trends_key()", {
  expect_equal(status_key("woothr"), "2023/woothr")
  expect_equal(status_key("woothr", "config.json"), "2023/woothr/config.json")
  expect_equal(
    status_key("woothr", "weekly", c("a.tif", "b.tif")),
    c("2023/woothr/weekly/a.tif", "2023/woothr/weekly/b.tif")
  )
  expect_equal(trends_key("woothr", "trends", "a.parquet"), "2022/woothr/trends/a.parquet")
})


test_that("resolve_species()", {
  expect_equal(resolve_species("woothr"), "woothr")
  expect_equal(resolve_species("Wood Thrush"), "woothr")
  expect_error(resolve_species("XXXX"), "does not correspond")
})


test_that("select_status_keys() default and flag-based selection", {
  default <- select_status_keys(keys)
  expect_true(all(grepl("config.json$|_abundance_|_proportion-population_", default)))
  expect_false(any(grepl("_occurrence_|_count_|ranges|regional_stats|pis|ppms", default)))

  expect_true(any(grepl("_occurrence_", select_status_keys(keys, download_occurrence = TRUE))))
  expect_true(any(grepl("_count_", select_status_keys(keys, download_count = TRUE))))
  expect_true(any(grepl("/ranges/", select_status_keys(keys, download_ranges = TRUE))))
  expect_true(any(grepl(
    "regional_stats.csv",
    select_status_keys(keys, download_regional = TRUE)
  )))
  expect_true(any(grepl("/pis/", select_status_keys(keys, download_pis = TRUE))))
  expect_true(any(grepl("/ppms/", select_status_keys(keys, download_ppms = TRUE))))

  # config is always retained, even when no flags are set
  none <- select_status_keys(keys, download_abundance = FALSE)
  expect_equal(none, keys[grepl("config.json$", keys)])
})


test_that("select_status_keys() download_all", {
  all_keys <- select_status_keys(keys, download_all = TRUE)
  expect_equal(sort(all_keys), sort(keys))
})


test_that("select_status_keys() pattern argument", {
  filtered <- select_status_keys(keys, download_all = TRUE, pattern = "27km")
  expect_true(all(grepl("config.json$|27km", filtered)))
  expect_true(any(grepl("27km", filtered)))

  expect_error(
    select_status_keys(keys, download_all = TRUE, pattern = "zzzzz"),
    "No files matched pattern"
  )
})


test_that("ensure_data_dir()", {
  tmp <- withr::local_tempdir()
  nested <- file.path(tmp, "a", "b")
  ensure_data_dir(nested)
  expect_true(dir.exists(nested))

  # a path nested under a file rather than a directory can't be created
  blocker <- file.path(tmp, "blocker")
  file.create(blocker)
  expect_error(ensure_data_dir(file.path(blocker, "nested")))
})


test_that("fetch_data() downloads example data and returns normalized paths", {
  tmp <- withr::local_tempdir()
  key <- status_key("yebsap-example", "config.json")
  local_path <- suppressMessages(fetch_data(key, path = tmp, show_progress = FALSE))
  expect_equal(local_path, normalizePath(file.path(tmp, key)))
  expect_true(file.exists(local_path))
})


test_that("fetch_data() reports existing files only when requested", {
  tmp <- withr::local_tempdir()
  key <- status_key("yebsap-example", "config.json")
  suppressMessages(fetch_data(key, path = tmp, show_progress = FALSE))

  expect_silent(fetch_data(key, path = tmp, show_progress = FALSE))
  expect_message(
    fetch_data(key, path = tmp, show_progress = FALSE, report_existing = TRUE),
    "already exists"
  )
})


test_that("fetch_data() re-downloads on force = TRUE", {
  tmp <- withr::local_tempdir()
  key <- status_key("yebsap-example", "config.json")
  suppressMessages(fetch_data(key, path = tmp, show_progress = FALSE))
  expect_message(
    fetch_data(key, path = tmp, force = TRUE, show_progress = TRUE),
    "Downloading file"
  )
})


test_that("fetch_data() errors, with hint, for a nonexistent key", {
  tmp <- withr::local_tempdir()
  bad_key <- status_key("yebsap-example", "does-not-exist.tif")
  expect_error(
    fetch_data(bad_key, path = tmp, show_progress = FALSE, hint = "some hint"),
    "some hint"
  )
})
