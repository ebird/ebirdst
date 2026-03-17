context("Data management")

skip_on_cran()


# format_size ----

test_that("format_size()", {
  expect_equal(ebirdst:::format_size(500), "500 B")
  expect_equal(ebirdst:::format_size(1500), "1.5 KB")
  expect_equal(ebirdst:::format_size(1.5e6), "1.5 MB")
  expect_equal(ebirdst:::format_size(2.3e9), "2.3 GB")
})



# ebirdst_data_inventory ----

test_that("ebirdst_data_inventory() returns empty tibble for missing path", {
  inv <- ebirdst_data_inventory("/nonexistent/path/xyz")
  expect_s3_class(inv, "ebirdst_inventory")
  expect_s3_class(inv, "tbl_df")
  expect_equal(nrow(inv), 0)
  expect_named(inv, c("species_code", "common_name", "scientific_name",
                       "version_year", "dataset", "n_files", "size_mb"))
  expect_type(inv$version_year, "integer")
  expect_type(inv$n_files, "integer")
  expect_type(inv$size_mb, "double")
  expect_type(inv$dataset, "character")
})

test_that("ebirdst_data_inventory() returns empty tibble for path with no year dirs", {
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "not_a_year"))
  inv <- ebirdst_data_inventory(tmp)
  expect_s3_class(inv, "tbl_df")
  expect_equal(nrow(inv), 0)
})

test_that("ebirdst_data_inventory() validates path argument", {
  expect_error(ebirdst_data_inventory(123))
  expect_error(ebirdst_data_inventory(c("a", "b")))
})

test_that("ebirdst_data_inventory() returns empty tibble for empty species directory", {
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "2023", "emptysp"), recursive = TRUE)
  inv <- ebirdst_data_inventory(tmp)
  expect_s3_class(inv, "tbl_df")
  expect_equal(nrow(inv), 0)
})

test_that("ebirdst_data_inventory() detects status data", {
  tmp <- withr::local_tempdir()
  sp_dir <- file.path(tmp, "2023", "woothr")
  dir.create(sp_dir, recursive = TRUE)
  writeLines("data", file.path(sp_dir, "config.json"))
  writeLines("data", file.path(sp_dir, "abundance.tif"))

  inv <- ebirdst_data_inventory(tmp)
  expect_equal(nrow(inv), 1)
  expect_equal(inv$dataset, "status")
  expect_equal(inv$n_files, 2L)
  expect_gt(inv$size_mb, 0)
})

test_that("ebirdst_data_inventory() detects trends data", {
  tmp <- withr::local_tempdir()
  trends_dir <- file.path(tmp, "2022", "woothr", "trends")
  dir.create(trends_dir, recursive = TRUE)
  writeLines("data", file.path(trends_dir, "woothr_trends.parquet"))

  inv <- ebirdst_data_inventory(tmp)
  expect_equal(nrow(inv), 1)
  expect_equal(inv$dataset, "trends")
  expect_equal(inv$n_files, 1L)
})

test_that("ebirdst_data_inventory() reports status and trends as separate rows", {
  tmp <- withr::local_tempdir()
  sp_dir <- file.path(tmp, "2023", "woothr")
  dir.create(file.path(sp_dir, "trends"), recursive = TRUE)
  writeLines("data", file.path(sp_dir, "config.json"))
  writeLines("data", file.path(sp_dir, "trends", "woothr_trends.parquet"))

  inv <- ebirdst_data_inventory(tmp)
  expect_equal(nrow(inv), 2)
  expect_setequal(inv$dataset, c("status", "trends"))
  expect_true(all(inv$species_code == "woothr"))
})

test_that("ebirdst_data_inventory() assigns Data Coverage common name", {
  tmp <- withr::local_tempdir()
  cov_dir <- file.path(tmp, "2023", "data_coverage")
  dir.create(cov_dir, recursive = TRUE)
  writeLines("data", file.path(cov_dir, "coverage.tif"))

  inv <- ebirdst_data_inventory(tmp)
  expect_equal(inv$common_name, "Data Coverage")
  expect_true(is.na(inv$scientific_name))
})

test_that("ebirdst_data_inventory() leaves names NA for unrecognized species", {
  tmp <- withr::local_tempdir()
  sp_dir <- file.path(tmp, "2023", "notaspecies")
  dir.create(sp_dir, recursive = TRUE)
  writeLines("data", file.path(sp_dir, "file.tif"))

  inv <- ebirdst_data_inventory(tmp)
  expect_true(is.na(inv$common_name))
  expect_true(is.na(inv$scientific_name))
})

test_that("ebirdst_data_inventory() returns data for downloaded example species", {
  inv <- ebirdst_data_inventory()
  expect_s3_class(inv, "ebirdst_inventory")
  expect_s3_class(inv, "tbl_df")
  # yebsap-example should be present (downloaded in setup.R)
  expect_true("yebsap-example" %in% inv$species_code)
  expect_named(inv, c("species_code", "common_name", "scientific_name",
                       "version_year", "dataset", "n_files", "size_mb"))
  # files should have been found
  yeb <- inv[inv$species_code == "yebsap-example", ]
  expect_true(all(yeb$n_files > 0))
  expect_true(all(yeb$size_mb > 0))
  # dataset should be status and/or trends
  expect_true(all(yeb$dataset %in% c("status", "trends")))
  # version_year should be a plausible 4-digit year
  expect_true(all(yeb$version_year >= 2021L))
  # sorted by version_year, species_code, dataset
  if (nrow(inv) > 1) {
    expect_identical(inv, dplyr::arrange(inv, .data$version_year,
                                         .data$species_code, .data$dataset))
  }
})


# print.ebirdst_inventory ----

test_that("print.ebirdst_inventory() header summarises species and packages", {
  tmp <- withr::local_tempdir()
  sp_dir <- file.path(tmp, "2023", "woothr")
  dir.create(sp_dir, recursive = TRUE)
  writeLines("data", file.path(sp_dir, "config.json"))

  out <- capture.output(print(ebirdst_data_inventory(tmp)))
  expect_match(out[1], "eBird Status and Trends data: 1 species, 1 package")
})

test_that("print.ebirdst_inventory() groups by year and dataset", {
  tmp <- withr::local_tempdir()
  # status data in 2023
  sp_dir <- file.path(tmp, "2023", "woothr")
  dir.create(sp_dir, recursive = TRUE)
  writeLines("data", file.path(sp_dir, "config.json"))
  # trends data in 2022
  trends_dir <- file.path(tmp, "2022", "woothr", "trends")
  dir.create(trends_dir, recursive = TRUE)
  writeLines("data", file.path(trends_dir, "woothr_trends.parquet"))

  out <- capture.output(print(ebirdst_data_inventory(tmp)))
  expect_match(out[1], "2 packages")
  expect_true(any(grepl("2022 Trends Data Products", out)))
  expect_true(any(grepl("2023 Status Data Products", out)))
  # species line should be indented under its group
  status_group <- which(grepl("2023 Status", out))
  trends_group <- which(grepl("2022 Trends", out))
  expect_true(any(grepl("^  Wood Thrush", out[status_group + 1])))
  expect_true(any(grepl("^  Wood Thrush", out[trends_group + 1])))
})

test_that("print.ebirdst_inventory() uses species code when common name is NA", {
  tmp <- withr::local_tempdir()
  sp_dir <- file.path(tmp, "2023", "notaspecies")
  dir.create(sp_dir, recursive = TRUE)
  writeLines("data", file.path(sp_dir, "file.tif"))

  out <- capture.output(print(ebirdst_data_inventory(tmp)))
  # no common name available; code used directly without parentheses
  expect_true(any(grepl("notaspecies:", out)))
  expect_false(any(grepl("NA", out)))
})

test_that("print.ebirdst_inventory() prints empty inventory cleanly", {
  out <- capture.output(print(ebirdst_data_inventory("/nonexistent/path/xyz")))
  expect_match(out[1], "0 species, 0 packages")
  expect_equal(length(out), 1)
})

test_that("print.ebirdst_inventory() returns x invisibly", {
  inv <- ebirdst_data_inventory("/nonexistent/path/xyz")
  result <- withVisible(print(inv))
  expect_false(result$visible)
  expect_identical(result$value, inv)
})


# ebirdst_delete ----

test_that("ebirdst_delete() validates arguments", {
  expect_error(ebirdst_delete(path = 123))
  expect_error(ebirdst_delete(path = c("a", "b")))
  expect_error(ebirdst_delete(force = NA))
  expect_error(ebirdst_delete(force = "yes"))
  expect_error(ebirdst_delete(species = 1L))
  expect_error(ebirdst_delete(year = "2023"))
  expect_error(ebirdst_delete(year = -1L))
})

test_that("ebirdst_delete() messages when no data found", {
  expect_message(
    result <- ebirdst_delete(path = "/nonexistent/path/xyz", force = TRUE),
    "No eBird Status and Trends data found"
  )
  expect_identical(result, character(0))
})

test_that("ebirdst_delete() messages when no matching species found", {
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "2023", "yebsap"), recursive = TRUE)
  file.create(file.path(tmp, "2023", "yebsap", "dummy.tif"))
  expect_message(
    result <- ebirdst_delete(species = "Wood Thrush", path = tmp,
                              force = TRUE),
    "No matching data found"
  )
  expect_identical(result, character(0))
})

test_that("ebirdst_delete() messages when no matching year found", {
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "2023", "yebsap"), recursive = TRUE)
  file.create(file.path(tmp, "2023", "yebsap", "dummy.tif"))
  expect_message(
    result <- ebirdst_delete(year = 2021L, path = tmp, force = TRUE),
    "No matching data found"
  )
  expect_identical(result, character(0))
})

test_that("ebirdst_delete() deletes with force = TRUE", {
  tmp <- withr::local_tempdir()
  sp_dir <- file.path(tmp, "2023", "yebsap")
  dir.create(sp_dir, recursive = TRUE)
  writeLines("data", file.path(sp_dir, "dummy.tif"))

  inv_before <- ebirdst_data_inventory(tmp)
  expect_equal(nrow(inv_before), 1)

  suppressMessages(
    deleted <- ebirdst_delete(path = tmp, force = TRUE)
  )
  expect_type(deleted, "character")
  expect_equal(length(deleted), 1)
  expect_false(dir.exists(sp_dir))
  # year directory should also be cleaned up when empty
  expect_false(dir.exists(file.path(tmp, "2023")))

  inv_after <- ebirdst_data_inventory(tmp)
  expect_equal(nrow(inv_after), 0)
})

test_that("ebirdst_delete() deletes both status and trends from one directory", {
  tmp <- withr::local_tempdir()
  sp_dir <- file.path(tmp, "2023", "woothr")
  dir.create(file.path(sp_dir, "trends"), recursive = TRUE)
  writeLines("data", file.path(sp_dir, "config.json"))
  writeLines("data", file.path(sp_dir, "trends", "woothr_trends.parquet"))

  inv_before <- ebirdst_data_inventory(tmp)
  expect_equal(nrow(inv_before), 2)

  suppressMessages(
    deleted <- ebirdst_delete(path = tmp, force = TRUE)
  )
  # only one directory deleted (the species dir containing both datasets)
  expect_equal(length(deleted), 1)
  expect_false(dir.exists(sp_dir))
})

test_that("ebirdst_delete() filters by species and year", {
  tmp <- withr::local_tempdir()
  # create two species in two years
  for (yr in c("2022", "2023")) {
    for (sp in c("yebsap", "woothr")) {
      d <- file.path(tmp, yr, sp)
      dir.create(d, recursive = TRUE)
      writeLines("data", file.path(d, "dummy.tif"))
    }
  }

  # delete only woothr in 2023
  suppressMessages(
    deleted <- ebirdst_delete(species = "woothr", year = 2023L,
                               path = tmp, force = TRUE)
  )
  expect_equal(length(deleted), 1)
  expect_false(dir.exists(file.path(tmp, "2023", "woothr")))
  # others should remain
  expect_true(dir.exists(file.path(tmp, "2022", "yebsap")))
  expect_true(dir.exists(file.path(tmp, "2022", "woothr")))
  expect_true(dir.exists(file.path(tmp, "2023", "yebsap")))
})

test_that("ebirdst_delete() warns on unrecognized species", {
  tmp <- withr::local_tempdir()
  sp_dir <- file.path(tmp, "2023", "yebsap")
  dir.create(sp_dir, recursive = TRUE)
  writeLines("data", file.path(sp_dir, "dummy.tif"))

  expect_warning(
    suppressMessages(
      ebirdst_delete(species = c("yebsap", "NOTASPECIES"), path = tmp,
                      force = TRUE)
    ),
    "Unrecognized species"
  )
})

test_that("ebirdst_delete() returns invisible character vector", {
  tmp <- withr::local_tempdir()
  sp_dir <- file.path(tmp, "2023", "yebsap")
  dir.create(sp_dir, recursive = TRUE)
  writeLines("data", file.path(sp_dir, "dummy.tif"))

  suppressMessages({
    result <- withVisible(ebirdst_delete(path = tmp, force = TRUE))
  })
  expect_false(result$visible)
  expect_type(result$value, "character")
})
