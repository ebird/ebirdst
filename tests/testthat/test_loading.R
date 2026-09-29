skip_on_cran()

test_that("load_config()", {
  p <- load_config("yebsap-example")
  expect_type(p, "list")
  expect_true(all(c("bins", "bins_seasonal", "srd_pred_year") %in% names(p)))

  expect_error(load_config("XXXX"))
})


test_that("load_config() downloads data on demand", {
  tmp <- withr::local_tempdir()
  p <- suppressMessages(load_config("yebsap-example", path = tmp))
  expect_type(p, "list")
  expect_true("srd_pred_year" %in% names(p))
  cfg <- file.path(
    tmp,
    ebirdst_version()[["status_version_year"]],
    "yebsap-example",
    "config.json"
  )
  expect_true(file.exists(cfg))
})


test_that("load_fac_map_parameters()", {
  p <- load_fac_map_parameters("yebsap-example")
  expect_type(p, "list")
  expect_named(
    p,
    c(
      "custom_projection",
      "fa_extent",
      "res",
      "fa_extent_projected",
      "weekly_bins",
      "weekly_labels",
      "seasonal_bins",
      "seasonal_labels"
    )
  )

  # check components
  # projection
  expect_type(terra::crs(p$custom_projection), "character")
  # extent
  expect_s4_class(p$fa_extent, "SpatExtent")
  # resolution
  expect_true(is.numeric(p$res))
  # projected extent
  expect_s4_class(p$fa_extent_projected, "SpatExtent")
  # bins
  expect_type(p$weekly_bins, "double")
  expect_type(p$seasonal_bins, "double")

  expect_error(load_fac_map_parameters("XXXX"))
})


test_that("list_available_pis()", {
  pis <- list_available_pis("yebsap-example")
  expect_s3_class(pis, "data.frame")
  expect_true(all(c("predictor", "rank_mean", "rank") %in% names(pis)))
  expect_equal(nrow(pis), 10)

  expect_error(list_available_pis("XXXX"))
})


test_that("load_pi()", {
  pi_occ <- load_pi(
    "yebsap-example",
    predictor = "gsw_c2_pland",
    response = "occurrence"
  )
  expect_s4_class(pi_occ, "SpatRaster")
  expect_equal(terra::nlyr(pi_occ), 52)

  pi_count <- load_pi(
    "yebsap-example",
    predictor = "mcd12q1_lccs1_c21_pland",
    response = "count"
  )
  expect_s4_class(pi_count, "SpatRaster")
  expect_equal(terra::nlyr(pi_count), 52)

  expect_error(load_pi("XXXX", predictor = "gsw_c2_pland"))
  expect_error(load_pi("yebsap-example", predictor = "gsw_c2_pland", response = "abundance"))
  expect_error(load_pi("yebsap-example", predictor = "elevation_250m_sd"))
})


test_that("load_ppm()", {
  ppm_pois <- load_ppm("yebsap-example", ppm = "abd_poisson_dev")
  expect_s4_class(ppm_pois, "SpatRaster")
  expect_equal(terra::nlyr(ppm_pois), 52)

  ppm_fs <- load_ppm("yebsap-example", ppm = "binary_f1")
  expect_s4_class(ppm_fs, "SpatRaster")
  expect_equal(terra::nlyr(ppm_fs), 52)

  expect_error(load_ppm("XXXX"))
  expect_error(load_ppm("yebsap-example", ppm = "pr_auc"))
  expect_error(load_ppm("yebsap-example", ppm = "elevation_250m_sd"))
})


test_that("load_pi() downloads data on demand", {
  tmp <- withr::local_tempdir()
  pi <- suppressMessages(load_pi(
    "yebsap-example",
    predictor = "gsw_c2_pland",
    response = "occurrence",
    path = tmp
  ))
  expect_s4_class(pi, "SpatRaster")
  tifs <- list.files(tmp, pattern = "gsw-c2-pland.*\\.tif$", recursive = TRUE)
  expect_length(tifs, 1)
})


test_that("load_ppm() downloads data on demand", {
  tmp <- withr::local_tempdir()
  ppm <- suppressMessages(load_ppm(
    "yebsap-example",
    ppm = "binary_f1",
    path = tmp
  ))
  expect_s4_class(ppm, "SpatRaster")
  tifs <- list.files(tmp, pattern = "binary-f1.*\\.tif$", recursive = TRUE)
  expect_length(tifs, 1)
})


test_that("load_regional_stats() downloads data on demand", {
  tmp <- withr::local_tempdir()
  stats <- suppressMessages(load_regional_stats("yebsap-example", path = tmp))
  expect_s3_class(stats, "data.frame")
  csv <- file.path(
    tmp,
    ebirdst_version()[["status_version_year"]],
    "yebsap-example",
    "regional_stats.csv"
  )
  expect_true(file.exists(csv))
})


test_that("list_available_pis() only downloads the rangewide csv", {
  tmp <- withr::local_tempdir()
  pis <- suppressMessages(list_available_pis("yebsap-example", path = tmp))
  expect_s3_class(pis, "data.frame")
  expect_equal(nrow(pis), 10)

  files <- list.files(tmp, recursive = TRUE)
  expect_true(any(grepl("pi_rangewide.csv$", files)))
  expect_false(any(grepl("\\.tif$", files)))
})


test_that("available_pi_predictors() uses the list of available data", {
  preds <- available_pi_predictors("yebsap-example", path = ebirdst_data_dir())
  expect_type(preds, "character")
  expect_true("gsw_c2_pland" %in% preds)
  # the other tifs stored alongside the pi rasters aren't predictors
  expect_false(any(grepl("folds|day-of-year", preds)))
})


test_that("available_pi_predictors() warns when it falls back to local files", {
  tmp <- withr::local_tempdir()
  pis_dir <- file.path(tmp, status_key("yebsap-example", "pis"))
  dir.create(pis_dir, recursive = TRUE)
  file.create(file.path(
    pis_dir,
    "yebsap-example_pi_occurrence_gsw-c2-pland_27km_2023.tif"
  ))

  local_mocked_bindings(
    list_object_keys = function(...) stop("no internet"),
    .package = "ebirdst"
  )
  expect_warning(
    preds <- available_pi_predictors("yebsap-example", path = tmp),
    "may be incomplete"
  )
  expect_equal(preds, "gsw_c2_pland")
})


test_that("available_pi_predictors() errors when there is nothing to go on", {
  tmp <- withr::local_tempdir()
  local_mocked_bindings(
    list_object_keys = function(...) stop("key has expired"),
    .package = "ebirdst"
  )
  # the underlying problem is surfaced rather than reported as no PI data
  expect_error(
    available_pi_predictors("yebsap-example", path = tmp),
    "key has expired"
  )
})


test_that("ebirdst_regional_stats() loads an existing file", {
  tmp <- withr::local_tempdir()
  version_year <- ebirdst_version()[["status_version_year"]]
  dir.create(file.path(tmp, version_year), recursive = TRUE)
  file <- file.path(
    tmp,
    version_year,
    sprintf("regional-stats_%s.parquet", version_year)
  )
  stats <- dplyr::tibble(species_code = "yebsap", total_pop_percent = 1)
  arrow::write_parquet(stats, file)

  loaded <- ebirdst_regional_stats(path = tmp)
  expect_s3_class(loaded, "tbl_df")
  expect_equal(loaded, stats)
})


test_that("ebirdst_regional_stats() validates arguments", {
  tmp <- withr::local_tempdir()
  expect_error(ebirdst_regional_stats(path = 1))
  expect_error(ebirdst_regional_stats(path = c(tmp, tmp)))
})


test_that("load_data_coverage() requires at least one valid week", {
  # weeks is required; these rasters are ~50 MB each so there is no sensible
  # default and nothing should be downloaded without an explicit request
  expect_error(load_data_coverage())
  expect_error(load_data_coverage(NULL))
  expect_error(load_data_coverage(character(0)))
  expect_error(load_data_coverage(NA_character_))
  expect_error(
    load_data_coverage(c("01-04", "01-05")),
    "weeks are invalid"
  )
})
