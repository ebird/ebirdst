context("Loading functions")

skip_on_cran()

test_that("load_config()", {
  p <- load_config("yebsap-example")
  expect_is(p, "list")
  expect_true(all(c("bins", "bins_seasonal", "srd_pred_year") %in% names(p)))

  expect_error(load_config("XXXX"))
})


test_that("load_config() downloads data on demand", {
  tmp <- withr::local_tempdir()
  p <- suppressMessages(load_config("yebsap-example", path = tmp))
  expect_is(p, "list")
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
  expect_is(p, "list")
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
  expect_is(terra::crs(p$custom_projection), "character")
  # extent
  expect_is(p$fa_extent, "SpatExtent")
  # resolution
  expect_is(p$res, c("numeric", "integer"))
  # projected extent
  expect_is(p$fa_extent_projected, "SpatExtent")
  # bins
  expect_is(p$weekly_bins, "numeric")
  expect_is(p$seasonal_bins, "numeric")

  expect_error(load_fac_map_parameters("XXXX"))
})


test_that("list_available_pis()", {
  pis <- list_available_pis("yebsap-example")
  expect_is(pis, "data.frame")
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
  expect_is(pi_occ, "SpatRaster")
  expect_equal(terra::nlyr(pi_occ), 52)

  pi_count <- load_pi(
    "yebsap-example",
    predictor = "mcd12q1_lccs1_c21_pland",
    response = "count"
  )
  expect_is(pi_count, "SpatRaster")
  expect_equal(terra::nlyr(pi_count), 52)

  expect_error(load_pi("XXXX", predictor = "gsw_c2_pland"))
  expect_error(load_pi("yebsap-example", predictor = "gsw_c2_pland", response = "abundance"))
  expect_error(load_pi("yebsap-example", predictor = "elevation_250m_sd"))
})


test_that("load_ppm()", {
  ppm_pois <- load_ppm("yebsap-example", ppm = "abd_poisson_dev")
  expect_is(ppm_pois, "SpatRaster")
  expect_equal(terra::nlyr(ppm_pois), 52)

  ppm_fs <- load_ppm("yebsap-example", ppm = "binary_f1")
  expect_is(ppm_fs, "SpatRaster")
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
  expect_is(pi, "SpatRaster")
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
  expect_is(ppm, "SpatRaster")
  tifs <- list.files(tmp, pattern = "binary-f1.*\\.tif$", recursive = TRUE)
  expect_length(tifs, 1)
})


test_that("load_regional_stats() downloads data on demand", {
  tmp <- withr::local_tempdir()
  stats <- suppressMessages(load_regional_stats("yebsap-example", path = tmp))
  expect_is(stats, "data.frame")
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
  expect_is(pis, "data.frame")
  expect_equal(nrow(pis), 10)

  files <- list.files(tmp, recursive = TRUE)
  expect_true(any(grepl("pi_rangewide.csv$", files)))
  expect_false(any(grepl("\\.tif$", files)))
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
  expect_is(loaded, "tbl_df")
  expect_equal(loaded, stats)
})


test_that("ebirdst_regional_stats() validates arguments", {
  tmp <- withr::local_tempdir()
  expect_error(ebirdst_regional_stats(path = 1))
  expect_error(ebirdst_regional_stats(path = c(tmp, tmp)))
})
