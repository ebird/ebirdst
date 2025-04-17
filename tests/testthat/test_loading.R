context("Loading functions")

skip_on_cran()

test_that("load_config()", {
  p <- load_config("yebsap-example")
  expect_is(p, "list")
  expect_true(all(c("bins", "bins_seasonal", "srd_pred_year") %in% names(p)))

  expect_error(load_config("Yellow Warbler"))
  expect_error(load_config("XXXX"))
})


test_that("load_fac_map_parameters()", {
  p <- load_fac_map_parameters("yebsap-example")
  expect_is(p, "list")
  expect_named(p, c("custom_projection", "fa_extent",
                    "res", "fa_extent_projected",
                    "weekly_bins", "weekly_labels",
                    "seasonal_bins", "seasonal_labels"))

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

  expect_error(load_fac_map_parameters("Yellow Warbler"))
  expect_error(load_fac_map_parameters("XXXX"))
})


test_that("list_available_pis()", {
  pis <- list_available_pis("yebsap-example")
  expect_is(pis, "data.frame")
  expect_true(all(c("predictor", "rank_mean", "rank") %in% names(pis)))
  expect_equal(nrow(pis), 10)

  expect_error(list_available_pis("Yellow Warbler"))
  expect_error(list_available_pis("XXXX"))
})


test_that("load_pi()", {
  pi_occ <- load_pi("yebsap-example",
                    predictor = "gsw_c2_pland",
                    response = "occurrence")
  expect_is(pi_occ, "SpatRaster")
  expect_equal(terra::nlyr(pi_occ), 52)

  pi_count <- load_pi("yebsap-example",
                      predictor = "mcd12q1_lccs1_c21_pland",
                      response = "count")
  expect_is(pi_count, "SpatRaster")
  expect_equal(terra::nlyr(pi_count), 52)

  expect_error(load_pi("Yellow Warbler"))
  expect_error(load_pi("XXXX"))
  expect_error(load_pi("yebsap-example", response = "abundance"))
  expect_error(load_pi("yebsap-example", predictor = "elevation_250m_sd"))
})


test_that("load_ppm()", {
  ppm_pois <- load_ppm("yebsap-example", ppm = "abd_poisson_dev")
  expect_is(ppm_pois, "SpatRaster")
  expect_equal(terra::nlyr(ppm_pois), 52)

  ppm_fs <- load_ppm("yebsap-example", ppm = "binary_f1")
  expect_is(ppm_fs, "SpatRaster")
  expect_equal(terra::nlyr(ppm_fs), 52)

  expect_error(load_ppm("Yellow Warbler"))
  expect_error(load_ppm("XXXX"))
  expect_error(load_ppm("yebsap-example", ppm = "pr_auc"))
  expect_error(load_ppm("yebsap-example", ppm = "elevation_250m_sd"))
})
