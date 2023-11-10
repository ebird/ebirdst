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
  expect_named(p, c("custom_projection", "fa_extent", "res", "fa_extent_sinu",
                    "weekly_bins", "weekly_labels",
                    "seasonal_bins", "seasonal_labels"))

  # check components
  # projection
  expect_is(terra::crs(p$custom_projection), "character")
  # extent
  expect_is(p$fa_extent, "SpatExtent")
  # resolution
  expect_is(p$res, c("numeric", "integer"))
  # sinusoidal extent
  expect_is(p$fa_extent_sinu, "SpatExtent")
  # bins
  expect_is(p$weekly_bins, "numeric")
  expect_is(p$seasonal_bins, "numeric")

  expect_error(load_fac_map_parameters("Yellow Warbler"))
  expect_error(load_fac_map_parameters("XXXX"))
})


test_that("list_available_pis()", {
  pis <- list_available_pis("yebsap-example")
  expect_is(pis, "data.frame")
  expect_true(all(c("predictor", "rangewide_rank", "rank") %in% names(pis)))
  expect_equal(nrow(pis), 5)

  expect_error(list_available_pis("Yellow Warbler"))
  expect_error(list_available_pis("XXXX"))
})


test_that("load_pi()", {
  pi_occ <- load_pi("yebsap-example",
                    predictor = "elevation_250m_median",
                    response = "occurrence")
  expect_is(pi_occ, "SpatRaster")
  expect_equal(terra::nlyr(pi_occ), 52)

  pi_count <- load_pi("yebsap-example",
                      predictor = "mcd12q1_lccs1_c22_pland",
                      response = "count")
  expect_is(pi_count, "SpatRaster")
  expect_equal(terra::nlyr(pi_count), 52)

  expect_error(load_pi("Yellow Warbler"))
  expect_error(load_pi("XXXX"))
  expect_error(load_pi("yebsap-example", response = "abundance"))
  expect_error(load_pi("yebsap-example", predictor = "elevation_250m_sd"))
})


test_that("load_ppm()", {
  ppm_spear <- load_ppm("yebsap-example", ppm = "abd_spearman")
  expect_is(ppm_spear, "SpatRaster")
  expect_equal(terra::nlyr(ppm_spear), 52)

  ppm_auc <- load_ppm("yebsap-example", ppm = "binary_pr_auc")
  expect_is(ppm_auc, "SpatRaster")
  expect_equal(terra::nlyr(ppm_auc), 52)

  expect_error(load_ppm("Yellow Warbler"))
  expect_error(load_ppm("XXXX"))
  expect_error(load_ppm("yebsap-example", ppm = "pr_auc"))
  expect_error(load_ppm("yebsap-example", ppm = "elevation_250m_sd"))
})
