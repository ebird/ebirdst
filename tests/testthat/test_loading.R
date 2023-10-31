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
