context("Loading trends data")

skip_on_cran()

test_that("load_trends()", {
  trends <- load_trends("yebsap-example")
  expect_s3_class(trends, "data.frame")
  cols <- c("species_code", "season", "start_year", "end_year", "start_date",
            "end_date", "srd_id", "longitude", "latitude", "abd", "abd_ppy",
            "abd_ppy_lower", "abd_ppy_upper", "abd_ppy_nonzero", "abd_trend",
            "abd_trend_lower", "abd_trend_upper")
  expect_equal(names(trends), cols)
  expect_gt(nrow(trends), 0)

  trends_folds <- load_trends("yebsap-example", fold_estimates = TRUE)
  cols <- c("species_code", "season", "fold", "srd_id", "latitude", "longitude",
            "abd", "abd_ppy")
  expect_equal(names(trends_folds), cols)
  expect_equal(nrow(trends_folds), 100 * nrow(trends))
})

test_that("convert_ppy_to_cumulative()", {
  # 0% per year -> 0 cumulative
  expect_equal(convert_ppy_to_cumulative(0, n_years = 10), 0)
  # 100% per year for 1 year -> 100% cumulative
  expect_equal(convert_ppy_to_cumulative(100, n_years = 1), 100)
  # known compounding value
  expect_equal(convert_ppy_to_cumulative(10, n_years = 10),
               100 * (1.1^10 - 1), tolerance = 1e-6)
  # vectorized, sign follows input
  result <- convert_ppy_to_cumulative(c(-50, 0, 50), n_years = 5)
  expect_length(result, 3L)
  expect_true(result[1] < 0)
  expect_equal(result[2], 0)
  expect_true(result[3] > 0)
  # non-integer n_years errors
  expect_error(convert_ppy_to_cumulative(10, n_years = 1.5))
})

test_that("rasterize_trends()", {
  trends <- load_trends("yebsap-example")

  # single layer
  r <- rasterize_trends(trends, layers = "abd_ppy")
  expect_s4_class(r, "SpatRaster")
  expect_equal(terra::nlyr(r), 1L)
  expect_equal(names(r), "abd_ppy")

  # multiple layers
  r_multi <- rasterize_trends(trends, layers = c("abd_ppy", "abd_ppy_lower"))
  expect_equal(terra::nlyr(r_multi), 2L)

  # trim = FALSE returns a larger (global) raster
  r_notrim <- rasterize_trends(trends, layers = "abd_ppy", trim = FALSE)
  expect_gt(terra::ncell(r_notrim), terra::ncell(r))

  # fold-level estimates cannot be rasterized
  trends_folds <- load_trends("yebsap-example", fold_estimates = TRUE)
  expect_error(rasterize_trends(trends_folds))

  # invalid layer name errors
  expect_error(rasterize_trends(trends, layers = "invalid_col"))
})

test_that("vectorize_trends()", {
  trends <- load_trends("yebsap-example")

  # points output
  pts <- vectorize_trends(trends, output = "points")
  expect_s3_class(pts, "sf")
  expect_equal(nrow(pts), nrow(trends))
  expect_equal(as.character(sf::st_geometry_type(pts, by_geometry = FALSE)),
               "POINT")

  # circles output (equal-area CRS required for sensible geometry)
  circles <- vectorize_trends(trends, output = "circles",
                              crs = "+proj=eqearth")
  expect_s3_class(circles, "sf")
  expect_equal(nrow(circles), nrow(trends))
})
