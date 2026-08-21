skip_on_cran()

test_that("load_trends()", {
  trends <- load_trends("yebsap-example")
  expect_s3_class(trends, "data.frame")
  cols <- c(
    "species_code",
    "season",
    "start_year",
    "end_year",
    "start_date",
    "end_date",
    "srd_id",
    "longitude",
    "latitude",
    "abd",
    "abd_ppy",
    "abd_ppy_lower",
    "abd_ppy_upper",
    "abd_ppy_nonzero",
    "abd_trend",
    "abd_trend_lower",
    "abd_trend_upper"
  )
  expect_equal(names(trends), cols)
  expect_gt(nrow(trends), 0)

  trends_folds <- load_trends("yebsap-example", fold_estimates = TRUE)
  cols <- c(
    "species_code",
    "season",
    "fold",
    "srd_id",
    "latitude",
    "longitude",
    "abd",
    "abd_ppy"
  )
  expect_equal(names(trends_folds), cols)
  expect_equal(nrow(trends_folds), 100 * nrow(trends))
})

test_that("load_trends() downloads data on demand", {
  tmp <- withr::local_tempdir()
  trends <- suppressMessages(load_trends("yebsap-example", path = tmp))
  expect_s3_class(trends, "data.frame")
  expect_gt(nrow(trends), 0)
})

test_that("load_trends() only downloads the requested parquet file", {
  # a fresh load of the non-fold estimates should not also fetch the fold
  # estimates or the model summary csv
  tmp <- withr::local_tempdir()
  suppressMessages(load_trends("yebsap-example", path = tmp))
  trends_dir <- file.path(
    tmp,
    ebirdst_version()[["trends_version_year"]],
    "yebsap-example",
    "trends"
  )
  expect_equal(list.files(trends_dir), c(
    "yebsap-example_breeding_ebird-trends_2022.parquet"
  ))

  # requesting the fold estimates next should only add that one file
  suppressMessages(load_trends("yebsap-example", path = tmp, fold_estimates = TRUE))
  expect_equal(sort(list.files(trends_dir)), sort(c(
    "yebsap-example_breeding_ebird-trends_2022.parquet",
    "yebsap-example_breeding_ebird-trends_folds_2022.parquet"
  )))
})

test_that("convert_ppy_to_cumulative()", {
  # 0% per year -> 0 cumulative
  expect_equal(convert_ppy_to_cumulative(0, n_years = 10), 0)
  # 100% per year for 1 year -> 100% cumulative
  expect_equal(convert_ppy_to_cumulative(100, n_years = 1), 100)
  # known compounding value
  expect_equal(
    convert_ppy_to_cumulative(10, n_years = 10),
    100 * (1.1^10 - 1),
    tolerance = 1e-6
  )
  # vectorized, sign follows input
  result <- convert_ppy_to_cumulative(c(-50, 0, 50), n_years = 5)
  expect_length(result, 3L)
  expect_true(result[1] < 0)
  expect_equal(result[2], 0)
  expect_true(result[3] > 0)
  # non-integer n_years errors
  expect_error(convert_ppy_to_cumulative(10, n_years = 1.5))
})

test_that("categorize()", {
  expect_equal(
    categorize(c(0.5, 1.5, 2.5), breaks = c(0, 1, 2, 3), labels = c(1, 2, 3)),
    c(1, 2, 3)
  )
  # values on the lowest break are included in the first category rather than
  # being dropped as missing
  expect_equal(categorize(0, breaks = c(0, 1, 2), labels = c(10, 20)), 10)
  # values outside the breaks have no category
  expect_equal(
    categorize(c(-1, 5), breaks = c(0, 1, 2), labels = c(10, 20)),
    c(NA_real_, NA_real_)
  )
  expect_equal(
    categorize(c(0.5, 1.5), breaks = c(0, 1, 2), labels = c("a", "b")),
    c("a", "b")
  )

  expect_error(categorize("a", breaks = c(0, 1), labels = 1))
  expect_error(categorize(1, breaks = c(0, 1), labels = c(1, 2)))
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
  expect_equal(
    as.character(sf::st_geometry_type(pts, by_geometry = FALSE)),
    "POINT"
  )

  # circles output (equal-area CRS required for sensible geometry)
  circles <- vectorize_trends(trends, output = "circles", crs = "+proj=eqearth")
  expect_s3_class(circles, "sf")
  expect_equal(nrow(circles), nrow(trends))
})

test_that("vectorize_trends() drops zero-abundance cells from circles output", {
  # spec1 mixes zero and nonzero abundance, spec2 has a single distinct
  # nonzero abundance value (one quantile bin), spec3 is entirely zero
  trends <- data.frame(
    species_code = c(rep("spec1", 4), rep("spec2", 3), rep("spec3", 2)),
    latitude = c(10, 11, 12, 13, 20, 21, 22, 30, 31),
    longitude = c(-80, -81, -82, -83, -90, -91, -92, -100, -101),
    abd = c(0, 1, 2, 0, 5, 5, 5, 0, 0)
  )

  circles <- vectorize_trends(trends, output = "circles", crs = "+proj=eqearth")
  expect_equal(nrow(circles), 5L)
  expect_true(all(circles$abd > 0))
  expect_false("spec3" %in% circles$species_code)

  # every remaining spec2 cell shares one abundance value, so all get the
  # maximum circle radius
  max_radius <- 0.99 * min(terra::res(trends_raster_template())) / 2
  spec2_radii <- circles$radii[circles$species_code == "spec2"]
  expect_equal(spec2_radii, rep(max_radius, 3))
})
