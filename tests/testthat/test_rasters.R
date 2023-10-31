context("Loading raster data")

skip_on_cran()

test_that("load_raster()", {
  # weekly
  abd <- load_raster("yebsap-example",
                     product = "abundance",
                     resolution = "27km")
  expect_is(abd, "SpatRaster")
  expect_equal(terra::nlyr(abd), 52)

  # check labellling
  expect_match(names(abd), "^[0-9]{4}-[0-9]{2}-[0-9]{2}")
  expect_is(as.Date(names(abd)), "Date")

  # seasonal
  abd <- load_raster("yebsap-example",
                     product = "abundance",
                     period = "seasonal",
                     resolution = "27km")
  expect_is(abd, "SpatRaster")
  expect_equal(terra::nlyr(abd), 4)
  expect_named(abd, c("breeding", "nonbreeding",
                      "prebreeding_migration", "postbreeding_migration"))
})

test_that("load_raster() error", {
  # weekly
  expect_error(load_raster("Yellow Warbler"))
  expect_error(load_raster("yebsap-example", product = "abndnce"))
  expect_error(load_raster("yebsap-example", resolution = "57km"))
  expect_error(load_raster("yebsap-example", path = "/bad/path/"))
})
