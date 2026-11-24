skip_on_cran()

# shared test region: bounding box of the example raster data
test_region <- local({
  rast <- suppressMessages(
    load_raster("yebsap-example", "abundance", resolution = "27km")
  )
  sf::st_as_sf(terra::as.polygons(terra::ext(rast), crs = terra::crs(rast)))
})

# helper: write a minimal regional-stats parquet for region_code tests
make_region_stats_parquet <- function(path, region_code = "US-TEST") {
  version_year <- ebirdst_version()[["status_version_year"]]
  dir_path <- file.path(path, version_year)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  file <- file.path(
    dir_path,
    sprintf("regional-stats_%s.parquet", version_year)
  )
  stats <- data.frame(
    species_code = rep("yebsap-example", 2L),
    region_type = rep("state", 2L),
    region_code = rep(region_code, 2L),
    region_name = rep("Test Region", 2L),
    continent_code = rep("NA", 2L),
    continent_name = rep("North America", 2L),
    season = c("breeding", "nonbreeding"),
    max_week_percent_pop = c(0.5, 0.05),
    total_pop_percent = c(0.4, 0.04),
    stringsAsFactors = FALSE
  )
  arrow::write_parquet(stats, file)
  return(invisible(file))
}


# polygon mode ----

test_that("stewardship_connections() errors on unknown species", {
  expect_error(
    stewardship_connections("NOTASPECIES", test_region, resolution = "27km"),
    "not modeled"
  )
})

test_that("stewardship_connections() errors on residents and failed-review species", {
  # ostric2 is a resident; tunbeg1 has breeding_quality == 0
  expect_error(
    stewardship_connections(
      c("ostric2", "tunbeg1"),
      test_region,
      resolution = "27km"
    ),
    "ostric2"
  )
  expect_error(
    stewardship_connections(
      c("ostric2", "tunbeg1"),
      test_region,
      resolution = "27km"
    ),
    "tunbeg1"
  )
})

test_that("stewardship_connections() errors on non-spatial origin_region", {
  expect_error(
    stewardship_connections("yebsap-example", list(x = 1)),
    class = "simpleError"
  )
})

test_that("stewardship_connections() errors on bad filename", {
  expect_error(
    stewardship_connections(
      "yebsap-example",
      test_region,
      filename = c("a.tif", "b.tif")
    ),
    class = "simpleError"
  )
})

test_that("stewardship_connections() errors when species missing for origin_region", {
  expect_error(
    stewardship_connections(origin_region = test_region, resolution = "27km"),
    "species.*required"
  )
})

test_that("stewardship_connections() integration smoke test with yebsap-example", {
  result <- suppressMessages(
    stewardship_connections(
      "yebsap-example",
      origin_region = test_region,
      destination_mask = test_region,
      resolution = "27km"
    )
  )
  expect_type(result, "list")
  expect_named(result, c("connection", "weights"))
  expect_s4_class(result$connection, "SpatRaster")
  expect_equal(terra::nlyr(result$connection), 1L)
  expect_s3_class(result$weights, "data.frame")
  expect_named(
    result$weights,
    c(
      "species_code",
      "common_name",
      "origin_season",
      "weight",
      "rel_pct_decrease"
    )
  )
  expect_equal(nrow(result$weights), 1L)
  expect_true(result$weights$weight >= 0 && result$weights$weight <= 1)
  expect_true("rel_pct_decrease" %in% names(result$weights))
})

test_that("stewardship_connections() errors when no species retained", {
  # min_origin_threshold above any plausible weight → no species retained
  expect_error(
    suppressMessages(
      stewardship_connections(
        "yebsap-example",
        origin_region = test_region,
        min_origin_threshold = 1.0,
        destination_mask = test_region,
        resolution = "27km"
      )
    ),
    "retained after applying"
  )
})

test_that("stewardship_connections() min_migration_threshold conditional default", {
  # default resolves to 0.85 without migration seasons, 0.5 with
  # we can confirm the function runs (or errors sensibly) in both cases
  # run with include_migration_seasons = TRUE to exercise the 0.5 default
  result_mig <- suppressMessages(
    stewardship_connections(
      "yebsap-example",
      origin_region = test_region,
      include_migration_seasons = TRUE,
      destination_mask = test_region,
      resolution = "27km"
    )
  )
  expect_s4_class(result_mig$connection, "SpatRaster")
})

test_that("stewardship_connections() accepts SpatVector origin_region", {
  sv_region <- terra::vect(test_region)
  result <- suppressMessages(
    stewardship_connections(
      "yebsap-example",
      origin_region = sv_region,
      destination_mask = test_region,
      resolution = "27km"
    )
  )
  expect_s4_class(result$connection, "SpatRaster")
})


# region_code mode ----

test_that("stewardship_connections() region_code mode returns correct structure", {
  tmp <- withr::local_tempdir()
  make_region_stats_parquet(tmp)
  result <- suppressMessages(
    stewardship_connections(
      region_code = "US-TEST",
      destination_mask = test_region,
      resolution = "27km",
      path = tmp
    )
  )
  expect_type(result, "list")
  expect_named(result, c("connection", "weights"))
  expect_s4_class(result$connection, "SpatRaster")
  expect_equal(terra::nlyr(result$connection), 1L)
  expect_s3_class(result$weights, "data.frame")
  expect_named(
    result$weights,
    c(
      "species_code",
      "common_name",
      "origin_season",
      "weight",
      "rel_pct_decrease"
    )
  )
  expect_equal(nrow(result$weights), 1L)
  # weight must equal total_pop_percent from the stub (0.4 for breeding)
  expect_equal(result$weights$weight, 0.4)
  expect_equal(result$weights$species_code, "yebsap-example")
})

test_that("stewardship_connections() region_code errors when code absent", {
  tmp <- withr::local_tempdir()
  make_region_stats_parquet(tmp)
  expect_error(
    stewardship_connections(
      region_code = "NOTAREGION",
      destination_mask = test_region,
      resolution = "27km",
      path = tmp
    ),
    "No data found for region_code"
  )
})

test_that("stewardship_connections() region_code + species subset narrows list", {
  tmp <- withr::local_tempdir()
  # add a second species to the stub so subsetting is testable
  version_year <- ebirdst_version()[["status_version_year"]]
  file <- file.path(
    tmp,
    version_year,
    sprintf("regional-stats_%s.parquet", version_year)
  )
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  stats <- data.frame(
    species_code = c(
      rep("yebsap-example", 2L),
      rep("woothr", 2L)
    ),
    region_type = rep("state", 4L),
    region_code = rep("US-TEST", 4L),
    region_name = rep("Test Region", 4L),
    continent_code = rep("NA", 4L),
    continent_name = rep("North America", 4L),
    season = rep(c("breeding", "nonbreeding"), 2L),
    max_week_percent_pop = c(0.5, 0.05, 0.3, 0.02),
    total_pop_percent = c(0.4, 0.04, 0.25, 0.015),
    stringsAsFactors = FALSE
  )
  arrow::write_parquet(stats, file)

  result <- suppressMessages(
    stewardship_connections(
      species = "yebsap-example",
      region_code = "US-TEST",
      destination_mask = test_region,
      resolution = "27km",
      path = tmp
    )
  )
  expect_equal(result$weights$species_code, "yebsap-example")
  expect_equal(nrow(result$weights), 1L)
})


# species_weights mode ----

test_that("stewardship_connections() species_weights mode returns correct structure", {
  wt <- data.frame(
    species_code = "yebsap-example",
    weight = 0.3,
    stringsAsFactors = FALSE
  )
  result <- suppressMessages(
    stewardship_connections(
      species_weights = wt,
      destination_mask = test_region,
      resolution = "27km"
    )
  )
  expect_type(result, "list")
  expect_named(result, c("connection", "weights"))
  expect_s4_class(result$connection, "SpatRaster")
  expect_equal(terra::nlyr(result$connection), 1L)
  expect_equal(nrow(result$weights), 1L)
  expect_equal(result$weights$weight, 0.3)
  expect_true(is.na(result$weights$rel_pct_decrease))
  expect_true(is.na(result$weights$origin_season))
})

test_that("stewardship_connections() species_weights errors on missing column", {
  expect_error(
    stewardship_connections(
      species_weights = data.frame(species_code = "yebsap-example"),
      destination_mask = test_region
    ),
    "weight"
  )
  expect_error(
    stewardship_connections(
      species_weights = data.frame(weight = 0.3),
      destination_mask = test_region
    ),
    "species_code"
  )
})

test_that("stewardship_connections() species_weights errors on negative weight", {
  expect_error(
    stewardship_connections(
      species_weights = data.frame(
        species_code = "yebsap-example",
        weight = -0.1
      ),
      destination_mask = test_region
    ),
    "non-negative"
  )
})

test_that("stewardship_connections() species_weights errors on duplicate species", {
  expect_error(
    stewardship_connections(
      species_weights = data.frame(
        species_code = c("yebsap-example", "yebsap-example"),
        weight = c(0.3, 0.2)
      ),
      destination_mask = test_region
    ),
    "duplicate"
  )
})


# mode validation ----

test_that("stewardship_connections() errors when no mode supplied", {
  expect_error(
    stewardship_connections(destination_mask = test_region),
    "Exactly one of"
  )
})

test_that("stewardship_connections() errors when multiple modes supplied", {
  wt <- data.frame(species_code = "yebsap-example", weight = 0.3)
  expect_error(
    stewardship_connections(
      origin_region = test_region,
      species_weights = wt
    ),
    "Exactly one of"
  )
})


# importance_breaks ----

test_that("importance_breaks() default returns numeric breaks", {
  r <- terra::rast(nrows = 10, ncols = 10, vals = runif(100, 0, 0.5))
  brks <- importance_breaks(r)
  expect_type(brks, "double")
  expect_length(brks, 12L) # n_bins + 2
  expect_equal(brks[1], 0)
  expect_true(all(diff(brks) >= 0)) # monotonic
})

test_that("importance_breaks() return_raster = TRUE gives classified SpatRaster", {
  r <- terra::rast(nrows = 10, ncols = 10, vals = runif(100, 0, 0.5))
  r_class <- importance_breaks(r, return_raster = TRUE)
  expect_s4_class(r_class, "SpatRaster")
  expect_equal(terra::nlyr(r_class), 1L)
  vals <- terra::values(r_class, na.rm = TRUE)
  expect_true(all(vals %in% seq_len(10)))
})

test_that("importance_breaks() validates arguments", {
  r <- terra::rast(nrows = 10, ncols = 10, vals = runif(100))
  expect_error(
    importance_breaks(r, lower_quantile = 1.5),
    class = "simpleError"
  )
  expect_error(importance_breaks(r, n_bins = 0), class = "simpleError")
  expect_error(
    importance_breaks(r, return_raster = "yes"),
    class = "simpleError"
  )
  # multi-layer raster
  r2 <- c(r, r)
  expect_error(importance_breaks(r2), class = "simpleError")
})
