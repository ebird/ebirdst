skip_on_cran()

# shared test region: bounding box of the example raster data
test_region <- local({
  rast <- suppressMessages(
    load_raster("yebsap-example", "abundance", resolution = "27km")
  )
  sf::st_as_sf(terra::as.polygons(terra::ext(rast), crs = terra::crs(rast)))
})

test_that("migration_chronology() basic output", {
  chron <- suppressMessages(
    migration_chronology("yebsap-example", test_region, resolution = "27km")
  )
  expect_s3_class(chron, "data.frame")
  expect_named(chron, c("species_code", "season", "week",
                        "proportion_population"))
  expect_equal(nrow(chron), 52L)
  expect_equal(unique(chron$species_code), "yebsap-example")
  expect_s3_class(chron$week, "Date")
  expect_true(is.numeric(chron$proportion_population))
  # seasons may include NA for weeks outside any season
  season_vals <- unique(chron$season[!is.na(chron$season)])
  expect_true(all(season_vals %in% c(
    "breeding", "nonbreeding",
    "prebreeding_migration", "postbreeding_migration",
    "resident"
  )))
})

test_that("migration_chronology() single polygon has no .region_number", {
  chron <- suppressMessages(
    migration_chronology("yebsap-example", test_region, resolution = "27km")
  )
  expect_false(".region_number" %in% names(chron))
})

test_that("migration_chronology() multi-polygon adds .region_number", {
  two_regions <- sf::st_make_grid(test_region, n = c(1, 2))
  two_regions <- sf::st_as_sf(two_regions)
  chron2 <- suppressMessages(
    migration_chronology("yebsap-example", two_regions, resolution = "27km")
  )
  expect_true(".region_number" %in% names(chron2))
  expect_equal(nrow(chron2), 104L)  # 52 weeks x 2 regions
  expect_setequal(unique(chron2$.region_number), c(1L, 2L))
})

test_that("migration_chronology() bind_regions attaches region columns", {
  region_named <- test_region
  region_named$myvar <- "test_value"
  chron <- suppressMessages(
    migration_chronology(
      "yebsap-example", region_named,
      resolution = "27km", bind_regions = TRUE
    )
  )
  expect_true("myvar" %in% names(chron))
  expect_equal(unique(chron$myvar), "test_value")
})

test_that("seasonal_proportion_population() basic output", {
  seas <- suppressMessages(
    seasonal_proportion_population(
      "yebsap-example", test_region,
      resolution = "27km"
    )
  )
  expect_s3_class(seas, "data.frame")
  expect_named(seas, c("species_code", "season", "proportion_population"))
  expect_true(nrow(seas) >= 1L)
  expect_equal(unique(seas$species_code), "yebsap-example")
  expect_true(is.numeric(seas$proportion_population))
  expect_true(all(seas$season %in% c(
    "breeding", "nonbreeding",
    "prebreeding_migration", "postbreeding_migration",
    "resident"
  )))
  # single polygon: no .region_number
  expect_false(".region_number" %in% names(seas))
})

test_that("seasonal_proportion_population() multi-polygon adds .region_number", {
  two_regions <- sf::st_as_sf(sf::st_make_grid(test_region, n = c(1, 2)))
  seas2 <- suppressMessages(
    seasonal_proportion_population(
      "yebsap-example", two_regions,
      resolution = "27km"
    )
  )
  expect_true(".region_number" %in% names(seas2))
  expect_setequal(unique(seas2$.region_number), c(1L, 2L))
})

test_that("migration_chronology() errors on unknown species", {
  expect_error(
    migration_chronology("NOTASPECIES", test_region, resolution = "27km"),
    "not found"
  )
})

test_that("migration_chronology() errors on non-spatial region input", {
  expect_error(
    migration_chronology("yebsap-example", list(x = 1)),
    class = "simpleError"
  )
})

test_that("migration_chronology() accepts SpatVector region", {
  sv_region <- terra::vect(test_region)
  chron <- suppressMessages(
    migration_chronology("yebsap-example", sv_region, resolution = "27km")
  )
  expect_s3_class(chron, "data.frame")
  expect_equal(nrow(chron), 52L)
})
