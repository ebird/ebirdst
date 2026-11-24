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
  expect_named(
    chron,
    c("species_code", "season", "week", "proportion_population")
  )
  expect_equal(nrow(chron), 52L)
  expect_equal(unique(chron$species_code), "yebsap-example")
  expect_s3_class(chron$week, "Date")
  expect_true(is.numeric(chron$proportion_population))
  # seasons may include NA for weeks outside any season
  season_vals <- unique(chron$season[!is.na(chron$season)])
  expect_true(all(
    season_vals %in%
      c(
        "breeding",
        "nonbreeding",
        "prebreeding_migration",
        "postbreeding_migration",
        "resident"
      )
  ))
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
  expect_equal(nrow(chron2), 104L) # 52 weeks x 2 regions
  expect_setequal(unique(chron2$.region_number), c(1L, 2L))
})

test_that("migration_chronology() bind_regions attaches region columns", {
  region_named <- test_region
  region_named$myvar <- "test_value"
  chron <- suppressMessages(
    migration_chronology(
      "yebsap-example",
      region_named,
      resolution = "27km",
      bind_regions = TRUE
    )
  )
  expect_true("myvar" %in% names(chron))
  expect_equal(unique(chron$myvar), "test_value")
})

test_that("seasonal_proportion_population() basic output", {
  seas <- suppressMessages(
    seasonal_proportion_population(
      "yebsap-example",
      test_region,
      resolution = "27km"
    )
  )
  expect_s3_class(seas, "data.frame")
  expect_named(seas, c("species_code", "season", "proportion_population"))
  expect_true(nrow(seas) >= 1L)
  expect_equal(unique(seas$species_code), "yebsap-example")
  expect_true(is.numeric(seas$proportion_population))
  expect_true(all(
    seas$season %in%
      c(
        "breeding",
        "nonbreeding",
        "prebreeding_migration",
        "postbreeding_migration",
        "resident"
      )
  ))
  # single polygon: no .region_number
  expect_false(".region_number" %in% names(seas))
})

test_that("seasonal_proportion_population() multi-polygon adds .region_number", {
  two_regions <- sf::st_as_sf(sf::st_make_grid(test_region, n = c(1, 2)))
  seas2 <- suppressMessages(
    seasonal_proportion_population(
      "yebsap-example",
      two_regions,
      resolution = "27km"
    )
  )
  expect_true(".region_number" %in% names(seas2))
  expect_setequal(unique(seas2$.region_number), c(1L, 2L))
})

test_that("migration_chronology() errors on unknown species", {
  expect_error(
    migration_chronology("NOTASPECIES", test_region, resolution = "27km"),
    "not modeled"
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

test_that("migration_chronology() reference_region computes relative proportions", {
  rast <- suppressMessages(
    load_raster(
      "yebsap-example",
      product = "proportion-population",
      period = "seasonal",
      resolution = "27km"
    )
  )
  trimmed_region <- sf::st_as_sf(
    terra::as.polygons(terra::ext(terra::trim(rast)), crs = terra::crs(rast))
  )
  two_regions <- sf::st_as_sf(sf::st_make_grid(trimmed_region, n = c(2, 1)))

  chron_ref <- suppressMessages(
    migration_chronology(
      "yebsap-example",
      two_regions,
      reference_region = trimmed_region,
      resolution = "27km"
    )
  )
  chron_reg <- suppressMessages(
    migration_chronology("yebsap-example", two_regions, resolution = "27km")
  )
  chron_tot <- suppressMessages(
    migration_chronology("yebsap-example", trimmed_region, resolution = "27km")
  )

  denom <- rep(chron_tot$proportion_population, times = 2L)
  expected <- ifelse(
    denom == 0 & chron_reg$proportion_population == 0,
    0,
    ifelse(denom == 0, NA_real_, chron_reg$proportion_population / denom)
  )
  expect_named(
    chron_ref,
    c(
      "species_code",
      "season",
      "week",
      "proportion_population",
      "proportion_population_reference",
      ".region_number"
    )
  )
  expect_equal(chron_ref$proportion_population, expected)
  expect_equal(chron_ref$proportion_population_reference, denom)

  # partitioning reference_region into two regions sums to 1 for non-zero weeks
  r1 <- chron_ref$proportion_population[chron_ref$.region_number == 1L]
  r2 <- chron_ref$proportion_population[chron_ref$.region_number == 2L]
  pos <- chron_tot$proportion_population > 0
  expect_true(any(pos))
  expect_true(any(r1[pos] > 0 & r1[pos] < 1))
  expect_equal(r1[pos] + r2[pos], rep(1, sum(pos)), tolerance = 1e-6)
})

test_that("seasonal_proportion_population() reference_region computes relative proportions", {
  rast <- suppressMessages(
    load_raster(
      "yebsap-example",
      product = "proportion-population",
      period = "seasonal",
      resolution = "27km"
    )
  )
  trimmed_region <- sf::st_as_sf(
    terra::as.polygons(terra::ext(terra::trim(rast)), crs = terra::crs(rast))
  )
  two_regions <- sf::st_as_sf(sf::st_make_grid(trimmed_region, n = c(2, 1)))

  seas_ref <- suppressMessages(
    seasonal_proportion_population(
      "yebsap-example",
      two_regions,
      reference_region = trimmed_region,
      resolution = "27km"
    )
  )
  seas_reg <- suppressMessages(
    seasonal_proportion_population(
      "yebsap-example",
      two_regions,
      resolution = "27km"
    )
  )
  seas_tot <- suppressMessages(
    seasonal_proportion_population(
      "yebsap-example",
      trimmed_region,
      resolution = "27km"
    )
  )

  denom <- rep(seas_tot$proportion_population, times = 2L)
  expected <- ifelse(
    denom == 0 & seas_reg$proportion_population == 0,
    0,
    ifelse(denom == 0, NA_real_, seas_reg$proportion_population / denom)
  )
  expect_named(
    seas_ref,
    c(
      "species_code",
      "season",
      "proportion_population",
      "proportion_population_reference",
      ".region_number"
    )
  )
  expect_equal(seas_ref$proportion_population, expected)
  expect_equal(seas_ref$proportion_population_reference, denom)

  # two halves sum to 1 for all positive seasons
  r1 <- seas_ref$proportion_population[seas_ref$.region_number == 1L]
  r2 <- seas_ref$proportion_population[seas_ref$.region_number == 2L]
  pos <- seas_tot$proportion_population > 0
  expect_true(any(pos))
  expect_true(all(r1[pos] > 0 & r1[pos] < 1))
  expect_equal(r1[pos] + r2[pos], rep(1, sum(pos)), tolerance = 1e-6)

  # region equal to reference_region gives 1 for all positive seasons
  seas_self <- suppressMessages(
    seasonal_proportion_population(
      "yebsap-example",
      trimmed_region,
      reference_region = trimmed_region,
      resolution = "27km"
    )
  )
  expect_named(
    seas_self,
    c(
      "species_code",
      "season",
      "proportion_population",
      "proportion_population_reference"
    )
  )
  expect_equal(seas_self$proportion_population[pos], rep(1, sum(pos)))
  expect_equal(
    seas_self$proportion_population_reference,
    seas_tot$proportion_population
  )
})

test_that("reference_region unions multiple and overlapping features and accepts SpatVector", {
  rast <- suppressMessages(
    load_raster(
      "yebsap-example",
      product = "proportion-population",
      period = "seasonal",
      resolution = "27km"
    )
  )
  trimmed_region <- sf::st_as_sf(
    terra::as.polygons(terra::ext(terra::trim(rast)), crs = terra::crs(rast))
  )
  two_regions <- sf::st_as_sf(sf::st_make_grid(trimmed_region, n = c(2, 1)))
  sub_region <- two_regions[1, ]

  # overlapping multi-feature reference_region covering the exact same area as trimmed_region
  overlapping_ref <- rbind(
    sf::st_sf(geometry = sf::st_geometry(trimmed_region)),
    sf::st_sf(geometry = sf::st_geometry(two_regions[1, ]))
  )

  seas_single <- suppressMessages(
    seasonal_proportion_population(
      "yebsap-example",
      sub_region,
      reference_region = trimmed_region,
      resolution = "27km"
    )
  )
  seas_multi <- suppressMessages(
    seasonal_proportion_population(
      "yebsap-example",
      sub_region,
      reference_region = two_regions,
      resolution = "27km"
    )
  )
  seas_overlap <- suppressMessages(
    seasonal_proportion_population(
      "yebsap-example",
      sub_region,
      reference_region = overlapping_ref,
      resolution = "27km"
    )
  )
  seas_sv <- suppressMessages(
    seasonal_proportion_population(
      "yebsap-example",
      sub_region,
      reference_region = terra::vect(overlapping_ref),
      resolution = "27km"
    )
  )

  expect_true(all(
    seas_single$proportion_population > 0 &
      seas_single$proportion_population < 1
  ))
  expect_equal(
    seas_multi$proportion_population,
    seas_single$proportion_population,
    tolerance = 1e-6
  )
  expect_equal(
    seas_overlap$proportion_population,
    seas_single$proportion_population,
    tolerance = 1e-6
  )
  expect_equal(
    seas_sv$proportion_population,
    seas_single$proportion_population,
    tolerance = 1e-6
  )
})

test_that("reference_region handles zero denominator correctly", {
  r <- suppressMessages(
    load_raster(
      "yebsap-example",
      product = "proportion-population",
      period = "seasonal",
      resolution = "27km"
    )
  )
  vals <- terra::values(r)
  # find a cell with 0 in at least one season (prefer 0 across all seasons)
  all_zero <- which(apply(vals, 1, function(x) all(!is.na(x) & x == 0)))
  any_zero <- which(apply(vals, 1, function(x) any(!is.na(x) & x == 0)))
  zero_cell <- if (length(all_zero) > 0) all_zero[1] else any_zero[1]
  skip_if(length(zero_cell) == 0 || is.na(zero_cell))

  xy <- terra::xyFromCell(r, zero_cell)
  res_xy <- terra::res(r)
  zero_poly <- sf::st_as_sf(terra::as.polygons(
    terra::ext(
      xy[1] - res_xy[1] / 4,
      xy[1] + res_xy[1] / 4,
      xy[2] - res_xy[2] / 4,
      xy[2] + res_xy[2] / 4
    ),
    crs = terra::crs(r)
  ))

  # both region and reference_region have 0 -> 0 (not NaN)
  seas_zero_both <- suppressMessages(
    seasonal_proportion_population(
      "yebsap-example",
      zero_poly,
      reference_region = zero_poly,
      resolution = "27km"
    )
  )
  zero_seasons <- which(vals[zero_cell, ] == 0)
  expect_equal(
    seas_zero_both$proportion_population[zero_seasons],
    rep(0, length(zero_seasons))
  )
  expect_equal(
    seas_zero_both$proportion_population_reference[zero_seasons],
    rep(0, length(zero_seasons))
  )
  expect_false(any(is.nan(seas_zero_both$proportion_population[zero_seasons])))

  # region has > 0 while reference_region has 0 -> NA_real_ (not Inf)
  seas_zero_ref <- suppressMessages(
    seasonal_proportion_population(
      "yebsap-example",
      test_region,
      reference_region = zero_poly,
      resolution = "27km"
    )
  )
  expect_true(all(is.na(seas_zero_ref$proportion_population[zero_seasons])))
  expect_equal(
    seas_zero_ref$proportion_population_reference[zero_seasons],
    rep(0, length(zero_seasons))
  )
  expect_false(any(is.infinite(seas_zero_ref$proportion_population[
    zero_seasons
  ])))
})

test_that("reference_region validates input type and non-emptiness", {
  expect_error(
    migration_chronology(
      "yebsap-example",
      test_region,
      reference_region = list(x = 1)
    ),
    class = "simpleError"
  )
  expect_error(
    seasonal_proportion_population(
      "yebsap-example",
      test_region,
      reference_region = list(x = 1)
    ),
    class = "simpleError"
  )
  empty_sf <- test_region[0, ]
  expect_error(
    migration_chronology(
      "yebsap-example",
      test_region,
      reference_region = empty_sf
    ),
    class = "simpleError"
  )
  expect_error(
    seasonal_proportion_population(
      "yebsap-example",
      test_region,
      reference_region = empty_sf
    ),
    class = "simpleError"
  )

  empty_sfc <- sf::st_geometry(test_region)[0]
  expect_error(
    migration_chronology("yebsap-example", empty_sfc),
    class = "simpleError"
  )
  expect_error(
    seasonal_proportion_population("yebsap-example", empty_sfc),
    class = "simpleError"
  )
  expect_error(
    migration_chronology(
      "yebsap-example",
      test_region,
      reference_region = empty_sfc
    ),
    class = "simpleError"
  )
})

test_that("migration_chronology() and seasonal_proportion_population() accept sfc inputs", {
  rast <- suppressMessages(
    load_raster(
      "yebsap-example",
      product = "proportion-population",
      period = "seasonal",
      resolution = "27km"
    )
  )
  trimmed_sf <- sf::st_as_sf(
    terra::as.polygons(terra::ext(terra::trim(rast)), crs = terra::crs(rast))
  )
  trimmed_sfc <- sf::st_geometry(trimmed_sf)
  two_sfc <- sf::st_make_grid(trimmed_sf, n = c(2, 1))
  two_sf <- sf::st_as_sf(two_sfc)

  # single-feature sfc region with bind_regions = TRUE (no geometry column attached)
  chron_sf <- suppressMessages(
    migration_chronology("yebsap-example", trimmed_sf, resolution = "27km")
  )
  chron_sfc <- suppressMessages(
    migration_chronology(
      "yebsap-example",
      trimmed_sfc,
      resolution = "27km",
      bind_regions = TRUE
    )
  )
  expect_equal(chron_sfc, chron_sf)

  # multi-feature sfc region and sfc reference_region
  seas_sf <- suppressMessages(
    seasonal_proportion_population(
      "yebsap-example",
      two_sf,
      reference_region = trimmed_sf,
      resolution = "27km"
    )
  )
  seas_sfc <- suppressMessages(
    seasonal_proportion_population(
      "yebsap-example",
      two_sfc,
      reference_region = trimmed_sfc,
      resolution = "27km"
    )
  )
  expect_equal(seas_sfc, seas_sf)
})
