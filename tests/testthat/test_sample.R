# synthetic observation data used across tests
set.seed(1)
n_obs <- 5000L
checklists <- data.frame(
  longitude = runif(n_obs, -90, -70),
  latitude = runif(n_obs, 35, 50),
  day_of_year = sample.int(365L, n_obs, replace = TRUE),
  year = sample(2016L:2020L, n_obs, replace = TRUE),
  obs = sample(c(0L, 1L), n_obs, replace = TRUE, prob = c(0.9, 0.1))
)

test_that("assign_to_grid() - projected XYT", {
  pts <- data.frame(
    x = runif(100),
    y = runif(100),
    t = sample.int(52L, 100L, replace = TRUE)
  )
  cells <- assign_to_grid(pts, res = c(0.1, 0.1, 7), jitter_grid = FALSE)
  expect_s3_class(cells, "data.frame")
  expect_equal(nrow(cells), 100L)
  expect_true("cell_xy" %in% names(cells))
  expect_true("cell_xyt" %in% names(cells))
  expect_false(anyNA(cells$cell_xy))
  expect_false(anyNA(cells$cell_xyt))
  expect_false(is.null(attr(cells, "grid_definition")))
})

test_that("assign_to_grid() - lonlat spatial-only", {
  pts <- data.frame(
    longitude = runif(50, -90, -70),
    latitude = runif(50, 35, 50)
  )
  cells <- assign_to_grid(
    pts,
    res = c(10000, 10000),
    is_lonlat = TRUE,
    jitter_grid = FALSE
  )
  expect_equal(nrow(cells), 50L)
  expect_true("cell_xy" %in% names(cells))
  expect_false("cell_xyt" %in% names(cells))
})

test_that("assign_to_grid() - grid_definition reuse", {
  # pts1 spans the full [0, 1] range so pts2, also drawn from [0, 1], is
  # guaranteed to fall inside pts1's bounding box
  pts1 <- data.frame(x = c(0, 1, runif(48)), y = c(0, 1, runif(48)))
  pts2 <- data.frame(x = runif(50), y = runif(50))
  cells1 <- assign_to_grid(pts1, res = c(0.2, 0.2), jitter_grid = FALSE)
  gd <- attr(cells1, "grid_definition")
  cells2 <- assign_to_grid(pts2, grid_definition = gd)
  expect_equal(attr(cells2, "grid_definition"), gd)
})

test_that("assign_to_grid() errors when points fall below the grid origin", {
  pts1 <- data.frame(x = runif(50, min = 1), y = runif(50, min = 1))
  cells1 <- assign_to_grid(pts1, res = c(0.2, 0.2), jitter_grid = FALSE)
  gd <- attr(cells1, "grid_definition")
  pts2 <- data.frame(x = c(0.5, -1), y = c(0.5, 0.5))
  expect_error(
    assign_to_grid(pts2, grid_definition = gd),
    "outside the provided spatial grid"
  )
})

test_that("grid_sample()", {
  set.seed(1)
  sampled <- grid_sample(checklists, jitter_grid = FALSE)
  expect_s3_class(sampled, "data.frame")
  expect_lt(nrow(sampled), nrow(checklists))
  # all input columns are preserved
  expect_true(all(names(checklists) %in% names(sampled)))

  # keep_cell_id adds .cell_id column
  set.seed(1)
  sampled_id <- grid_sample(
    checklists,
    keep_cell_id = TRUE,
    jitter_grid = FALSE
  )
  expect_true(".cell_id" %in% names(sampled_id))

  # empty input returns empty output with the same column structure
  empty <- checklists[0, ]
  out <- grid_sample(empty)
  expect_equal(nrow(out), 0L)
  expect_equal(names(out), names(checklists))
})

test_that("grid_sample() spatial-only grid", {
  # a 2-element res drops the time dimension, so it thins far more than a
  # spacetime grid at the same spatial resolution
  set.seed(1)
  sampled <- grid_sample(
    checklists,
    res = c(100000, 100000),
    jitter_grid = FALSE
  )
  set.seed(1)
  sampled_xyt <- grid_sample(
    checklists,
    res = c(100000, 100000, 30),
    jitter_grid = FALSE
  )
  expect_s3_class(sampled, "data.frame")
  expect_true(all(names(checklists) %in% names(sampled)))
  expect_lt(nrow(sampled), nrow(sampled_xyt))

  # naming only the spatial coordinates gives the same result
  set.seed(1)
  sampled_coords2 <- grid_sample(
    checklists,
    coords = c("longitude", "latitude"),
    res = c(100000, 100000),
    jitter_grid = FALSE
  )
  expect_equal(sampled_coords2, sampled)

  # cell ids have no time component
  sampled_id <- grid_sample(
    checklists,
    res = c(100000, 100000),
    jitter_grid = FALSE,
    keep_cell_id = TRUE
  )
  expect_true(all(lengths(strsplit(sampled_id$.cell_id, "-")) == 2L))

  # a spacetime res needs a temporal coordinate
  expect_error(
    grid_sample(checklists, coords = c("longitude", "latitude")),
    "spatial-only sampling"
  )
})

test_that("grid_sample_stratified()", {
  set.seed(1)
  sampled <- grid_sample_stratified(checklists, jitter_grid = FALSE)
  expect_s3_class(sampled, "data.frame")
  expect_lt(nrow(sampled), nrow(checklists))
  expect_true(all(names(checklists) %in% names(sampled)))

  # maximum_ss is respected
  set.seed(1)
  sampled_max <- grid_sample_stratified(
    checklists,
    maximum_ss = 500L,
    jitter_grid = FALSE
  )
  expect_lte(nrow(sampled_max), 500L)

  # a maximum larger than the sample is a no-op rather than being subsampled
  set.seed(1)
  sampled_high_max <- grid_sample_stratified(
    checklists,
    maximum_ss = nrow(checklists) * 2L,
    jitter_grid = FALSE
  )
  expect_equal(nrow(sampled_high_max), nrow(sampled))
})

test_that("grid_sample_stratified() spatial-only grid", {
  res_xy <- c(100000, 100000)

  set.seed(1)
  sampled <- grid_sample_stratified(
    checklists,
    res = res_xy,
    jitter_grid = FALSE
  )
  set.seed(1)
  sampled_xyt <- grid_sample_stratified(
    checklists,
    res = c(res_xy, 30),
    jitter_grid = FALSE
  )
  expect_true(all(names(checklists) %in% names(sampled)))
  expect_lt(nrow(sampled), nrow(sampled_xyt))

  # a unified spatial-only grid yields cell ids with no time component
  sampled_id <- grid_sample_stratified(
    checklists,
    res = res_xy,
    unified_grid = TRUE,
    keep_cell_id = TRUE,
    jitter_grid = FALSE
  )
  expect_true(all(lengths(strsplit(sampled_id$.cell_id, "-")) == 2L))

  # detection oversampling also works on a spatial-only grid
  set.seed(1)
  sampled_cc <- grid_sample_stratified(
    checklists,
    res = res_xy,
    min_detection_probability = 0.3,
    jitter_grid = FALSE
  )
  expect_gte(mean(sampled_cc$obs > 0), 0.29)

  # with no strata, coords are passed through to grid_sample()
  set.seed(1)
  no_strata <- grid_sample_stratified(
    checklists,
    coords = c("longitude", "latitude"),
    res = res_xy,
    by_year = FALSE,
    case_control = FALSE,
    jitter_grid = FALSE
  )
  set.seed(1)
  expect_equal(
    no_strata,
    grid_sample(
      checklists,
      coords = c("longitude", "latitude"),
      res = res_xy,
      jitter_grid = FALSE
    )
  )
})

test_that("grid_sample_stratified() unified grid uses res passed via ...", {
  set.seed(1)
  coarse <- grid_sample_stratified(
    checklists,
    unified_grid = TRUE,
    keep_cell_id = TRUE,
    res = c(100000, 100000, 30)
  )
  set.seed(1)
  fine <- grid_sample_stratified(
    checklists,
    unified_grid = TRUE,
    keep_cell_id = TRUE,
    res = c(3000, 3000, 7)
  )
  expect_lt(length(unique(coarse$.cell_id)), length(unique(fine$.cell_id)))
})

test_that("grid_sample_stratified() oversamples detections to reach min_detection_probability", {
  set.seed(1)
  x <- checklists
  x$obs <- sample(c(0L, 1L), nrow(x), replace = TRUE, prob = c(0.98, 0.02))

  # the grid assignment for the oversampling pool doesn't change across
  # resampling iterations, only the random draw does, so it's computed once
  # rather than on every one of up to 25 iterations; with by_year = FALSE the
  # initial stratified sample splits by case control class only (2 calls),
  # and the oversampling setup adds exactly 1 more call, regardless of how
  # many iterations it takes to reach the target
  original_grid_cell_id <- grid_cell_id
  call_count <- 0L
  local_mocked_bindings(
    grid_cell_id = function(...) {
      call_count <<- call_count + 1L
      return(original_grid_cell_id(...))
    }
  )

  sampled <- grid_sample_stratified(
    x,
    by_year = FALSE,
    min_detection_probability = 0.3,
    jitter_grid = FALSE
  )

  expect_equal(call_count, 3L)
  expect_gte(mean(sampled$obs > 0), 0.29)
})

test_that("grid_sample_stratified() errors on missing values in stratifying columns", {
  x <- checklists
  x$year[1] <- NA
  expect_error(grid_sample_stratified(x, jitter_grid = FALSE), "missing values")

  x2 <- checklists
  x2$island <- "a"
  x2$island[1] <- NA
  expect_error(
    grid_sample_stratified(x2, sample_by = "island", jitter_grid = FALSE),
    "missing values"
  )
})

test_that("grid_sample_stratified() errors on missing values in obs_column", {
  x <- checklists
  x$obs[1] <- NA
  expect_error(
    grid_sample_stratified(x, case_control = TRUE, jitter_grid = FALSE),
    "missing values"
  )
})

test_that("grid_sample_stratified() validates cell_quantile_cap", {
  expect_error(grid_sample_stratified(checklists, cell_quantile_cap = 0))
  expect_error(grid_sample_stratified(checklists, cell_quantile_cap = -0.1))
  expect_error(grid_sample_stratified(checklists, cell_quantile_cap = 1.5))
  expect_error(grid_sample_stratified(checklists, cell_quantile_cap = c(0.5, 0.6)))
  # a value of 1 is a valid no-op, not an error
  expect_no_error(grid_sample_stratified(checklists, cell_quantile_cap = 1))
})

test_that("grid_sample_stratified() cell_quantile_cap is a no-op backwards compatible default", {
  set.seed(1)
  baseline <- grid_sample_stratified(checklists, jitter_grid = FALSE)
  set.seed(1)
  with_null <- grid_sample_stratified(
    checklists,
    cell_quantile_cap = NULL,
    jitter_grid = FALSE
  )
  set.seed(1)
  with_one <- grid_sample_stratified(
    checklists,
    cell_quantile_cap = 1,
    jitter_grid = FALSE
  )
  expect_equal(with_null, baseline)
  expect_equal(with_one, baseline)
})

test_that("grid_sample_stratified() cell_quantile_cap reduces over-sampled sites", {
  # inject a chronically over-sampled site (e.g. a feeder) into the data
  extra <- data.frame(
    longitude = -80,
    latitude = 40,
    day_of_year = sample.int(365L, 60L, replace = TRUE),
    year = sample(2016L:2020L, 60L, replace = TRUE),
    obs = sample(c(0L, 1L), 60L, replace = TRUE, prob = c(0.5, 0.5))
  )
  x <- rbind(checklists, extra)

  set.seed(1)
  uncapped <- grid_sample_stratified(x, jitter_grid = FALSE)
  set.seed(1)
  capped <- grid_sample_stratified(
    x,
    cell_quantile_cap = 0.5,
    jitter_grid = FALSE
  )

  expect_lte(nrow(capped), nrow(uncapped))
  n_site_uncapped <- sum(uncapped$longitude == -80 & uncapped$latitude == 40)
  n_site_capped <- sum(capped$longitude == -80 & capped$latitude == 40)
  expect_lt(n_site_capped, n_site_uncapped)
})

test_that("grid_sample_stratified() cell_quantile_cap retains rare sample_by levels", {
  # inject a chronically over-sampled site, with a single observation of a
  # rare island level hiding among it, to confirm the cap doesn't drop it
  x <- checklists
  x$island <- "common"
  extra <- data.frame(
    longitude = -80,
    latitude = 40,
    day_of_year = sample.int(365L, 60L, replace = TRUE),
    year = sample(2016L:2020L, 60L, replace = TRUE),
    obs = sample(c(0L, 1L), 60L, replace = TRUE, prob = c(0.5, 0.5)),
    island = c("rare_island", rep("common", 59))
  )
  x <- rbind(x, extra)

  set.seed(1)
  capped <- grid_sample_stratified(
    x,
    sample_by = "island",
    cell_quantile_cap = 0.5,
    jitter_grid = FALSE
  )

  expect_true("rare_island" %in% capped$island)
})

test_that("grid_sample_stratified() cell_quantile_cap can drop rare year levels", {
  # a chronically over-sampled site with one observation per year, alongside
  # a rare island level; background checklist data supplies enough occupied
  # cells for cap_cells_by_quantile()'s np >= 10 guard to trust the quantile,
  # so the site is thinned, dropping most years in the process, while the
  # rare island level is still always retained
  busy <- data.frame(
    longitude = -80,
    latitude = 40,
    day_of_year = 100,
    year = 2000:2019,
    island = c("rare_island", rep("common", 19))
  )
  quiet <- data.frame(
    longitude = seq(-79, -70, length.out = 10),
    latitude = seq(41, 50, length.out = 10),
    day_of_year = 100,
    year = 2000:2009,
    island = "common"
  )
  bg <- checklists
  bg$island <- "common"
  x <- rbind(busy, quiet, bg[, names(busy)])

  set.seed(1)
  capped <- grid_sample_stratified(
    x,
    sample_by = "island",
    cell_quantile_cap = 0.5,
    case_control = FALSE,
    jitter_grid = FALSE
  )
  busy_capped <- capped[capped$longitude == -80 & capped$latitude == 40, ]

  expect_true("rare_island" %in% busy_capped$island)
  expect_lt(length(unique(busy_capped$year)), length(unique(busy$year)))
})

# 25 filler cells with per-cell counts 1:25, plus a quiet and a busy cell;
# the filler spread guarantees >= 10 occupied cells exceed the resulting
# cap, satisfying cap_cells_by_quantile()'s np >= 10 guard on trusting the
# quantile estimate
make_ramp_cells <- function(quiet_n, busy_n, filler = seq_len(25)) {
  counts <- c(filler, quiet_n, busy_n)
  return(data.frame(
    x = rep(seq_along(counts) * 100, times = counts),
    y = 1,
    t = 1
  ))
}
quiet_x <- 2600
busy_x <- 2700

test_that("cap_cells_by_quantile() caps over-sampled cells without case control", {
  # quiet cell (2 obs) stays below the cap, busy cell (30 obs) is trimmed
  sampled <- make_ramp_cells(quiet_n = 2, busy_n = 30)

  set.seed(1)
  capped <- cap_cells_by_quantile(
    sampled,
    prob = 0.5,
    coords = c("x", "y", "t"),
    res_xy = c(10, 10),
    case_control = FALSE
  )

  expect_true(all(names(sampled) %in% names(capped)))
  # busy cell trimmed down to the cap, quiet cell left untouched
  expect_equal(sum(capped$x == busy_x), 13L)
  expect_equal(sum(capped$x == quiet_x), 2L)
})

test_that("cap_cells_by_quantile() caps detections and non-detections independently", {
  # detections ramp to a busy cell of 30 obs, non-detections ramp (filler
  # counts doubled) to a busy cell of 60 obs, so each class gets its own cap
  true_cells <- make_ramp_cells(quiet_n = 2, busy_n = 30)
  true_cells$.detected <- TRUE
  false_cells <- make_ramp_cells(
    quiet_n = 3,
    busy_n = 60,
    filler = seq_len(25) * 2
  )
  false_cells$.detected <- FALSE
  sampled <- rbind(true_cells, false_cells)

  set.seed(1)
  capped <- cap_cells_by_quantile(
    sampled,
    prob = 0.5,
    coords = c("x", "y", "t"),
    res_xy = c(10, 10),
    case_control = TRUE
  )

  # busy cell trimmed independently, and to a different cap, for each class
  expect_equal(sum(capped$x == busy_x & capped$.detected), 13L)
  expect_equal(sum(capped$x == busy_x & !capped$.detected), 26L)
  # quiet cell is below the cap in both classes and left untouched
  expect_equal(sum(capped$x == quiet_x & capped$.detected), 2L)
  expect_equal(sum(capped$x == quiet_x & !capped$.detected), 3L)
})

test_that("cap_cells_by_quantile() retains rare sample_by levels", {
  # busy cell (30 obs): one obs carries a rare island, another a rare year;
  # quiet cell (2 obs) stays below the cap and is left untouched
  sampled <- make_ramp_cells(quiet_n = 2, busy_n = 30)
  sampled$island <- "common"
  sampled$year <- 2000L
  busy_rows <- which(sampled$x == busy_x)
  sampled$island[busy_rows[1]] <- "rare_island"
  sampled$year[busy_rows[2]] <- 2050L

  set.seed(1)
  capped <- cap_cells_by_quantile(
    sampled,
    prob = 0.5,
    coords = c("x", "y", "t"),
    res_xy = c(10, 10),
    case_control = FALSE,
    sample_by = c("island", "year")
  )

  # the cap is still enforced since the protected rows fit within it
  expect_equal(sum(capped$x == busy_x), 13L)
  expect_true("rare_island" %in% capped$island)
  expect_true(2050L %in% capped$year)
})

test_that("cap_cells_by_quantile() keeps all protected rows even if that exceeds the cap", {
  # busy cell (30 obs): 15 of which have a distinct rare island level
  sampled <- make_ramp_cells(quiet_n = 2, busy_n = 30)
  sampled$island <- "common"
  busy_rows <- which(sampled$x == busy_x)
  sampled$island[busy_rows[seq_len(15)]] <- paste0("island", seq_len(15))

  set.seed(1)
  capped <- cap_cells_by_quantile(
    sampled,
    prob = 0.5,
    coords = c("x", "y", "t"),
    res_xy = c(10, 10),
    case_control = FALSE,
    sample_by = "island"
  )

  # cap_n for the busy cell is 13, but 15 rare islands must all be retained
  expect_equal(sum(capped$x == busy_x), 15L)
  expect_true(all(paste0("island", seq_len(15)) %in% capped$island))
})

test_that("cap_cells_by_quantile() handles empty input", {
  sampled <- data.frame(x = numeric(0), y = numeric(0), t = numeric(0))
  capped <- cap_cells_by_quantile(
    sampled,
    prob = 0.5,
    coords = c("x", "y", "t"),
    res_xy = c(10, 10),
    case_control = FALSE
  )
  expect_equal(nrow(capped), 0L)
  expect_equal(names(capped), names(sampled))
})
