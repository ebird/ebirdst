context("Grid sampling")

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
  pts1 <- data.frame(x = runif(50), y = runif(50))
  pts2 <- data.frame(x = runif(50), y = runif(50))
  cells1 <- assign_to_grid(pts1, res = c(0.2, 0.2), jitter_grid = FALSE)
  gd <- attr(cells1, "grid_definition")
  cells2 <- assign_to_grid(pts2, grid_definition = gd)
  expect_equal(attr(cells2, "grid_definition"), gd)
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
  expect_true(all(unique(x$year) %in% capped$year))
})

test_that("cap_cells_by_quantile() caps over-sampled cells without case control", {
  # one busy cell (20 obs) and one quiet cell (3 obs)
  sampled <- data.frame(
    x = c(rep(1, 20), rep(100, 3)),
    y = c(rep(1, 20), rep(100, 3)),
    t = 1
  )

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
  expect_equal(sum(capped$x == 1), 12L)
  expect_equal(sum(capped$x == 100), 3L)
})

test_that("cap_cells_by_quantile() caps detections and non-detections independently", {
  # busy cell (x = 1): 10 detections, 15 non-detections
  # quiet cell (x = 100): 1 detection, 1 non-detection
  sampled <- data.frame(
    x = c(rep(1, 10), rep(1, 15), rep(100, 1), rep(100, 1)),
    y = c(rep(1, 10), rep(1, 15), rep(100, 1), rep(100, 1)),
    t = 1,
    .detected = c(
      rep(TRUE, 10),
      rep(FALSE, 15),
      rep(TRUE, 1),
      rep(FALSE, 1)
    )
  )

  set.seed(1)
  capped <- cap_cells_by_quantile(
    sampled,
    prob = 0.5,
    coords = c("x", "y", "t"),
    res_xy = c(10, 10),
    case_control = TRUE
  )

  # busy cell trimmed independently for each class
  expect_equal(sum(capped$x == 1 & capped$.detected), 6L)
  expect_equal(sum(capped$x == 1 & !capped$.detected), 8L)
  # quiet cell is untouched, a single detection is never dropped
  expect_equal(sum(capped$x == 100 & capped$.detected), 1L)
  expect_equal(sum(capped$x == 100 & !capped$.detected), 1L)
})

test_that("cap_cells_by_quantile() retains rare sample_by levels", {
  # busy cell (x = 1): 20 obs, one carrying a rare island and one carrying a
  # rare year; quiet cell (x = 100): 3 common obs, left untouched
  sampled <- data.frame(
    x = c(rep(1, 20), rep(100, 3)),
    y = c(rep(1, 20), rep(100, 3)),
    t = 1,
    island = c(
      "rare_island",
      "common",
      rep("common", 18),
      rep("common", 3)
    ),
    year = c(2000L, 2050L, rep(2000L, 18), rep(2000L, 3))
  )

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
  expect_equal(sum(capped$x == 1), 12L)
  expect_true("rare_island" %in% capped$island)
  expect_true(2050L %in% capped$year)
})

test_that("cap_cells_by_quantile() keeps all protected rows even if that exceeds the cap", {
  # busy cell (x = 1): 20 obs, 15 of which have a distinct rare island level
  sampled <- data.frame(
    x = c(rep(1, 20), rep(100, 3)),
    y = c(rep(1, 20), rep(100, 3)),
    t = 1,
    island = c(paste0("island", seq_len(15)), rep("common", 5), rep("common", 3))
  )

  set.seed(1)
  capped <- cap_cells_by_quantile(
    sampled,
    prob = 0.5,
    coords = c("x", "y", "t"),
    res_xy = c(10, 10),
    case_control = FALSE,
    sample_by = "island"
  )

  # cap_n for the busy cell is 12, but 15 rare islands must all be retained
  expect_equal(sum(capped$x == 1), 15L)
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
