context("Grid sampling")

# synthetic observation data used across tests
set.seed(1)
n_obs <- 5000L
checklists <- data.frame(
  longitude   = runif(n_obs, -90, -70),
  latitude    = runif(n_obs, 35, 50),
  day_of_year = sample.int(365L, n_obs, replace = TRUE),
  year        = sample(2016L:2020L, n_obs, replace = TRUE),
  obs         = sample(c(0L, 1L), n_obs, replace = TRUE, prob = c(0.9, 0.1))
)

test_that("assign_to_grid() - projected XYT", {
  pts <- data.frame(x = runif(100), y = runif(100),
                    t = sample.int(52L, 100L, replace = TRUE))
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
  pts <- data.frame(longitude = runif(50, -90, -70),
                    latitude  = runif(50, 35, 50))
  cells <- assign_to_grid(pts, res = c(10000, 10000),
                          is_lonlat = TRUE, jitter_grid = FALSE)
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
  sampled_id <- grid_sample(checklists, keep_cell_id = TRUE, jitter_grid = FALSE)
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
  sampled_max <- grid_sample_stratified(checklists, maximum_ss = 500L,
                                        jitter_grid = FALSE)
  expect_lte(nrow(sampled_max), 500L)
})
