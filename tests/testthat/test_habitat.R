context("Habitat associations")

skip_on_cran()

e <- ebirdst_extent(c(xmin = -86, xmax = -85, ymin = 42, ymax = 43))
habitat <- ebirdst_habitat(path, ext = e)

test_that("ebirdst_habitat", {
  # expected
  expect_is(habitat, "ebirdst_habitat")
  expect_is(habitat, "data.frame")
  expect_named(habitat, c("predictor", "week", "importance",
                          "prob_pos_slope", "direction"))
  expect_is(habitat$predictor, "character")
  expect_is(habitat$week, "character")
  expect_is(habitat$importance, "numeric")
  expect_is(habitat$direction, "numeric")
  expect_length(unique(habitat$week), 52)

  # invalid inputs
  expect_error(ebirdst_habitat("/invalid/path", ext = e))
})

test_that("ebirdst_habitat extent", {
  # temporal extent ignored
  e <- ebirdst_extent(c(xmin = -86, xmax = -85, ymin = 42, ymax = 43),
                      t = c(0.5, 0.75))
  h2 <- ebirdst_habitat(path, ext = e)
  expect_identical(habitat, h2)

  # spatial extent required
  expect_error(ebirdst_habitat(path))
})

test_that("ebirdst_habitat pland_only", {
  predictors <- unique(habitat$predictor)
  expect_true(all(!stringr::str_detect(predictors, "(ed|sd)$")))
})

test_that("plot ebirdst_habitat", {
  suppressWarnings({g <- plot(habitat)})
  expect_is(g, "gg")
})
