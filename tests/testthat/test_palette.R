context("Color palettes")

test_that("ebirdst_palettes", {
  expect_is(ebirdst_palettes(n = 10), "character")
  expect_length(ebirdst_palettes(n = 10), 10)
  expect_length(ebirdst_palettes(n = 10, type = "breeding"), 10)
  expect_length(ebirdst_palettes(n = 10, type = "trends"), 10)
  expect_match(ebirdst_palettes(n = 10, type = "nonbreeding"),
               "#[0-9A-F]{6}")
})

test_that("abundance_pallete throws warning and matches ebirdst_palletes", {
  expect_warning(p <- abundance_palette(10, "weekly"), regexp = "is deprecated")
  expect_equal(p, ebirdst_palettes(10, "weekly"))
})
