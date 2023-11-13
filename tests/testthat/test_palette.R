context("Color palettes")

test_that("ebirdst_palettes", {
  expect_is(ebirdst_palettes(n = 10), "character")
  expect_length(ebirdst_palettes(n = 10), 10)
  expect_length(ebirdst_palettes(n = 10, type = "breeding"), 10)
  expect_length(ebirdst_palettes(n = 10, type = "trends"), 10)
  expect_match(ebirdst_palettes(n = 10, type = "nonbreeding"),
               "#[0-9A-F]{6}")
})
