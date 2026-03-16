context("Color palettes")

test_that("ebirdst_palettes", {
  expect_type(ebirdst_palettes(n = 10), "character")
  expect_length(ebirdst_palettes(n = 10), 10)
  expect_length(ebirdst_palettes(n = 10, type = "breeding"), 10)
  expect_length(ebirdst_palettes(n = 10, type = "trends"), 10)
  expect_match(ebirdst_palettes(n = 10, type = "nonbreeding"),
               "#[0-9A-F]{6}")
  # all remaining types return n colors
  for (type in c("migration", "postbreeding_migration",
                 "prebreeding_migration", "year_round")) {
    expect_length(ebirdst_palettes(n = 10, type = type), 10)
  }
  # invalid type errors
  expect_error(ebirdst_palettes(n = 10, type = "invalid"))
  # n must be >= 1
  expect_error(ebirdst_palettes(n = 0))
})

test_that("abundance_palette throws warning and matches ebirdst_palettes", {
  expect_warning(p <- abundance_palette(10, "weekly"), regexp = "is deprecated")
  expect_equal(p, ebirdst_palettes(10, "weekly"))
})
