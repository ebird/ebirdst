context("Utility functions")

skip_on_cran()

test_that("get_species empty vector", {
  expect_equal(get_species(character(0)), character(0))
})

test_that("get_species NA_character_", {
  expect_equal(get_species(NA_character_), NA_character_)
})
