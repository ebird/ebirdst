context("Utility functions")

skip_on_cran()

test_that("get_species()", {
  expect_equal(get_species("Wood Thrush"), "woothr")
  expect_equal(get_species("Yellow-bellied Sapsucker"), "yebsap")
  expect_equal(get_species("yebsap-example"), "yebsap-example")
  expect_equal(get_species(character(0)), character(0))
  expect_equal(get_species(NA_character_), NA_character_)
  expect_equal(get_species("aakspa1"), NA_character_)
})
