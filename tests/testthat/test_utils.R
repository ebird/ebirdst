context("Utility functions")

skip_on_cran()

test_that("get_species()", {
  expect_equal(get_species("Wood Thrush"), "woothr")
  expect_equal(get_species("Yellow-bellied Sapsucker"), "yebsap")
  expect_equal(get_species("yebsap-example"), "yebsap-example")
  expect_equal(get_species(character(0)), character(0))
  expect_equal(get_species(NA_character_), NA_character_)
  expect_equal(get_species("aakspa1"), NA_character_)
  # case-insensitive lookup
  expect_equal(get_species("wood thrush"), "woothr")
})

test_that("calculate_mcc_f1()", {
  skip_if_not_installed("PresenceAbsence")

  obs  <- c(1, 1, 0, 0)
  pred <- c(1, 0, 0, 1)
  r <- calculate_mcc_f1(obs > 0, pred > 0)
  expect_equal(r$f1, 0.5)
  expect_equal(r$mcc, 0)

  # perfect predictions
  r_perfect <- calculate_mcc_f1(obs > 0, obs > 0)
  expect_equal(r_perfect$f1, 1)
  expect_equal(r_perfect$mcc, 1)

  # degenerate: no positive predictions -> 0, not NaN
  r_nopred <- calculate_mcc_f1(c(TRUE, FALSE), c(FALSE, FALSE))
  expect_equal(r_nopred$f1, 0)
  expect_equal(r_nopred$mcc, 0)

  # degenerate: no positive observations -> 0, not NaN
  r_noobs <- calculate_mcc_f1(c(FALSE, FALSE), c(TRUE, FALSE))
  expect_equal(r_noobs$f1, 0)
  expect_equal(r_noobs$mcc, 0)
})

test_that("date_to_st_week()", {
  # Jan 4 is the center of week 1; Jan 1 also falls in week 1
  expect_equal(date_to_st_week(as.Date("2022-01-04")), 1L)
  expect_equal(date_to_st_week(as.Date("2022-01-01")), 1L)
  # Dec 31 should be week 52
  expect_equal(date_to_st_week(as.Date("2022-12-31")), 52L)
  # all dates in a year map to 1:52
  dates <- seq(as.Date("2022-01-01"), as.Date("2022-12-31"), by = 1)
  weeks <- date_to_st_week(dates)
  expect_true(all(weeks >= 1L & weeks <= 52L))
  # vectorized: returns same length as input
  d <- as.Date(c("2022-01-04", "2022-06-15", "2022-12-31"))
  expect_length(date_to_st_week(d), 3L)
  # version 2021 returns integers in 1:52
  expect_type(date_to_st_week(as.Date("2020-06-15"), version = 2021), "integer")
  # invalid version errors
  expect_error(date_to_st_week(as.Date("2022-01-01"), version = 2020))
})
