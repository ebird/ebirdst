skip_on_cran()

test_that("assign_weeks_to_seasons()", {
  seasons <- assign_weeks_to_seasons("yebsap-example", min_quality = 3)
  expect_length(seasons, 52L)
  expect_type(seasons, "character")
  # all four seasons of the example species score a 3
  expect_setequal(
    unique(seasons),
    c(
      "breeding",
      "nonbreeding",
      "prebreeding_migration",
      "postbreeding_migration"
    )
  )
  # a lower quality threshold can only add weeks, never remove them
  relaxed <- assign_weeks_to_seasons("yebsap-example", min_quality = 1)
  expect_equal(relaxed, seasons)

  # invalid arguments
  expect_error(assign_weeks_to_seasons("yebsap-example", min_quality = 0))
  expect_error(assign_weeks_to_seasons("yebsap-example", min_quality = 4))
  expect_error(assign_weeks_to_seasons("yebsap-example", min_quality = 2.5))
  expect_error(assign_weeks_to_seasons("Yellow-bellied Sapsuckr"))
})

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

test_that("is_integer(), is_count(), and is_flag()", {
  expect_true(is_integer(1))
  expect_true(is_integer(c(-2, 0, 3)))
  expect_true(is_integer(.Machine$integer.max))
  expect_false(is_integer(1.5))
  expect_false(is_integer(NA_integer_))
  expect_false(is_integer(Inf))
  expect_false(is_integer("1"))
  # values beyond integer range are rejected without warning about coercion
  expect_silent(expect_false(is_integer(1e10)))
  expect_silent(expect_false(is_integer(-1e10)))

  expect_true(is_count(0))
  expect_false(is_count(-1))
  expect_false(is_count(c(1, 2)))
  expect_silent(expect_false(is_count(1e10)))

  expect_true(is_flag(TRUE))
  expect_false(is_flag(NA))
  expect_false(is_flag(c(TRUE, FALSE)))
  expect_false(is_flag(1))
})

test_that("calculate_mcc_f1()", {
  skip_if_not_installed("PresenceAbsence")

  obs <- c(1, 1, 0, 0)
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
