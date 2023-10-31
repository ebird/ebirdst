context("Loading trends data")

skip_on_cran()

test_that("load_trends()", {
  trends <- load_trends("yebsap-example")
  expect_is(trends, "data.frame")
  cols <- c("species_code", "season", "start_year", "end_year", "start_date",
            "end_date", "srd_id", "longitude", "latitude", "abd", "abd_ppy",
            "abd_ppy_lower", "abd_ppy_upper", "abd_ppy_nonzero", "abd_trend",
            "abd_trend_lower", "abd_trend_upper")
  expect_equal(names(trends), cols)
  expect_gt(nrow(trends), 0)

  trends_folds <- load_trends("yebsap-example", fold_estimates = TRUE)
  cols <- c("species_code", "season", "fold", "srd_id", "latitude", "longitude",
            "abd", "abd_ppy")
  expect_equal(names(trends_folds), cols)
  expect_equal(nrow(trends_folds), 100 * nrow(trends))
})
