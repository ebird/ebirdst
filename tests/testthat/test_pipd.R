context("PIs and PDs")

skip_on_cran()

test_that("load_pis", {
  pis <- load_pis(path)
  expect_is(pis, "data.frame")
  expect_gt(nrow(pis), 0)

  # columns
  expect_named(pis, c("stixel_id", "latitude", "longitude", "day_of_year",
                      "predictor", "importance"))

  # predictors
  expect_true(all(sort(unique(pis$predictor)) ==
                    sort(ebirdst_predictors$predictor)))

  # return sf object
  pis <- load_pis(path, return_sf = TRUE)
  expect_is(pis, "sf")
  expect_gt(nrow(pis), 0)
  expect_named(pis, c("stixel_id", "day_of_year",
                      "predictor", "importance", "geometry"))
  expect_true(all(sort(unique(pis$predictor)) ==
                    sort(ebirdst_predictors$predictor)))

  # invalid path
  expect_error(load_pis("/invalid/path"))
})

test_that("load_pds", {
  pds <- load_pds(path)
  expect_is(pds, "data.frame")
  expect_gt(nrow(pds), 0)

  # columns
  expect_named(pds, c("stixel_id", "latitude", "longitude", "day_of_year",
                      "predictor", "predictor_value", "response"))

  # return sf object
  pds <- load_pds(path, return_sf = TRUE)
  expect_is(pds, "sf")
  expect_gt(nrow(pds), 0)
  expect_named(pds, c("stixel_id", "day_of_year",
                      "predictor", "predictor_value", "response",
                      "geometry"))

  # invalid path
  expect_error(load_pis("/invalid/path"))
})

test_that("subset pis/pds", {
  # data frame subset
  pis <- load_pis(path)
  pds <- load_pds(path)
  e <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45),
                      t = c(0.5, 0.6))
  pi_ss <- ebirdst_subset(pis, e)
  pd_ss <- ebirdst_subset(pds, e)

  # expectations
  expect_is(pi_ss, "data.frame")
  expect_is(pd_ss, "data.frame")
  expect_lt(nrow(pi_ss), nrow(pis))
  expect_lt(nrow(pd_ss), nrow(pds))
  expect_equal(ncol(pis), ncol(pi_ss))
  expect_equal(ncol(pds), ncol(pd_ss))
  expect_named(pi_ss, names(pis))
  expect_named(pd_ss, names(pds))

  # data frame subset
  pis <- load_pis(path, return_sf = TRUE)
  pds <- load_pds(path, return_sf = TRUE)
  pi_ss <- ebirdst_subset(pis, e)
  pd_ss <- ebirdst_subset(pds, e)

  # expectations
  expect_is(pi_ss, "sf")
  expect_is(pd_ss, "sf")
  expect_lt(nrow(pi_ss), nrow(pis))
  expect_lt(nrow(pd_ss), nrow(pds))
  expect_equal(ncol(pis), ncol(pi_ss))
  expect_equal(ncol(pds), ncol(pd_ss))
  expect_named(pi_ss, names(pis))
  expect_named(pd_ss, names(pds))
})

test_that("plot_pis", {
  pis <- load_pis(path)
  e <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45),
                      t = c(0.5, 0.6))

  # expectations
  top_pi <- plot_pis(pis = pis, ext = e, plot = FALSE)
  expect_is(top_pi, "data.frame")
  expect_true(all(top_pi$predictor %in% ebirdst_predictors$lc_class_label))

  # check names and grouping work correctly
  top_pi <- plot_pis(pis = pis, ext = e,
                     pretty_names = FALSE, by_cover_class = TRUE,
                     plot = FALSE)
  expect_true("mcd12q1_lccs1_fs_c14" %in% top_pi$predictor)
  top_pi <- plot_pis(pis = pis, ext = e,
                     pretty_names = TRUE, by_cover_class = FALSE,
                     plot = FALSE)
  expect_true("Deciduous Broadleaf Forests PLAND" %in% top_pi$predictor)
  top_pi <- plot_pis(pis = pis, ext = e,
                     pretty_names = FALSE, by_cover_class = FALSE,
                     plot = FALSE)
  expect_true("mcd12q1_lccs1_fs_c14_1500_pland" %in% top_pi$predictor)

  # checking length
  top_pi <- plot_pis(pis = pis, ext = e, n_top_pred = 10, plot = FALSE)
  expect_length(unique(top_pi$predictor), 10)
  top_pi <- plot_pis(pis = pis, ext = e, n_top_pred = 10,
                     by_cover_class = FALSE, plot = FALSE)
  expect_length(unique(top_pi$predictor), 10)

  # not enough num_top_preds
  expect_error(plot_pis(pis = pis, ext = e, n_top_pred = 1,
                        plot = FALSE))

  # broken input
  expect_error(plot_pis(pis = pis))
  expect_error(plot_pis(pis = 1:10, ext = e))
  expect_error(plot_pis(pis = data.frame(), ext = e))
})

test_that("plot_pds", {
  pds <- load_pds(path)
  e <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45),
                      t = c(0.5, 0.6))

  # expectations
  pd_smooth <- plot_pds(pds, "solar_noon_diff_mid", ext = e, n_bs = 5,
                        plot = FALSE)
  expect_is(pd_smooth, "data.frame")
  expect_equal(nrow(pd_smooth), 25)
  expect_named(pd_smooth, c("x", "pd_median", "pd_lower", "pd_upper"))

  # broken input
  expect_error(plot_pds(pds, "invalid_pred", ext = e))
  expect_error(plot_pds(pds = 1:10, ext = e))
  expect_error(plot_pds(pds = data.frame(), ext = e))

  # missing temporal info
  nt <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45))
  expect_error(plot_pds(pds, "solar_noon_diff", ext = nt, n_bs = 5))
})
