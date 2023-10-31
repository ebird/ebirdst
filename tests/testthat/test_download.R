context("Data download")

skip_on_cran()
skip_if_offline()

# only run this test if an ebirdst access key is present
key <- Sys.getenv("EBIRDST_KEY")
skip_if_not(!is.na(key) && key != "" && nchar(key) > 0,
            message = "Missing ebirdst access key")

# status ----

test_that("ebirdst_download_status()", {
  suppressMessages({
    files <- ebirdst_download_status("leafly", dry_run = TRUE)
  })

  expect_true(any(grepl("config.json$", files)))
  expect_true(any(grepl("proportion-population", files)))
  expect_true(any(grepl("abundance", files)))
  expect_false(any(grepl("web_download", files)))
  expect_false(any(grepl("leafly2", files)))
  # occurrence and count not downloaded by default
  expect_false(any(grepl("occurrence", files)))
  expect_false(any(grepl("count", files)))
  # pis and ppms not downloaded by defaul
  expect_false(any(grepl("pis", files)))
  expect_false(any(grepl("ppms", files)))
  # trends not downloaded
  expect_false(any(grepl("\\.parquet$", files)))
  expect_false(any(grepl("trends", files)))

  suppressMessages({
    files_all <- ebirdst_download_status("leafly",
                                         download_all = TRUE,
                                         dry_run = TRUE)
  })
  expect_true(any(grepl("config.json$", files_all)))
  expect_true(any(grepl("occurrence", files_all)))
  expect_true(any(grepl("count", files_all)))
  expect_true(any(grepl("ppms", files_all)))
  expect_true(any(grepl("pis", files_all)))
  expect_true(any(grepl("regional_stats.csv", files_all)))
  expect_false(any(grepl("web_download", files_all)))
  expect_false(any(grepl("leafly2", files_all)))
  # trends not downloaded
  expect_false(any(grepl("\\.parquet$", files_all)))
  expect_false(any(grepl("trends", files_all)))
})

test_that("ebirdst_download_status() pattern argument", {
  # only download low resolution
  suppressMessages({
    files <- ebirdst_download_status("leafly",
                                     pattern = "_27km_",
                                     dry_run = TRUE)
  })
  expect_true(any(grepl("config.json$", files)))
  expect_true(any(grepl("27km", files)))
  expect_false(any(grepl("9km", files)))
  expect_false(any(grepl("3km", files)))
})

test_that("ebirdst_download_status() missing species", {
  # only download low resolution
  suppressMessages({
    files <- ebirdst_download_status("leafly",
                                     pattern = "_27km_",
                                     dry_run = TRUE)
  })
  expect_error(ebirdst_download_status("XXXX"))
  expect_error(ebirdst_download_status("Wod Thrush"))
  expect_error(ebirdst_download_status("aakspa1"))
})


# trends ----

test_that("ebirdst_download_trends()", {
  # download one species
  path1 <- ebirdst_download_trends("Pomatorhinus musicus",
                                   show_progress = FALSE)
  # download multiple species
  path2 <- ebirdst_download_trends(c("taibap1", "Taiwan Barbet"),
                                   show_progress = FALSE)

  paths <- file.path(c(path1, path2), "trends")
  expect_true(all(dir.exists(paths)))

  n_trends_files <- sapply(paths, function(x) length(list.files(x)))
  n_trends_files <- unname(n_trends_files)
  expect_equal(n_trends_files, rep_len(3L, length.out = length(paths)))
})

test_that("ebirdst_download_trends() missing species", {
  # missing from s&t
  expect_error(ebirdst_download_trends("XXXX"))
  expect_error(ebirdst_download_trends(c("Wod Thrush", "aakspa1")))

  # missing from trends
  expect_error(ebirdst_download_trends(c("amewoo", "Vultur gryphus")))
})


test_that("get_species_path()", {
  expect_true(dir.exists(get_species_path("yebsap-example")))
  expect_error(get_species_path("XXXXX"))
  expect_error(get_species_path("XXXXX", check_downloaded = FALSE))
  expect_false(dir.exists(get_species_path("Yellow Warbler",
                                           check_downloaded = FALSE)))
})
