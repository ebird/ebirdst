skip_on_cran()

test_that("set_ebirdst_access_key() stores the key in ~/.Renviron", {
  tmp_home <- withr::local_tempdir()
  withr::local_envvar(HOME = tmp_home, EBIRDST_KEY = NA)

  renv_path <- suppressMessages(set_ebirdst_access_key("XXXXXX"))
  expect_equal(
    normalizePath(renv_path),
    normalizePath(file.path(tmp_home, ".Renviron"))
  )
  expect_true(any(grepl("EBIRDST_KEY='XXXXXX'", readLines(renv_path))))
  expect_equal(Sys.getenv("EBIRDST_KEY"), "XXXXXX")
})


test_that("set_ebirdst_access_key() errors if the key is already set", {
  tmp_home <- withr::local_tempdir()
  withr::local_envvar(HOME = tmp_home, EBIRDST_KEY = NA)

  suppressMessages(set_ebirdst_access_key("XXXXXX"))
  expect_error(set_ebirdst_access_key("YYYYYY"), "already set")
  suppressMessages(set_ebirdst_access_key("YYYYYY", overwrite = TRUE))
  expect_equal(Sys.getenv("EBIRDST_KEY"), "YYYYYY")
})


test_that("set_ebirdst_access_key() warns when a project-level .Renviron exists", {
  tmp_home <- withr::local_tempdir()
  tmp_project <- withr::local_tempdir()
  withr::local_envvar(HOME = tmp_home, EBIRDST_KEY = NA)
  withr::local_dir(tmp_project)
  file.create(".Renviron")

  expect_warning(
    suppressMessages(set_ebirdst_access_key("XXXXXX")),
    "project-level \\.Renviron"
  )
})


test_that("set_ebirdst_access_key() doesn't warn when the working directory is the home directory", {
  tmp_home <- withr::local_tempdir()
  withr::local_envvar(HOME = tmp_home, EBIRDST_KEY = NA)
  withr::local_dir(tmp_home)

  expect_no_warning(suppressMessages(set_ebirdst_access_key("XXXXXX")))
})
