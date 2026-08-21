skip_on_cran()

# version years are taken from ebirdst_version() rather than hardcoded so these
# fixtures don't silently stop testing what they're meant to at the next data
# release
status_year <- ebirdst_version()[["status_version_year"]]
trends_year <- ebirdst_version()[["trends_version_year"]]

# these keys are independent of any real download, so no network or access
# key is required to test them
keys <- file.path(
  status_year,
  "woothr",
  c(
    "config.json",
    sprintf("weekly/woothr_abundance_median_3km_%s.tif", status_year),
    sprintf(
      "weekly/woothr_proportion-population_median_3km_%s.tif",
      status_year
    ),
    sprintf("weekly/woothr_occurrence_median_3km_%s.tif", status_year),
    sprintf("weekly/woothr_count_median_3km_%s.tif", status_year),
    sprintf("ranges/woothr_range_smooth_27km_%s.gpkg", status_year),
    "regional_stats.csv",
    sprintf("pis/woothr_pi_occurrence_gsw_c2_pland_27km_%s.tif", status_year),
    sprintf("ppms/woothr_ppm_binary-f1_mean_27km_%s.tif", status_year)
  )
)

test_that("status_key() and trends_key()", {
  expect_equal(status_key("woothr"), paste0(status_year, "/woothr"))
  expect_equal(
    status_key("woothr", "config.json"),
    paste0(status_year, "/woothr/config.json")
  )
  expect_equal(
    status_key("woothr", "weekly", c("a.tif", "b.tif")),
    paste0(status_year, "/woothr/weekly/", c("a.tif", "b.tif"))
  )
  expect_equal(
    trends_key("woothr", "trends", "a.parquet"),
    paste0(trends_year, "/woothr/trends/a.parquet")
  )
})


test_that("list_object_keys() reads the bundled list for the example data", {
  status_keys <- list_object_keys("yebsap-example", dataset = "status")
  expect_type(status_keys, "character")
  expect_true(length(status_keys) > 0)
  expect_true(all(grepl("yebsap-example", status_keys)))
  expect_true(any(grepl("config.json$", status_keys)))
  # the web_download folder isn't part of the data products
  expect_false(any(grepl("web_download", status_keys)))

  trends_keys <- list_object_keys("yebsap-example", dataset = "trends")
  expect_true(any(grepl("/trends/", trends_keys)))

  expect_error(list_object_keys(c("woothr", "yebsap")))
  expect_error(list_object_keys("yebsap-example", dataset = "XXXX"))
})


test_that("object_key_url() builds example and API urls", {
  example_key <- status_key("yebsap-example", "config.json")
  expect_match(object_key_url(example_key), "^https://raw.githubusercontent.com/")
  expect_match(object_key_url(example_key), paste0(example_key, "$"))

  # building an api url requires an access key, which the example data doesn't
  skip_if(Sys.getenv("EBIRDST_KEY") == "", "no access key available")
  api_url <- object_key_url(status_key("woothr", "config.json"))
  expect_match(api_url, "^https://st-download.ebird.org/v1/fetch\\?objKey=")
  expect_match(api_url, "&key=")

  # a mix of example and non-example keys keeps each url with its own key
  urls <- object_key_url(c(example_key, status_key("woothr", "config.json")))
  expect_length(urls, 2)
  expect_match(urls[1], "raw.githubusercontent.com")
  expect_match(urls[2], "st-download.ebird.org")
})


test_that("resolve_species()", {
  expect_equal(resolve_species("woothr"), "woothr")
  expect_equal(resolve_species("Wood Thrush"), "woothr")
  expect_error(resolve_species("XXXX"), "were not modeled")
  # all unrecognized species are listed, not just the first
  expect_error(resolve_species(c("woothr", "XXXX", "YYYY")), "XXXX, YYYY")
})


test_that("select_status_keys() default and flag-based selection", {
  default <- select_status_keys(keys)
  expect_true(all(grepl("config.json$|_abundance_|_proportion-population_", default)))
  expect_false(any(grepl("_occurrence_|_count_|ranges|regional_stats|pis|ppms", default)))

  expect_true(any(grepl("_occurrence_", select_status_keys(keys, download_occurrence = TRUE))))
  expect_true(any(grepl("_count_", select_status_keys(keys, download_count = TRUE))))
  expect_true(any(grepl("/ranges/", select_status_keys(keys, download_ranges = TRUE))))
  expect_true(any(grepl(
    "regional_stats.csv",
    select_status_keys(keys, download_regional = TRUE)
  )))
  expect_true(any(grepl("/pis/", select_status_keys(keys, download_pis = TRUE))))
  expect_true(any(grepl("/ppms/", select_status_keys(keys, download_ppms = TRUE))))

  # config is always retained, even when no flags are set
  none <- select_status_keys(keys, download_abundance = FALSE)
  expect_equal(none, keys[grepl("config.json$", keys)])
})


test_that("select_status_keys() download_all", {
  all_keys <- select_status_keys(keys, download_all = TRUE)
  expect_equal(sort(all_keys), sort(keys))
})


test_that("select_status_keys() pattern argument", {
  filtered <- select_status_keys(keys, download_all = TRUE, pattern = "27km")
  expect_true(all(grepl("config.json$|27km", filtered)))
  expect_true(any(grepl("27km", filtered)))

  expect_error(
    select_status_keys(keys, download_all = TRUE, pattern = "zzzzz"),
    "No files matched pattern"
  )
})


test_that("ensure_data_dir()", {
  tmp <- withr::local_tempdir()
  nested <- file.path(tmp, "a", "b")
  ensure_data_dir(nested)
  expect_true(dir.exists(nested))

  # a path nested under a file rather than a directory can't be created
  blocker <- file.path(tmp, "blocker")
  file.create(blocker)
  expect_error(ensure_data_dir(file.path(blocker, "nested")))
})


test_that("fetch_data() downloads example data and returns normalized paths", {
  tmp <- withr::local_tempdir()
  key <- status_key("yebsap-example", "config.json")
  local_path <- suppressMessages(fetch_data(key, path = tmp, show_progress = FALSE))
  expect_equal(local_path, normalizePath(file.path(tmp, key)))
  expect_true(file.exists(local_path))
})


test_that("fetch_data() reports existing files only when requested", {
  tmp <- withr::local_tempdir()
  key <- status_key("yebsap-example", "config.json")
  suppressMessages(fetch_data(key, path = tmp, show_progress = FALSE))

  expect_silent(fetch_data(key, path = tmp, show_progress = FALSE))
  expect_message(
    fetch_data(key, path = tmp, show_progress = FALSE, report_existing = TRUE),
    "already exists"
  )
})


test_that("fetch_data() re-downloads on force = TRUE", {
  tmp <- withr::local_tempdir()
  key <- status_key("yebsap-example", "config.json")
  suppressMessages(fetch_data(key, path = tmp, show_progress = FALSE))
  expect_message(
    fetch_data(key, path = tmp, force = TRUE, show_progress = TRUE),
    "Downloading file"
  )
})


test_that("fetch_data() errors, with hint, for a nonexistent key", {
  tmp <- withr::local_tempdir()
  bad_key <- status_key("yebsap-example", "does-not-exist.tif")
  expect_error(
    fetch_data(bad_key, path = tmp, show_progress = FALSE, hint = "some hint"),
    "some hint"
  )
})


test_that("http_url() only downgrades https urls", {
  expect_equal(
    http_url("https://st-download.ebird.org/v1"),
    "http://st-download.ebird.org/v1"
  )
  expect_equal(
    http_url("http://st-download.ebird.org/v1"),
    "http://st-download.ebird.org/v1"
  )
  expect_equal(
    http_url(c("https://a.org/x", "http://b.org/y")),
    c("http://a.org/x", "http://b.org/y")
  )
})


test_that("redact_access_key() removes the key from messages", {
  withr::local_envvar(EBIRDST_KEY = "")

  expect_equal(
    redact_access_key(
      "cannot open URL 'https://st-download.ebird.org/v1/fetch?objKey=2023/categr1/config.json&key=mrt5f8dhtp3g': HTTP status was '500 Internal Server Error'"
    ),
    "cannot open URL 'https://st-download.ebird.org/v1/fetch?objKey=2023/categr1/config.json&key=<redacted>': HTTP status was '500 Internal Server Error'"
  )
  # the listing url takes the key as the first query parameter
  expect_equal(
    redact_access_key("URL 'https://x/list-obj/2023/categr1?key=abc123' failed"),
    "URL 'https://x/list-obj/2023/categr1?key=<redacted>' failed"
  )
  # objKey isn't mistaken for the key itself, and is left intact
  expect_equal(
    redact_access_key("?objKey=2023/woothr/config.json"),
    "?objKey=2023/woothr/config.json"
  )
  expect_equal(redact_access_key("nothing to redact"), "nothing to redact")
  expect_equal(redact_access_key(""), "")

  # the key value is also matched literally, wherever it appears
  withr::local_envvar(EBIRDST_KEY = "mrt5f8dhtp3g")
  expect_equal(
    redact_access_key("the key mrt5f8dhtp3g is invalid"),
    "the key <redacted> is invalid"
  )
})


test_that("try_url() reports a redacted reason for the failure", {
  withr::local_envvar(EBIRDST_KEY = "")

  attempt <- try_url(
    stop("cannot open URL 'https://x/fetch?objKey=y&key=secret123': HTTP status was '500 Internal Server Error'")
  )
  expect_false(grepl("secret123", attempt$reason, fixed = TRUE))
  expect_match(attempt$reason, "key=<redacted>", fixed = TRUE)
  expect_match(attempt$reason, "500 Internal Server Error")

  # nothing to report when the request succeeded
  expect_equal(try_url(42L)$reason, "")
})


test_that("try_url() distinguishes http status errors from failed connections", {
  ok <- try_url(42L)
  expect_equal(ok$value, 42L)
  expect_false(ok$http_status)

  # the server responded, so the object simply isn't available
  not_found <- try_url(
    stop("cannot open URL 'x': HTTP status was '404 Not Found'")
  )
  expect_null(not_found$value)
  expect_true(not_found$http_status)

  # a connection level failure also reports "status was", but without "HTTP "
  unreachable <- try_url(stop("URL 'x': status was 'Couldn't resolve host name'"))
  expect_null(unreachable$value)
  expect_false(unreachable$http_status)

  # warnings are captured for classification, then muffled
  expect_silent({
    warned <- try_url({
      warning("cannot open URL 'x': HTTP status was '404 Not Found'")
      NULL
    })
  })
  expect_true(warned$http_status)
})


test_that("check_object_keys() rejects keys that escape the data directory", {
  expect_equal(check_object_keys(keys), keys)
  expect_true(withVisible(check_object_keys(keys))$visible == FALSE)

  expect_error(check_object_keys("2023/woothr/../../evil.tif"), "invalid names")
  expect_error(check_object_keys("../evil.tif"), "invalid names")
  expect_error(check_object_keys("2023\\woothr\\..\\evil.tif"), "invalid names")
  expect_error(check_object_keys("/etc/passwd"), "invalid names")
  expect_error(check_object_keys("\\\\server\\share"), "invalid names")
  expect_error(check_object_keys("C:/windows/evil.tif"), "invalid names")
  expect_error(check_object_keys(""), "invalid names")

  # a key may still contain dots, just not a .. path segment
  expect_equal(
    check_object_keys("2023/woothr/a..b.tif"),
    "2023/woothr/a..b.tif"
  )

  expect_error(check_object_keys(character(0)))
  expect_error(check_object_keys(NA_character_))
  expect_error(check_object_keys(1))
})


test_that("fetch_data() rejects keys that escape the data directory", {
  tmp <- withr::local_tempdir()
  expect_error(
    fetch_data("2023/woothr/../../evil.tif", path = tmp, show_progress = FALSE),
    "invalid names"
  )
})


test_that("partial_download_path() and is_partial_download()", {
  expect_equal(partial_download_path("2023/woothr/a.tif"), "2023/woothr/a.tif.part")
  expect_true(is_partial_download("2023/woothr/a.tif.part"))
  expect_false(is_partial_download("2023/woothr/a.tif"))
  expect_equal(is_partial_download(c("a.tif", "a.tif.part")), c(FALSE, TRUE))
})


test_that("download_files() copies a completed download into place", {
  tmp_dir <- withr::local_tempdir()
  key <- status_key("yebsap-example", "config.json")
  dest <- file.path(tmp_dir, "config.json")

  result <- download_files(
    object_key_url(key),
    dest,
    key,
    show_progress = FALSE
  )
  expect_true(result$success)
  expect_false(result$not_found)
  expect_true(file.exists(dest))
  # the temporary file used during the download isn't left behind
  expect_false(file.exists(paste0(dest, ".part")))
})

test_that("download_files() overwrites an existing destination file", {
  # a forced re-download targets a destination that already exists; renaming
  # onto an existing file fails on windows, so this must succeed via a copy
  tmp_dir <- withr::local_tempdir()
  key <- status_key("yebsap-example", "config.json")
  dest <- file.path(tmp_dir, "config.json")
  writeLines("stale content", dest)

  result <- download_files(
    object_key_url(key),
    dest,
    key,
    show_progress = FALSE
  )
  expect_true(result$success)
  expect_false(identical(readLines(dest), "stale content"))
})


test_that("download_files() leaves nothing behind when a download fails", {
  tmp_dir <- withr::local_tempdir()
  dest <- file.path(tmp_dir, "nothing.tif")

  # nothing is listening, so the connection fails before any bytes arrive
  result <- download_files(
    "http://127.0.0.1:1/nothing.tif",
    dest,
    "nothing.tif",
    show_progress = FALSE
  )
  expect_false(result$success)
  expect_false(result$not_found)
  expect_length(list.files(tmp_dir, all.files = TRUE, no.. = TRUE), 0)
})


test_that("download_files() flags data the server says doesn't exist", {
  tmp_dir <- withr::local_tempdir()
  key <- status_key("yebsap-example", "does-not-exist.tif")
  dest <- file.path(tmp_dir, "does-not-exist.tif")

  result <- download_files(
    object_key_url(key),
    dest,
    key,
    show_progress = FALSE
  )
  expect_false(result$success)
  expect_true(result$not_found)
  expect_length(list.files(tmp_dir, all.files = TRUE, no.. = TRUE), 0)
})


test_that("download_files() doesn't fall back to http after a partial transfer", {
  original <- ebirdst_env$api_base_url
  withr::defer(ebirdst_env$api_base_url <- original)
  tmp_dir <- withr::local_tempdir()
  dest <- file.path(tmp_dir, "partial.tif")

  # a truncated transfer errors but leaves a partial file behind, which means the
  # server responded, so the access key must not be re-sent over http
  n_attempts <- 0L
  local_mocked_bindings(
    try_url = function(expr) {
      n_attempts <<- n_attempts + 1L
      writeBin(raw(10), paste0(dest, ".part"))
      return(list(
        value = NULL,
        http_status = FALSE,
        reason = "Transferred a partial file"
      ))
    }
  )
  result <- download_files(
    paste0(api_base_url(), "/fetch?objKey=x&key=bogus"),
    dest,
    "partial.tif",
    show_progress = FALSE
  )
  expect_equal(n_attempts, 1L)
  expect_false(result$success)
  expect_equal(api_base_url(), original)
  # the partial file is discarded rather than passed off as a complete download
  expect_length(list.files(tmp_dir, all.files = TRUE, no.. = TRUE), 0)
})

test_that("download_files() retries a transient connection failure", {
  tmp_dir <- withr::local_tempdir()
  dest <- file.path(tmp_dir, "retry.tif")

  # the first attempt fails with nothing at all coming back, which is treated
  # as a transient blip and retried rather than given up on immediately
  n_attempts <- 0L
  local_mocked_bindings(
    try_url = function(expr) {
      n_attempts <<- n_attempts + 1L
      if (n_attempts < 2L) {
        return(list(value = NULL, http_status = FALSE, reason = "connection reset"))
      }
      writeBin(raw(10), paste0(dest, ".part"))
      return(list(value = 0L, http_status = FALSE, reason = ""))
    }
  )
  result <- download_files(
    "http://127.0.0.1:1/retry.tif",
    dest,
    "retry.tif",
    show_progress = FALSE
  )
  expect_equal(n_attempts, 2L)
  expect_true(result$success)
  expect_true(file.exists(dest))
})

test_that("download_files() gives up after repeated transient failures", {
  tmp_dir <- withr::local_tempdir()
  dest <- file.path(tmp_dir, "retry.tif")

  n_attempts <- 0L
  local_mocked_bindings(
    try_url = function(expr) {
      n_attempts <<- n_attempts + 1L
      return(list(value = NULL, http_status = FALSE, reason = "connection reset"))
    }
  )
  result <- download_files(
    "http://127.0.0.1:1/retry.tif",
    dest,
    "retry.tif",
    show_progress = FALSE
  )
  # the initial attempt plus 2 retries, then no more
  expect_equal(n_attempts, 3L)
  expect_false(result$success)
})


test_that("fetch_data() distinguishes a failed download from missing data", {
  tmp <- withr::local_tempdir()
  key <- status_key("yebsap-example", "config.json")
  local_mocked_bindings(
    object_key_url = function(keys) {
      return(rep("http://127.0.0.1:1/x", length(keys)))
    }
  )
  expect_error(
    fetch_data(key, path = tmp, show_progress = FALSE, hint = "some hint"),
    "failed to download"
  )
})


test_that("fetch_data() never reports the access key in a download error", {
  tmp <- withr::local_tempdir()
  withr::local_envvar(EBIRDST_KEY = "secret123")
  key <- status_key("yebsap-example", "config.json")

  # nothing is listening on this port, so the download fails with a message
  # naming the url, which carries the access key
  local_mocked_bindings(
    object_key_url = function(keys) {
      return(sprintf("http://127.0.0.1:1/fetch?objKey=%s&key=secret123", keys))
    }
  )
  err <- expect_error(fetch_data(key, path = tmp, show_progress = FALSE))
  msg <- conditionMessage(err)
  expect_false(grepl("secret123", msg, fixed = TRUE))
  expect_match(msg, "key=<redacted>", fixed = TRUE)
})


test_that("fetch_data() keeps the existing file when a forced re-download fails", {
  tmp <- withr::local_tempdir()
  key <- status_key("yebsap-example", "config.json")
  local_path <- suppressMessages(fetch_data(
    key,
    path = tmp,
    show_progress = FALSE
  ))
  original <- readLines(local_path, warn = FALSE)

  local_mocked_bindings(
    object_key_url = function(keys) {
      return(rep("http://127.0.0.1:1/x", length(keys)))
    }
  )
  expect_error(
    fetch_data(key, path = tmp, force = TRUE, show_progress = FALSE),
    "failed to download"
  )
  expect_true(file.exists(local_path))
  expect_equal(readLines(local_path, warn = FALSE), original)
})


test_that("download_files() stays on https when the server sends a response", {
  # restore the session cached base url in case the api can't be reached
  original <- ebirdst_env$api_base_url
  withr::defer(ebirdst_env$api_base_url <- original)

  # a bad key gets an http status error rather than a failed connection, so
  # the access key must not be re-sent over an unencrypted connection
  src <- paste0(
    api_base_url(),
    "/fetch?objKey=",
    status_key("woothr", "does-not-exist.tif"),
    "&key=bogus"
  )
  suppressMessages(download_files(
    src,
    withr::local_tempfile(),
    "does-not-exist.tif",
    show_progress = FALSE
  ))
  expect_equal(api_base_url(), original)
})
