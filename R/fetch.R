# this file contains the internal machinery shared by the ebirdst_download_*()
# functions in download.R and the on-demand downloads performed by the
# load_*() functions in load.R. the local path for a downloaded file is always
# its object key (e.g. "2023/woothr/config.json") appended to the data
# directory, and the API can fetch a single object directly by key, so
# fetch_data() is the one function anything in the package needs to call to
# make sure a set of files exist locally

# internal ----

# session-cached API base url; some VPNs block https to the download API, so
# a fallback to http is cached here once a request over http is known to have
# succeeded, to avoid re-probing on every request. the access key is passed in
# the query string, so the downgrade is only ever cached on success and only
# for a connection-level failure (see try_url())
ebirdst_env <- new.env(parent = emptyenv())
ebirdst_env$api_base_url <- "https://st-download.ebird.org/v1"

api_base_url <- function() {
  return(ebirdst_env$api_base_url)
}

http_url <- function(url) {
  return(sub("^https://", "http://", url))
}

use_http_fallback <- function() {
  ebirdst_env$api_base_url <- http_url(ebirdst_env$api_base_url)
  return(invisible(ebirdst_env$api_base_url))
}


# the access key is passed to the API in the query string of the request url, and
# both download.file() and read_json() name that url in the conditions they
# raise. those messages get pasted into bug reports and emails, so the key has to
# be stripped out of anything the package passes on to the user. the query
# parameter is matched rather than the key itself so that this works even when no
# key is set locally, e.g. for the example data; the key is then also matched
# literally in case it ever appears somewhere the query string pattern doesn't
redact_access_key <- function(x) {
  redacted <- stringr::str_replace_all(
    x,
    "([?&])key=[^&'\"\\s]*",
    "\\1key=<redacted>"
  )

  # Sys.getenv() is used directly because get_ebirdst_access_key() errors when
  # no key is set, and redacting must never itself be a point of failure
  access_key <- Sys.getenv("EBIRDST_KEY")
  if (nzchar(access_key)) {
    redacted <- stringr::str_replace_all(
      redacted,
      stringr::fixed(access_key),
      "<redacted>"
    )
  }

  return(redacted)
}


# attempt to access a url, returning the value of `expr` (NULL on failure)
# alongside a flag indicating whether the failure was an http status error.
# an http status error means the server was reached and responded, so the
# object simply isn't available; any other failure (dns, tls, proxy, timeout)
# is a connection-level problem and is the only case where retrying over http
# could help. this distinction matters because the access key travels in the
# query string, so http must never be probed for a request that already got a
# response over https
try_url <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    tryCatch(
      expr,
      error = function(e) {
        messages <<- c(messages, conditionMessage(e))
        return(NULL)
      }
    ),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  # note that a connection-level failure reports "status was '<reason>'"
  # without the "HTTP " prefix, so this matches responses only
  http_status <- any(stringr::str_detect(messages, "HTTP status"))

  # the underlying message is the only clue as to why a request failed, so it's
  # worth reporting, but only ever redacted
  reason <- redact_access_key(paste(unique(messages), collapse = "; "))

  return(list(value = value, http_status = http_status, reason = reason))
}


# object keys are appended to the data directory to give the local path of a
# downloaded file, and most of them come from the API listing rather than from
# the user, so check that none of them could write outside the data directory
# before using one as a path
check_object_keys <- function(keys) {
  stopifnot(is.character(keys), length(keys) >= 1, !anyNA(keys))

  segments <- strsplit(keys, "[/\\\\]")
  traversal <- vapply(segments, function(x) any(x == ".."), logical(1))
  absolute <- stringr::str_detect(keys, "^([/\\\\]|[A-Za-z]:)")
  invalid <- keys[traversal | absolute | keys == ""]
  if (length(invalid) > 0) {
    stop(
      "The following data files have invalid names:\n  ",
      paste(invalid, collapse = "\n  ")
    )
  }

  return(invisible(keys))
}


# downloads are written to a temporary file with this suffix alongside their
# destination and only moved into place once complete, so the suffix is defined
# here rather than inline in download_files(): ebirdst_data_inventory() needs it
# to recognize and ignore a partial download left behind by a session that was
# killed mid-transfer
partial_suffix <- ".part"

partial_download_path <- function(path) {
  return(paste0(path, partial_suffix))
}

is_partial_download <- function(path) {
  return(stringr::str_ends(path, stringr::fixed(partial_suffix)))
}


# resolve a species name/code to its eBird species code, raising an error for
# any species not modeled by status and trends; all internal callers that need
# a valid code should use this so the error message is consistent
resolve_species <- function(species) {
  species_code <- get_species(species)
  if (anyNA(species_code)) {
    stop(
      "The following species were not modeled by eBird Status and Trends. ",
      "Consult ebirdst_runs for a complete list of available species.\n  ",
      paste(species[is.na(species_code)], collapse = ", ")
    )
  }
  return(species_code)
}


# create the data directory if it doesn't already exist
ensure_data_dir <- function(path) {
  if (!dir.exists(path)) {
    created <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(created)) {
      stop("Unable to create data directory: ", path)
    }
  }
  return(invisible(path))
}


# build object keys within the status or trends data package for a species;
# vectorizes over the last argument, e.g. status_key("woothr", "weekly", files)
status_key <- function(species_code, ...) {
  version_year <- ebirdst_version()[["status_version_year"]]
  return(paste(version_year, species_code, ..., sep = "/"))
}

trends_key <- function(species_code, ...) {
  version_year <- ebirdst_version()[["trends_version_year"]]
  return(paste(version_year, species_code, ..., sep = "/"))
}


# request the object listing for a species from a given API base url; returns
# the result of try_url(), so the caller can tell a missing listing apart from
# an unreachable server
read_object_list <- function(version_year, species_code, base_url) {
  key <- get_ebirdst_access_key()
  list_obj_url <- stringr::str_glue(
    "{base_url}/list-obj/{version_year}/",
    "{species_code}?key={key}"
  )
  return(try_url(jsonlite::read_json(list_obj_url, simplifyVector = TRUE)))
}


# list all object keys available for a species, for callers that don't
# already know the exact key(s) they want: flag/pattern-based selection in
# ebirdst_download_status()/ebirdst_download_trends(), and PI availability in
# list_available_pis()
list_object_keys <- function(species_code, dataset = c("status", "trends")) {
  stopifnot(
    is.character(species_code),
    length(species_code) == 1,
    !is.na(species_code)
  )
  dataset <- match.arg(dataset)

  version_year <- ebirdst_version()[[paste0(dataset, "_version_year")]]
  is_example <- (species_code == "yebsap-example")

  if (is_example) {
    fl <- system.file(
      "extdata",
      paste0("example-data_file-list_", dataset, ".txt"),
      package = "ebirdst"
    )
    keys <- readLines(fl)
  } else {
    attempt <- read_object_list(version_year, species_code, api_base_url())
    keys <- attempt$value

    # some vpns block https to the download api, so retry over http, but only
    # if https failed to connect at all rather than returning a response
    retry_http <- is.null(keys) &&
      !attempt$http_status &&
      startsWith(api_base_url(), "https://")
    if (retry_http) {
      attempt <- read_object_list(
        version_year,
        species_code,
        http_url(api_base_url())
      )
      keys <- attempt$value
      # only cache the downgrade now that it's known to work
      if (!is.null(keys)) {
        use_http_fallback()
      }
    }

    if (is.null(keys)) {
      reason <- if (nzchar(attempt$reason)) {
        paste0("\nThe following error occurred:\n  ", attempt$reason)
      } else {
        ""
      }
      stop(
        "Cannot access Status and Trends data URL. Ensure that you have ",
        "a working internet connection and a valid API key for the ",
        "Status and Trends data. Note that the API keys expire after ",
        "6 months, so you may need to update your key. ",
        "Visit https://ebird.org/st/request",
        reason
      )
    }

    # remove web_download folder
    web_down <- stringr::str_detect(dirname(keys), pattern = "web_download")
    keys <- keys[!web_down]

    # remove additional species caused by bug in API, e.g. leafly will also
    # return leafly2
    only_target <- stringr::str_detect(
      keys,
      pattern = paste0("/", species_code, "/")
    )
    keys <- keys[only_target]
  }

  if (length(keys) == 0) {
    stop("No data found for species ", species_code)
  }

  return(keys)
}


# select which object keys should be downloaded based on the download_* flags
# and an optional filename pattern; the selection logic used by
# ebirdst_download_status()
select_status_keys <- function(
  keys,
  download_abundance = TRUE,
  download_occurrence = FALSE,
  download_count = FALSE,
  download_ranges = FALSE,
  download_regional = FALSE,
  download_pis = FALSE,
  download_ppms = FALSE,
  download_all = FALSE,
  pattern = NULL
) {
  # always download config file
  dl <- stringr::str_detect(keys, pattern = "config.json$")
  if (download_abundance || download_all) {
    # add abundance
    dl <- stringr::str_detect(keys, "\\_abundance\\_") | dl
    # add proportion of population
    dl <- stringr::str_detect(keys, "\\_proportion-population\\_") | dl
  }
  if (download_occurrence || download_all) {
    # add occurrence
    dl <- stringr::str_detect(keys, "\\_occurrence\\_") | dl
  }
  if (download_count || download_all) {
    # add count
    dl <- stringr::str_detect(keys, "\\_count\\_") | dl
  }
  if (download_ranges || download_all) {
    # add ranges
    dl <- stringr::str_detect(keys, "/ranges/") | dl
  }
  if (download_regional || download_all) {
    # add regional summary stats
    dl <- stringr::str_ends(keys, "regional_stats.csv") | dl
  }
  if (download_pis || download_all) {
    # add pis
    dl <- stringr::str_detect(keys, "/pis/") | dl
  }
  if (download_ppms || download_all) {
    # add ppms
    dl <- stringr::str_detect(keys, "/ppms/") | dl
  }
  keys <- keys[dl]

  # apply pattern
  if (!is.null(pattern)) {
    stopifnot(is.character(pattern), length(pattern) == 1, !is.na(pattern))
    pat_match <- stringr::str_detect(basename(keys), pattern = pattern)
    if (!any(pat_match)) {
      stop("No files matched pattern")
    }

    # always download config file
    is_config <- stringr::str_detect(basename(keys), pattern = "config.json$")
    keys <- keys[pat_match | is_config]
  }

  return(keys)
}


# build the source download url for a set of object keys
object_key_url <- function(keys) {
  is_example <- stringr::str_detect(keys, "yebsap-example")
  urls <- character(length(keys))

  if (any(is_example)) {
    example_url <- paste0(
      "https://raw.githubusercontent.com/",
      "ebird/ebirdst_example-data/main/",
      "example-data/"
    )
    urls[is_example] <- paste0(example_url, keys[is_example])
  }
  if (!all(is_example)) {
    key <- get_ebirdst_access_key()
    urls[!is_example] <- stringr::str_glue(
      "{api_base_url()}/fetch?objKey={keys[!is_example]}",
      "&key={key}"
    )
  }

  return(urls)
}


# ensure the local files for a set of object keys exist, downloading any that
# are missing (or all of them, if force = TRUE); returns the normalized local
# paths. every download in the package funnels through here. `hint` is
# appended to the error raised if a requested key can't be found, and
# `report_existing` controls whether "already downloaded" messages are shown
# (used by the ebirdst_download_*() functions, but not by on-demand loads,
# which should stay silent when the requested data is already cached)
fetch_data <- function(
  keys,
  path,
  force = FALSE,
  show_progress = interactive(),
  hint = NULL,
  report_existing = FALSE
) {
  check_object_keys(keys)
  ensure_data_dir(path)
  dest_paths <- file.path(path, keys)
  exists <- file.exists(dest_paths)

  if (!isTRUE(force) && all(exists)) {
    if (report_existing) {
      message("Data already exists, use force = TRUE to re-download.")
    }
    return(invisible(normalizePath(dest_paths)))
  }
  if (!isTRUE(force) && any(exists) && report_existing) {
    message(
      "Some files already exist, only downloading new files. ",
      "Use force = TRUE to re-download all files."
    )
  }

  to_fetch <- if (isTRUE(force)) keys else keys[!exists]
  fetch_dest <- file.path(path, to_fetch)

  # create necessary directories
  dirs <- unique(dirname(fetch_dest))
  for (d in dirs) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }

  result <- download_files(
    object_key_url(to_fetch),
    fetch_dest,
    to_fetch,
    show_progress = show_progress
  )

  # a download that failed for any reason other than the data not existing, e.g.
  # a dropped connection, isn't something the caller can fix by requesting
  # different data, so it gets its own error rather than the hint below. any
  # local copy of these files is left as it was
  failed <- !result$success & !result$not_found
  if (any(failed)) {
    detail <- ifelse(
      is.na(result$reason[failed]) | result$reason[failed] == "",
      to_fetch[failed],
      paste0(to_fetch[failed], ": ", result$reason[failed])
    )
    stop(
      "The following files failed to download:\n  ",
      paste(detail, collapse = "\n  "),
      "\nThis is usually a temporary problem, check your internet connection ",
      "and try again."
    )
  }

  missing <- keys[!file.exists(dest_paths)]
  if (length(missing) > 0) {
    msg <- paste0(
      "The requested data could not be found:\n  ",
      paste(missing, collapse = "\n  ")
    )
    if (!is.null(hint)) {
      stop(msg, "\n", hint)
    }
    stop(msg)
  }

  return(invisible(normalizePath(dest_paths)))
}


# download files from src urls to local destination paths, returning the outcome
# for each file: `success`, `not_found` for the files the server responded to
# with an http status error, meaning the data simply isn't there as opposed to
# the download failing for some other reason, and `reason`, the redacted message
# from the failed attempt (NA where the download succeeded). fetch_data() needs
# these to report a useful error. `keys` is used only to report progress
#
# each file is downloaded to a temporary file alongside its destination and only
# copied into place, overwriting any existing destination, once the transfer
# has completed, because download.file() leaves a partial file behind when a
# transfer is cut short part way, and deletes any existing destination file
# when it fails. file.copy(overwrite = TRUE) is used rather than file.rename()
# because rename() fails when the destination already exists on windows,
# which is exactly the case for a forced re-download. every temporary file is
# removed on any exit from this function, including an error or interrupt
#
# if https can't be reached at all, retry once over http in case it's being
# blocked (e.g. by a VPN), caching the fallback for the rest of the session only
# once it's known to work. a connection-level failure that persists after that
# is often just a transient blip (e.g. on a large batch download), so it gets a
# couple more retries with a short backoff before the file is given up on
download_files <- function(src, dest, keys, show_progress) {
  n_files <- length(src)
  old_timeout <- getOption("timeout")
  options(timeout = max(3000, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  tmp <- partial_download_path(dest)
  on.exit(unlink(tmp), add = TRUE)

  success <- rep(FALSE, n_files)
  not_found <- rep(FALSE, n_files)
  reason <- rep(NA_character_, n_files)

  for (i in seq_len(n_files)) {
    if (show_progress) {
      message(stringr::str_glue(
        "  Downloading file {i} of {n_files}: ",
        "{basename(keys[i])}"
      ))
    }
    attempt <- try_url(
      utils::download.file(src[i], tmp[i], quiet = TRUE, mode = "wb")
    )
    ok <- identical(attempt$value, 0L)

    # an http status, or a partial file, means the server responded, so only a
    # failure that left nothing at all behind is a connection-level problem
    # worth retrying over http
    retry_http <- !ok &&
      !attempt$http_status &&
      !file.exists(tmp[i]) &&
      stringr::str_starts(src[i], "https://st-download")
    if (retry_http) {
      attempt <- try_url(
        utils::download.file(
          http_url(src[i]),
          tmp[i],
          quiet = TRUE,
          mode = "wb"
        )
      )
      ok <- identical(attempt$value, 0L)
      # only cache the downgrade, and apply it to the files still to come,
      # once it's known to work
      if (ok) {
        use_http_fallback()
        is_api <- stringr::str_starts(src, "https://st-download")
        src[is_api] <- http_url(src[is_api])
      }
    }

    # a connection-level failure (nothing at all came back) is often transient,
    # so retry the same url a couple more times with a short backoff rather
    # than giving up on the file immediately
    retries <- 0L
    while (
      !ok && !attempt$http_status && !file.exists(tmp[i]) && retries < 2L
    ) {
      retries <- retries + 1L
      Sys.sleep(retries)
      attempt <- try_url(
        utils::download.file(src[i], tmp[i], quiet = TRUE, mode = "wb")
      )
      ok <- identical(attempt$value, 0L)
    }

    if (ok) {
      success[i] <- file.copy(tmp[i], dest[i], overwrite = TRUE)
      unlink(tmp[i])
    } else {
      not_found[i] <- attempt$http_status
      reason[i] <- attempt$reason
      unlink(tmp[i])
    }
  }

  return(invisible(list(
    success = success,
    not_found = not_found,
    reason = reason
  )))
}
