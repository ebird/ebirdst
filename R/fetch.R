# this file contains the internal machinery shared by the ebirdst_download_*()
# functions in download.R and the on-demand downloads performed by the
# load_*() functions in load.R. the local path for a downloaded file is always
# its object key (e.g. "2023/woothr/config.json") appended to the data
# directory, and the API can fetch a single object directly by key, so
# fetch_data() is the one function anything in the package needs to call to
# make sure a set of files exist locally

# internal ----

# session-cached API base url; some VPNs block https to the download API, so
# a fallback to http is cached here once discovered so it isn't re-probed on
# every request
ebirdst_env <- new.env(parent = emptyenv())
ebirdst_env$api_base_url <- "https://st-download.ebird.org/v1"

api_base_url <- function() {
  return(ebirdst_env$api_base_url)
}

use_http_fallback <- function() {
  ebirdst_env$api_base_url <- sub(
    "^https://",
    "http://",
    ebirdst_env$api_base_url
  )
  return(invisible(ebirdst_env$api_base_url))
}


# resolve a species name/code to its eBird species code; mirrors the
# validation in get_species_path() but doesn't require path to already exist
resolve_species <- function(species) {
  species_code <- get_species(species)
  if (anyNA(species_code)) {
    stop(
      paste(species[is.na(species_code)], collapse = ", "),
      " does not correspond to a valid Status and Trends species."
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
    key <- get_ebirdst_access_key()
    list_obj_url <- stringr::str_glue(
      "{api_base_url()}/list-obj/{version_year}/",
      "{species_code}?key={key}"
    )
    keys <- tryCatch(
      suppressWarnings({
        jsonlite::read_json(list_obj_url, simplifyVector = TRUE)
      }),
      error = function(e) NULL
    )
    if (is.null(keys)) {
      # try http instead in case of ssl issues on vpn
      use_http_fallback()
      list_obj_url <- stringr::str_glue(
        "{api_base_url()}/list-obj/{version_year}/",
        "{species_code}?key={key}"
      )
      keys <- tryCatch(
        suppressWarnings({
          jsonlite::read_json(list_obj_url, simplifyVector = TRUE)
        }),
        error = function(e) NULL
      )
      if (is.null(keys)) {
        stop(
          "Cannot access Status and Trends data URL. Ensure that you have ",
          "a working internet connection and a valid API key for the ",
          "Status and Trends data. Note that the API keys expire after ",
          "6 months, so you may need to update your key. ",
          "Visit https://ebird.org/st/request"
        )
      }
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

  download_files(
    object_key_url(to_fetch),
    fetch_dest,
    to_fetch,
    show_progress = show_progress
  )

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


# download files from src urls to local destination paths; on failure, retry
# once over http in case https is being blocked (e.g. by a VPN), caching the
# fallback for the rest of the session if it succeeds. a file that still
# can't be downloaded after the retry is simply left missing on disk, so
# fetch_data() can report it (with its caller-specific hint) rather than
# failing here with a generic message. `keys` is used only to report progress
download_files <- function(src, dest, keys, show_progress) {
  n_files <- length(src)
  old_timeout <- getOption("timeout")
  options(timeout = max(3000, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  for (i in seq_len(n_files)) {
    if (show_progress) {
      message(stringr::str_glue(
        "  Downloading file {i} of {n_files}: ",
        "{basename(keys[i])}"
      ))
    }
    dl_response <- tryCatch(
      suppressWarnings(
        utils::download.file(src[i], dest[i], quiet = TRUE, mode = "wb")
      ),
      error = function(e) 1L
    )
    if (
      dl_response != 0 && stringr::str_starts(src[i], "https://st-download")
    ) {
      use_http_fallback()
      src[i:n_files] <- sub("^https://", "http://", src[i:n_files])
      tryCatch(
        suppressWarnings(
          utils::download.file(src[i], dest[i], quiet = TRUE, mode = "wb")
        ),
        error = function(e) 1L
      )
    }
  }

  return(invisible(n_files))
}


# check that the geotiff driver is installed; required to load any of the
# raster data products
check_gtiff_support <- function() {
  drv <- terra::gdal(drivers = TRUE)
  drv <- drv$name[stringr::str_detect(drv$can, "read")]
  if (!"GTiff" %in% drv) {
    stop(
      "GDAL does not have GeoTIFF support. GeoTIFF support is required to ",
      "load Status and Trends raster data."
    )
  }
  return(invisible(TRUE))
}
